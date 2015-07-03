#include <stdio.h>
#include <ctype.h>
#include <stdlib.h>
#include <string.h>
#include "../libs/knhx.h"
#include <time.h>
#include <math.h>
#include <gsl/gsl_rng.h>
#include <gsl/gsl_randist.h>
#include <gsl/gsl_sf_gamma.h>

int propose_new_b(int *b,int num_branches, int * b_star);
int count_phenos( int set[], int parents[], int b[], int pheno[], int **counts);
double log_likelihood(int ** counts,int num_sec);
double log_lambda_likelihood(double lambda);
double log_b_likelihood(double lambda, int *b, double *b_length);
double propose_new_lambda(double old_lambda);
double log_likelihood_all(int **counts, int num_sec, double lambda, int *b, double *b_lengths);
double log_likelihood_lambda_constant(int **counts, int num_sec, int *b, double lambda,double *b_length);
double log_likelihood_b_constant(double lambda,int *b, double *b_length);



int number_branches;
int number_leaves;
int number_phenotypes;
double T, log_T;
gsl_rng *r;


int main(int argc, char *argv[]){
    knhx1_t *tree;
    int error;
    int i,j,k, mcmc_counter, changed_branch;
    kstring_t str;
    char *newick_str = NULL;
    int bytes_read;
    FILE * fp;
    char **names;
    int *temp_phenos;
    int *phenos;
    int *parents;
    double *branches_len ;
    int **leaves_under;
    int *n_leaves_under;
    int *b, *b_star;
    int *sections, num_sec;
    int **counts,temp_counter=0;
    double old_log_likelihood, lambda,temp;
    double proposal_log_likelihood, proposal_lambda;
    unsigned long int *b_counts;

    r = gsl_rng_alloc(gsl_rng_mt19937);
    gsl_rng_set(r,time(NULL));
    mcmc_counter = atoi(argv[3]);

    newick_str = get_newick_from_file(argv[1]);
    printf("%s\n",newick_str);

    tree = kn_parse(newick_str,&number_branches, &error);
    number_leaves = get_number_leaves(tree, number_branches);
    number_phenotypes = get_pheno_from_file(argv[2], &names, &temp_phenos, number_leaves);
    set_pheno_in_tree(tree, number_branches, number_leaves, names, temp_phenos);
    get_tree_data(tree, number_branches,number_leaves, &parents, &branches_len,&phenos);
    get_leaves_under(parents, &leaves_under, &n_leaves_under, number_leaves, number_branches);

    for(i = 0; i<number_branches; i++){ 
          printf("[%3d]\t%3d\t%3d\t%4g", i, parents[i], n_leaves_under[i], branches_len[i]);
          for (j = 0; j < n_leaves_under[i]; ++j)
          printf("\t%d", leaves_under[i][j]);
          putchar('\n');
    }

    if ((b_counts = calloc(number_branches, sizeof(unsigned long int))) == NULL){
        fprintf(stderr,"Out of memory. Could not allocate enough memory initializing b_counts.\n");
        exit(EXIT_FAILURE);
    }

    /* let's set memory for b and initialize it. */
    if ((b = calloc(number_branches, sizeof(int))) == NULL){
        fprintf(stderr,"Out of memory. Could not allocate enough memory initializing b.\n");
        exit(EXIT_FAILURE);
    }
    b[number_branches-1] = 1; /* all set to zero apart from the root which should always be 1 */

    if ((b_star = calloc(number_branches, sizeof(int))) == NULL){
        fprintf(stderr,"Out of memory. Could not allocate enough memory initializing b_star.\n");
        exit(EXIT_FAILURE);
    }
    b_star[number_branches-1] = 1; /* all set to zero apart from the root which should always be 1 */

    /* let's set memory for sections and intialize it according the the intial b. */
    if ((sections = calloc(number_leaves, sizeof(int))) == NULL){
        fprintf(stderr,"Out of memory. Could not allocate enough memory initializing sections.\n");
        exit(EXIT_FAILURE);
    }

    if ((counts = malloc(number_branches *sizeof(*counts))) == NULL){
        fprintf(stderr,"Out of memory. Could not allocate enough memory initializing counts.\n");
        exit(EXIT_FAILURE);
    }
    for(i = 0; i < number_branches; i++)
        if ((counts[i] = calloc(number_phenotypes, sizeof(int))) == NULL){
            fprintf(stderr,"Out of memory. Could not allocate enough memory initializing counts[i].\n");
            exit(EXIT_FAILURE);
        }
    num_sec = count_phenos(sections, parents, b, phenos, counts);
    /* let's calculate the total length of branches on the tree. */
    T = 0.0;
    for(i = 0; i < number_branches-1; i++)
        T += branches_len[i];
    log_T = log(T);
    lambda = 1 / T;
    /*    b[30] = 1;*/

    /*    temp = log_likelihood(counts,num_sec);
          printf("The value of log_likelihood output is:%e\n",temp);
          temp = log_lambda_likelihood(lambda);
          printf("The value of log_lambda_likelihood output is:%e\n",temp);
          temp = log_b_likelihood(lambda, b, branches_len);
          printf("The value of log_b_likelihood output is:%e\n",temp);
          temp = propose_new_lambda(lambda);
          printf("The value of proposed lambda is:%e\n",temp);*/
    old_log_likelihood = log_likelihood_all(counts,num_sec, lambda, b, branches_len);
    for(i = 0;i <num_sec;i++)
            for (j = 0; j<number_phenotypes;j++)
                counts[i][j] = 0;

/*    printf("The value of the likelihood is:%e\n",old_log_likelihood);*/
    /* let's do the mcmc now. */
    for(i = 0; i<mcmc_counter; i++){
        changed_branch = propose_new_b(b,number_branches, b_star);
        
/*        printf("proposed changed branch is: %d\n",changed_branch);
        for(j = 0; j<number_branches; j++)
            printf("%d\t%d\n",j,b_star[j]);
        printf("\n");*/
        
        num_sec = count_phenos(sections, parents, b_star, phenos, counts);

        /*for (j = 0; j<num_sec; j++)
            printf("counts are as follows:%d\t%d\n",counts[j][0],counts[j][1]);*/

        proposal_log_likelihood = log_likelihood_all(counts,num_sec,lambda,b_star,branches_len);
        if (gsl_rng_uniform(r) <(exp(proposal_log_likelihood - old_log_likelihood))){
            b[changed_branch] = b_star[changed_branch];
            old_log_likelihood = proposal_log_likelihood;
            temp_counter++;
        }
        for(j = 0;j <num_sec;j++)
            for (k = 0; k<number_phenotypes;k++)
                counts[j][k] = 0;
        for(j = 0;j <number_branches; j++)
            b_counts[j] += b[j];

        proposal_lambda = propose_new_lambda(lambda);
        if (proposal_lambda > 0.0){
            num_sec = count_phenos(sections, parents, b, phenos, counts);
            proposal_log_likelihood = log_likelihood_all(counts,num_sec,proposal_lambda,b,branches_len);
            if (gsl_rng_uniform(r) <(exp(proposal_log_likelihood - old_log_likelihood))){
                lambda = proposal_lambda;
                old_log_likelihood = proposal_log_likelihood;
            }
            for(j = 0;j <num_sec;j++)
                for (k = 0; k<number_phenotypes;k++)
                     counts[j][k] = 0;
        }
        /* we need to put something in here later */
        
    }

    printf("counter for acceptance is: %d\n",temp_counter);
    printf("posterior values are as follows.\n");
    for(i = 0; i<number_branches; i++)
        printf("[%d]\t%e\n",i,((double) b_counts[i])/mcmc_counter);

    /*for(i = 0; i<number_branches; i++)
      printf("b are: [%d]\t%d\n",i,b[i]);
      for(i = 0; i<number_leaves; i++)
      printf("phenos are: [%d]\t%d\n",i,phenos[i]);
      for(i = 0; i<number_branches;i++)
      printf("parents are: [%d]\t%d\n",i,parents[i]);
      for(i = 0; i<number_leaves;i++)
      printf("sections are: [%d]\t%d\n",i,sections[i]);



      printf("There are %d sections on this tree given b.\n",num_sec);
      for (i = 0; i<num_sec; i++)
      printf("counts are as follows:%d\t%d\n",counts[i][0],counts[i][1]);
      */

    /*log_likelihood_old = log_likelihood_all(counts, num_sec, lambda, b,branches_len);
    */   



    /* Now I have the parents, branchs_lens, phenos, leaves_under and n_leaves_under, I will start doing the MCMC.*/


    /* my debugging code.*/

   /*    for(i = 0; i<number_leaves; i++){
          printf("%s\t%d\n",names[i],temp_phenos[i]);
          }

          printf("number of nodes on the tree is %d\n",number_branches);
          for(i = 0; i<number_branches; i++){
          knhx1_t *p = tree + i;
          printf("[%3d] %10s\t%3d\t%3d\t%4g", i, p->name, p->parent, p->n, p->d);
          for (j = 0; j < p->n; ++j)
          printf("\t%d", p->child[j]);
          putchar('\n');
          }
          printf("\n");

          for(i = 0; i<number_branches; i++){ 
          printf("[%3d]\t%3d\t%3d\t%4g", i, parents[i], n_leaves_under[i], branchs_len[i]);
          for (j = 0; j < n_leaves_under[i]; ++j)
          printf("\t%d", leaves_under[i][j]);
          putchar('\n');
          }

          for (i = 0; i<number_leaves; i++){
          printf("name is: %d, pheno is: %d\n",i,phenos[i]);
          }


          for(i = 0; i< number_branches; i++){
          if(isleaf(tree+i))
          printf("node %d is a leaf.\n",(tree+i)->index);
          }
          printf("total number of leaves on this tree is %d.\n",get_number_leaves(tree,number_branches));
          */

    return 0;
}

/**************************************
 * given par and b and pheno arrays, we determine the number of sections on the tree and then
 * count the number of each phenotype in each section. Here we assume binary phenotypes.
 * We just simply start at the bottom of the tree (leaves) and go up until we reach a branch with a change point on it. We then set the leaves cluster number.
 * In addition we use an array of pointers that point to array of counts for each section. If a branch has a change point on it and it is reached by a leaf then 
 * we set the pointer to point counts and modify the counts as appropriate.
 */
int count_phenos(int set[], int par[], int b[], int pheno[], int **counts)
{
    int i, j, p;
    /*int set[N] = {0};*/
    int **code;
    /*    printf("number_branches: %d and number_leaves: %d.\n",number_branches,number_leaves);*/
    if ((code = malloc(number_branches *sizeof(int *))) == NULL){
        fprintf(stderr,"Out of memory. Could not assign memory when setting code.");
        exit(EXIT_FAILURE);
    }
    for (i = 0; i<number_branches; i++)
        code[i] = NULL;

    j = 0;
    for (i = 0; i<number_leaves; i++)
    {
        p = i;
        while(!b[p])
            p = par[p];
        /*printf("At leaf %d and parent is %d.\n",i,p);*/

        set[i] = p;
        if (!code[p])
        {
            code[p] = &counts[j][0];
            j++;
        }
        code[p][pheno[i]]++;
    }
    return j;
}

/* an alternative for the count can be that when we add a change point, only the leaves below that node will be affected.
 * We go throught those nodes (pre-stored). For each one we go back on the tree, if we get to something else before them, then nothing changes.
 * if we get to the change point, then we take one off the right  pheno for the section that the leaf was on before and add one to the new section.
 * If none of the leaves get to the change point, then there is no change at all. keep a list of sections where the counts are changing. at the end if the counts for the section
 * are all zero, then set it to NULL and reduce the counter.
 *
 * When we remove a change point, again it either has an affect or not on the sections. Again we only look at the leaves under the branch with change points.
 * for each leaf we go back. Only go back up if the leaf was part of the section for the affected branch. if it is not part of the section, nothing has changed. if it is then go up, if we reach a branch with change point then take one of the right pheno for the section and add one to the pheno for the new section. If section is new then point to it. make sure that the counts for the old section is zero and then set the pointer to NULL.
 */
int propose_new_b(int *b,int num_branches, int *b_star)
{
    int j;
    int i = gsl_rng_uniform_int(r,num_branches-1);
/*    i = 23;*/
    for(j = 0; j < num_branches; j++)
        b_star[j] = b[j];
    b_star[i] = ! b_star[i];
    return i;
}

/* function to calculate the log of the likelihood and prior.
*/
double log_likelihood(int ** counts,int num_sec)
{
    double res = 0.0;
    int i,j, sec_total;
    for(i = 0; i<num_sec; i++){
        sec_total = 0;
        for (j = 0; j<number_phenotypes; j++){
            res += gsl_sf_lngamma(1.0 + counts[i][j]);
            sec_total += counts[i][j];
        }
        res -= gsl_sf_lngamma(number_phenotypes + sec_total);
    }
    res += (gsl_sf_lngamma(number_phenotypes) * num_sec);
    return res;
}


/* function to calculate the log likelihood of the lambda parameter. */

double log_lambda_likelihood(double lambda)
{
    double res;
    res = (-T * lambda) + log_T;
    return res;
}

/* add the function to calculate the log of the likelihood for the prior on the branches of the tree.
 * One thing to notice is that we don't need to do all this calculations as at each step only one branch changes and
 * we can calculate that one branch change and take off and add the changes to get the new value */
double log_b_likelihood(double lambda, int *b, double *b_length)
{
    int i;
    double res = 0.0;
    for(i = 0; i < number_branches-1; i++) /* this has to be number_branches-1, as number_branches includes the root, we want non root branches. */
        if(b[i])
            res += log( 1.0 - exp(-lambda * b_length[i]));
        else
            res += -lambda * b_length[i];
    return res;
}

double propose_new_lambda(double old_lambda)
{
    double new_lambda;
    double sigma = 1.0;
    new_lambda = old_lambda + gsl_ran_gaussian(r,sigma);
    return new_lambda;
}

/* sum of all three segments of the likelihood. */
double log_likelihood_all(int **counts, int num_sec, double lambda, int *b, double *b_lengths)
{
    double res;
    res = log_likelihood(counts, num_sec) + log_lambda_likelihood(lambda) + log_b_likelihood(lambda, b,b_lengths);
    return res;
}

/*log likelihood when only b has changed. */
double log_likelihood_lambda_constant(int **counts, int num_sec, int *b, double lambda,double *b_length)
{
    double res;
    res = log_likelihood(counts, num_sec) + log_b_likelihood(lambda, b,b_length);
    return res;
}
/* log likelihood when only lambda has changed. */
double log_likelihood_b_constant(double lambda,int *b, double *b_length)
{
    double res;
    res = log_lambda_likelihood(lambda) + log_b_likelihood(lambda, b,b_length);
    return res;
}
