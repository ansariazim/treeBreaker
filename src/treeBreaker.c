/***************************************************************************
 *   Copyright (C) 2015 by Azim Ansari and Xavier Didelot                  *
 *   ansari.azim@gmail.com and xavier.didelot@gmail.com                    *
 *                                                                         *
 *   This program is free software; you can redistribute it and/or modify  *
 *   it under the terms of the GNU General Public License as published by  *
 *   the Free Software Foundation; either version 2 of the License, or     *
 *   (at your option) any later version.                                   *
 *                                                                         *
 *   This program is distributed in the hope that it will be useful,       *
 *   but WITHOUT ANY WARRANTY; without even the implied warranty of        *
 *   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the         *
 *   GNU General Public License for more details.                          *
 *                                                                         *
 *   You should have received a copy of the GNU General Public License     *
 *   along with this program; if not, write to the                         *
 *   Free Software Foundation, Inc.,                                       *
 *   59 Temple Place - Suite 330, Boston, MA  02111-1307, USA.             *
 ***************************************************************************/

#include <stdio.h>
#include <float.h>
#include <ctype.h>
#include <stdlib.h>
#include <unistd.h>
#include <string.h>
#include <stdbool.h>
#include "../libs/knhx.h"
#include <time.h>
#include <math.h>
#include <gsl/gsl_rng.h>
#include <gsl/gsl_randist.h>
#include <gsl/gsl_sf_gamma.h>

static const char * help=
"\
        Usage: treeBreaker [OPTIONS] inputfile_tree inputfile_phenotype outputfile\n\
        \n\
        Options:\n\
        -x NUM      Sets the number of iterations after burn-in (default is 500000)\n\
        -y NUM      Sets the number of burn-in iterations (default is 500000)\n\
        -z NUM      Sets the number of iterations between samples (default is 1000)\n\
        -S NUM      Sets the seed for the random number generator to NUM\n\
        -v          Verbose mode\n\
        \n";

int propose_new_b(int *b,int num_branches, int * b_star);
int count_phenos(int **code, int par[], int b[], int pheno[], int **counts);
double log_likelihood(int ** counts,int num_sec);
double log_lambda_prior(double lambda);
double log_b_prior(double lambda, int *b, double *b_length);
double propose_new_lambda(double old_lambda);
double log_posterior(int **counts, int num_sec, double lambda, int *b, double *b_lengths);
double log_likelihood_lambda_constant(int **counts, int num_sec, int *b, double lambda,double *b_length);
double log_likelihood_b_constant(double lambda,int *b, double *b_length);
void set_posterior(unsigned long int *b_counts, unsigned long int denominator, knhx1_t *tree);
double m0_propose_lambda(void);
void m0_propose_b(int b_star[], double branches_len[], double lambda_star);
double calculate_log_evidence_model_0(int pheno[]);

int number_branches;
int number_leaves;
int number_phenotypes;
double T, log_T, inv_T;
gsl_rng *r;

int main(int argc, char *argv[]){
    knhx1_t *tree;
    int error;
    int i,j,k, mcmc_counter, changed_branch;
    kstring_t str;
    char *newick_str = NULL;
    FILE *fp;
    char **names;
    int *temp_phenos;
    int *phenos;
    int *parents;
    double *branches_len ;
    int **leaves_under;
    int *n_leaves_under;
    int *b, *b_star;
    int *sections, num_sec, recording_counter;
    int **counts,temp_counter=0;
    double  lambda, model_0_log_evidence, mh_model_0_all;
    unsigned long int *b_counts, denominator, model_0_counter=0, model_1_counter=0;
    int model_state;
    int postburn=500000;
    int burn=500000;
    int thin=1000;
    bool seeded=false;
    unsigned int seed=0;
    bool verbose=false;
    bool root_binary_flag = false;
    int c,node1_index, node2_index, long_node, short_node;
    double prob_move_0_to_1 = 0.5, prob_move_1_to_0 = 0.05;
    double log_prob_move_0_to_1 = log(prob_move_0_to_1), log_prob_move_1_to_0 = log(prob_move_1_to_0);
    int **code;
    double proposal_log_likelihood, proposal_log_b_prior, proposal_log_lambda_prior, lambda_star, old_log_likelihood, old_log_b_prior, old_log_lambda_prior;

    while ((c = getopt (argc, argv, "x:y:z:S:v")) != -1)
        switch (c)
        {
            case('x'):postburn=atoi(optarg);break;
            case('y'):burn=atoi(optarg);break;
            case('z'):thin=atoi(optarg);break;
            case('S'):seeded=true;seed=atoi(optarg);break;
            case('v'):verbose=true;break;
            default:fprintf(stderr,"Syntax error.\n%s",help);exit(EXIT_FAILURE);
        }
    if (argc-optind!=3) {fprintf(stderr,"Syntax error.\n%s",help);exit(EXIT_FAILURE);}
    r = gsl_rng_alloc(gsl_rng_mt19937);
    if (seeded==true) gsl_rng_set(r,seed); else gsl_rng_set(r,time(NULL));
    mcmc_counter = postburn+burn;
    recording_counter = burn;

    newick_str = get_newick_from_file(argv[optind++]);
    if (verbose) printf("%s\n",newick_str);

    tree = kn_parse(newick_str,&number_branches, &error);
    number_leaves = get_number_leaves(tree, number_branches);
    number_phenotypes = get_pheno_from_file(argv[optind++], &names, &temp_phenos, number_leaves);
    set_pheno_in_tree(tree, number_branches, number_leaves, names, temp_phenos);
    get_tree_data(tree, number_branches,number_leaves, &parents, &branches_len,&phenos);
    get_leaves_under(parents, &leaves_under, &n_leaves_under, number_leaves, number_branches);

    /* the way to have an unrooted tree is as follows:
     * check to make sure that the root node has only two nodes under it (binary).
     * set the length of the branch for the smaller of the two branches to zero.
     * set the lenght of the longer branch so that it equals to the sum of the length of the two branches.
     * at every iteration of mcmc check to see if the branch with zero length has been set and if that is the case then
     * reject the move and stay at the previous state. */

    if ((tree+number_branches-1)->n == 2){
        if (verbose) printf("The root has two branches under it.\n");
        root_binary_flag = true;
        node1_index = ((tree+(tree+number_branches-1)->child[0])->index);
        node2_index = ((tree+(tree+number_branches-1)->child[1])->index);
        if (branches_len[node1_index] >= branches_len[node2_index]){
            long_node = node1_index;
            short_node = node2_index;
        }else{
            long_node = node2_index;
            short_node = node1_index;
        }
        branches_len[long_node] += branches_len[short_node];
        branches_len[short_node] = 0;
        if (verbose) printf("The long branch now is: %f and short branch is: %f.\n",branches_len[long_node],branches_len[short_node]  );
    }

    char * output_filename = argv[optind++];

    if (verbose) for(i = 0; i<number_branches; i++){ 
        printf("[%3d]\t%3d\t%3d\t%4g", i, parents[i], n_leaves_under[i], branches_len[i]);
        for (j = 0; j < n_leaves_under[i]; ++j)
            printf("\t%d", leaves_under[i][j]);
        putchar('\n');
    }

    if ((code = malloc(number_branches *sizeof(int *))) == NULL){
        fprintf(stderr,"Out of memory. Could not assign memory when setting code.");
        exit(EXIT_FAILURE);
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
    num_sec = count_phenos(code, parents, b, phenos, counts);
    /* let's calculate the total length of branches on the tree. */
    T = 0.0;
    for(i = 0; i < number_branches-1; i++)
        T += branches_len[i];
    log_T = log(T);
    inv_T = 1 / T;
    model_0_log_evidence = calculate_log_evidence_model_0(phenos); /* calculate the evidence for model 0 */
    mh_model_0_all = model_0_log_evidence + log_prob_move_0_to_1;

    /* we always start in model 0 */
    model_state = 0;
    lambda = 0;
    old_log_likelihood = model_0_log_evidence;
    for(i = 0;i <num_sec;i++)
        for (j = 0; j<number_phenotypes;j++)
            counts[i][j] = 0;

    denominator = 0;
    if ((fp = fopen(output_filename,"w")) == NULL){
        fprintf(stderr,"Cannot open file %s.\n",output_filename);
        exit(EXIT_FAILURE);
    }

    /* let's do the mcmc now. */
    for(i = 0; i<mcmc_counter; i++){
        if (mcmc_counter>50 && (i)%(mcmc_counter/50)==0){
            printf("\b\b\b\b\b# %3d%%",i*100/mcmc_counter);
            fflush(0);
        }
        if (i+1==mcmc_counter){
            printf("\b\b\b\b\b# 100%%\n");
            fflush(0);
        }
        switch( model_state ){
            case 0: /* if in model 0 then: */
                if (gsl_rng_uniform(r) > prob_move_0_to_1) { /* if propose to stay in model 0 */
                    model_0_counter++;
                }else{ /* if propose to go to model 1 */
                    lambda_star = m0_propose_lambda();
                    m0_propose_b(b_star, branches_len, lambda_star); /* I have to make sure that I don't propose change points on the branch under root that has length zero.*/
                    num_sec = count_phenos(code, parents, b_star, phenos, counts);
                    proposal_log_likelihood = log_likelihood(counts, num_sec);
                    if (gsl_rng_uniform(r) <(exp(proposal_log_likelihood + log_prob_move_0_to_1 - mh_model_0_all))){ /* move to model 1 accepted */
                        lambda = lambda_star;
                        for(j = 0; j < number_branches-1; j++)
                            b[j] = b_star[j];
                        old_log_likelihood = proposal_log_likelihood;
                        old_log_b_prior = log_b_prior(lambda,b,branches_len);
                        old_log_lambda_prior = log_lambda_prior(lambda);
                        model_1_counter++;
                        model_state = 1;
                    }else{ /* move to 1 rejected, stay in model 0 */
                        model_0_counter++;
                    }
                    for(j = 0;j <num_sec;j++)
                        for (k = 0; k<number_phenotypes;k++)
                            counts[j][k] = 0;
                }
                break;
            case 1: /* if in model 1 then */
                if (gsl_rng_uniform(r) > prob_move_1_to_0) { /* if proposed to stay in model 1: */
                    model_1_counter++;
                    changed_branch = propose_new_b(b,number_branches, b_star);
                    if (!(root_binary_flag && (branches_len[changed_branch] < DBL_EPSILON))){
                        num_sec = count_phenos(code, parents, b_star, phenos, counts);
                        proposal_log_likelihood = log_likelihood(counts,num_sec);
                        proposal_log_b_prior = log_b_prior(lambda,b_star,branches_len);
                        if (gsl_rng_uniform(r) <(exp(proposal_log_likelihood + proposal_log_b_prior  - old_log_likelihood - old_log_b_prior ))){
                            b[changed_branch] = b_star[changed_branch];
                            old_log_likelihood = proposal_log_likelihood;
                            old_log_b_prior = proposal_log_b_prior;
                        }
                    }

                    lambda_star = propose_new_lambda(lambda);
                    if (lambda_star > 0.0){
                        proposal_log_b_prior = log_b_prior(lambda_star, b, branches_len);
                        proposal_log_lambda_prior = log_lambda_prior(lambda_star);
                        if (gsl_rng_uniform(r) <(exp(proposal_log_lambda_prior + proposal_log_b_prior - old_log_lambda_prior -old_log_b_prior))){
                            lambda = lambda_star;
                            old_log_b_prior = proposal_log_b_prior;
                            old_log_lambda_prior = proposal_log_lambda_prior;
                        }
                    }
                    for(j = 0;j <num_sec;j++)
                        for (k = 0; k<number_phenotypes;k++)
                            counts[j][k] = 0;
                }else{ /* if proposed to go from model 1 to model 0 */
                    if (gsl_rng_uniform(r) <(exp(mh_model_0_all - old_log_likelihood - log_prob_move_1_to_0))){
                        model_state = 0;
                        model_0_counter++;
                        lambda = 0.0; /* in model 0, lambda =0 */
                        for(j = 0; j < number_branches-1; j++) /* in model 0 there are no change points on the tree.*/
                            b[j] = 0;
                    }else{
                        model_1_counter++;
                    }

                }
                break;
        }

        if (i == recording_counter){
            for(j = 0;j <number_branches; j++)
                b_counts[j] += b[j];
            recording_counter += thin;
            denominator++;
            fprintf(fp,"#\t");
            for(k = 0; k < number_branches-1; k++)
                fprintf(fp,"%i\t",b[k]);
            fprintf(fp,"%i\t",b[number_branches-1]);
            fprintf(fp,"%g\n",lambda);
        }
    }

    set_posterior(b_counts, denominator, tree);
    str.l = str.m = 0; str.s = 0;
    kn_format(tree, number_branches - 1, &str);
    printf("sum of the model counters is %lu. Model 0: %lu and Model 1: %lu.\n",model_0_counter + model_1_counter,model_0_counter,model_1_counter);
    printf("Bayes factor of model 1 to 0 is:%f.\n",((double) model_1_counter) / model_0_counter);

    if (verbose) printf("Writing output file.\n");
    fprintf(fp,"%s\n",str.s);
    fclose(fp);
    free(str.s);

    /*if (verbose) printf("Counter for acceptance is: %d\n",temp_counter);*/
    if (verbose){
        for(i = 0; i<number_branches; i++)
            printf("[%d]\t%f\n",i,((double) b_counts[i])/denominator);
    }

    free(code);
    free(b_counts);
    free(b);
    free(b_star);
    free(sections);
    for(i = 0; i<number_branches; i++)
        free(counts[i]);
    free(counts);
    free(newick_str);
    for (i = 0; i<number_leaves;i++) free(names[i]);
    free(names);
    free(temp_phenos);
    free(parents);
    free(branches_len);
    free(phenos);
    for (i=0;i<number_branches;i++) free(leaves_under[i]);
    free(leaves_under);
    free(n_leaves_under);

    return 0;
}

/**************************************
 * given par and b and pheno arrays, we determine the number of sections on the tree and then
 * count the number of each phenotype in each section. Here we assume binary phenotypes.
 * We just simply start at the bottom of the tree (leaves) and go up until we reach a branch with a change point on it. We then set the leaves cluster number.
 * In addition we use an array of pointers that point to array of counts for each section. If a branch has a change point on it and it is reached by a leaf then 
 * we set the pointer to point counts and modify the counts as appropriate.
 */
int count_phenos(int **code, int par[], int b[], int pheno[], int **counts)
{
    int i, j, p;
    /*int set[N] = {0};*/
    for (i = 0; i<number_branches; i++)
        code[i] = NULL;

    j = 0;
    for (i = 0; i<number_leaves; i++)
    {
        p = i;
        while(!b[p])
            p = par[p];
        /*printf("At leaf %d and parent is %d.\n",i,p);*/

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

double log_lambda_prior(double lambda)
{
    double res;
    res = (-T * lambda) + log_T;
    return res;
}

/* add the function to calculate the log of the likelihood for the prior on the branches of the tree.
 * One thing to notice is that we don't need to do all this calculations as at each step only one branch changes and
 * we can calculate that one branch change and take off and add the changes to get the new value */
double log_b_prior(double lambda, int *b, double *b_length)
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
double log_posterior(int **counts, int num_sec, double lambda, int *b, double *b_lengths)
{
    double res;
    res = log_likelihood(counts, num_sec) + log_lambda_prior(lambda) + log_b_prior(lambda, b,b_lengths);
    return res;
}

/*log likelihood when only b has changed. */
double log_likelihood_lambda_constant(int **counts, int num_sec, int *b, double lambda,double *b_length)
{
    double res;
    res = log_likelihood(counts, num_sec) + log_b_prior(lambda, b,b_length);
    return res;
}
/* log likelihood when only lambda has changed. */
double log_likelihood_b_constant(double lambda,int *b, double *b_length)
{
    double res;
    res = log_lambda_prior(lambda) + log_b_prior(lambda, b,b_length);
    return res;
}

/* set the posterior field in the node structrue.*/
void set_posterior(unsigned long int *b_counts, unsigned long int denominator, knhx1_t *tree)
{
    int i;
    knhx1_t *p;
    for(i = 0; i<number_branches; i++){
        p = tree + i;
        p->posterior = ((double) b_counts[p->index]) / denominator;
    }
}
/* calculate the evidence [p(D|m_1)] for model r0.
 * under model 0 there are no change points on the tree. */
double calculate_log_evidence_model_0(int pheno[]){
    double res;
    int i, j, counter;
    res = gsl_sf_lngamma(number_phenotypes) - gsl_sf_lngamma(number_leaves + number_phenotypes);
    for (i = 0 ; i<number_phenotypes; i++){
        counter = 0;
        for (j = 0; j < number_leaves; j++){
            if (pheno[j] == i)
                counter += 1;
        }
        res += gsl_sf_lngamma(counter + 1);
    }
    return res;
}
/* propose lambda from the prior on lambda when in model 0. */
double m0_propose_lambda(void){
    return gsl_ran_exponential(r,inv_T);
}

/*propose b from the prior on b when in model 0 
 * Here I need to check to make sure that branches are not of length zero.
 * I have set the length of one of the branches under the root to zero and
 * exp(0) is not defined. This could cause a problem on some machines. */
void m0_propose_b(int b_star[], double branches_len[], double lambda_star){
    int i;
    for (i = 0; i < number_branches-1; i++){
        if ( (branches_len[i] > DBL_EPSILON) && (gsl_rng_uniform(r) > exp(-lambda_star * branches_len[i])))
            b_star[i] = 1;
        else
            b_star[i] = 0;
    }
    b_star[number_branches-1] = 1;
}
