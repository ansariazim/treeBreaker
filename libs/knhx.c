#include <ctype.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include "knhx.h"
#define BUFF_SIZE 100 /* buffer size when reading the text file */

typedef struct {
    int error, n, max;
    knhx1_t *node;
} knaux_t;

static int min(const int *arr, int length) {
    int i;
    int minimum = arr[0];
    for (i = 1; i < length; ++i) {
        if (minimum > arr[i]) {
            minimum = arr[i];
        }
    }
    return minimum;
}

static int max(const int *arr, int length) {
    int i;
    int maximum = arr[0];
    for (i = 1; i < length; ++i) {
        if (maximum < arr[i]) {
            maximum = arr[i];
        }
    }
    return maximum;
}

static int find_first(const int item, const int *arr, int length) {
    int i;
    for (i = 1; i < length; ++i) {
        if (item == arr[i]) {
            return i;
        }
    }
    return -1;
}

static inline char *add_node(const char *s, knaux_t *aux, int x)
{
    char *p, *nbeg, *nend = 0;
    knhx1_t *r;
    if (aux->n == aux->max) {
        aux->max = aux->max? aux->max<<1 : 8;
        aux->node = (knhx1_t*)realloc(aux->node, sizeof(knhx1_t) * aux->max);
    }
    r = aux->node + (aux->n++);
    r->n = x; r->parent = -1;
    for (p = (char*)s, nbeg = p, r->d = -1.0; *p && *p != ',' && *p != ')'; ++p) {
        if (*p == '[') {
            if (nend == 0) nend = p;
            do ++p; while (*p && *p != ']');
            if (*p == 0) {
                aux->error |= KNERR_BRACKET;
                break;
            }
        } else if (*p == ':') {
            if (nend == 0) nend = p;
            r->d = strtod(p + 1, &p);
            --p;
        } else if (!isgraph(*p)) if (nend == 0) nend = p;
    }
    if (nend == 0) nend = p;
    if (nend != nbeg) {
        r->name = (char*)calloc(nend - nbeg + 1, 1);
        strncpy(r->name, nbeg, nend - nbeg);
    } else r->name = strdup("");
    return p;
}

/*function defined by me to get things ready for getting the par, children and distance arrays for each node.
 * This function goes through all the nodes and sets the tree->index such that leaves are marked from 0 to n-1 and 
 * the internal nodes are marked from n to 2n-2*. this function is called by kn_parse function. We don need to use it 
 * by ourselves.*/

static void set_index(knhx1_t *tree, int n_nodes){
    knhx1_t *p ;
    int i, leaf_count,internal_node_count, n_leaves;
    n_leaves = get_number_leaves(tree, n_nodes);
    leaf_count = 0;
    internal_node_count = n_leaves;
    for (i = 0; i < n_nodes; i++) {
        p = tree + i;
        if (isleaf(p)){
            p->index = leaf_count;
            leaf_count++;
        }else{
            p->index = internal_node_count;
            internal_node_count++;
        }
    }
}

knhx1_t *kn_parse(const char *nhx, int *_n, int *_error)
{
    char *p;
    int *stack, top, max;
    knaux_t *aux;
    knhx1_t *ret;

#define __push_back(y) do {										\
    if (top == max) {										\
        max = max? max<<1 : 16;								\
        stack = (int*)realloc(stack, sizeof(int) * max);	\
    }														\
    stack[top++] = (y);										\
} while (0)													\

    stack = 0; top = max = 0;
    p = (char*)nhx;
    aux = (knaux_t*)calloc(1, sizeof(knaux_t));
    while (*p) {
        while (*p && !isgraph(*p)) ++p;
        if (*p == 0) break;
        if (*p == ',') ++p;
        else if (*p == '(') {
            __push_back(-1);
            ++p;
        } else if (*p == ')') {
            int x = aux->n, m, i;
            for (i = top - 1; i >= 0; --i)
                if (stack[i] < 0) break;
            m = top - 1 - i;
            p = add_node(p + 1, aux, m);
            aux->node[x].child = (int*)calloc(m, sizeof(int));
            for (i = top - 1, m = m - 1; m >= 0; --m, --i) {
                aux->node[x].child[m] = stack[i];
                aux->node[stack[i]].parent = x;
            }
            top = i;
            __push_back(x);
        } else {
            __push_back(aux->n);
            p = add_node(p, aux, 0);
        }
    }
*_n = aux->n;
*_error = aux->error;
ret = aux->node;
free(aux); free(stack);
set_index(ret,*_n);
return ret;
}

#ifndef kroundup32
#define kroundup32(x) (--(x), (x)|=(x)>>1, (x)|=(x)>>2, (x)|=(x)>>4, (x)|=(x)>>8, (x)|=(x)>>16, ++(x))
#endif

static inline int kputsn(const char *p, int l, kstring_t *s)
{
    if (s->l + l + 1 >= s->m) {
        s->m = s->l + l + 2;
        kroundup32(s->m);
        s->s = (char*)realloc(s->s, s->m);
    }
    memcpy(s->s + s->l, p, l);
    s->l += l; s->s[s->l] = 0;
    return l;
}

static inline int kputc(int c, kstring_t *s)
{
    if (s->l + 1 >= s->m) {
        s->m = s->l + 2;
        kroundup32(s->m);
        s->s = (char*)realloc(s->s, s->m);
    }
    s->s[s->l++] = c; s->s[s->l] = 0;
    return c;
}

static void format_node_recur(const knhx1_t *node, const knhx1_t *p, kstring_t *s, char *numbuf)
{
    if (p->n) {
        int i;
        kputc('(', s);
        for (i = 0; i < p->n; ++i) {
            if (i) kputc(',', s);
            format_node_recur(node, &node[p->child[i]], s, numbuf);
        }
        kputc(')', s);
        if (p->name) kputsn(p->name, strlen(p->name), s);
        if (p->d >= 0) {
            sprintf(numbuf, ":%g", p->d);
            kputsn(numbuf, strlen(numbuf), s);
        }
    } else kputsn(p->name, strlen(p->name), s);
}

void kn_format(const knhx1_t *node, int root, kstring_t *s) // TODO: get rid of recursion
{
    char numbuf[128];
    format_node_recur(node, &node[root], s, numbuf);
}

int isleaf( knhx1_t *tree){
    return !(tree->n);
}

int get_number_leaves(knhx1_t *tree, int n_nodes){
    int n = 0, i;
    for(i = 0; i<n_nodes;i++)
        if (isleaf(tree+i))
            n++;
    return n;
}

/* remove white space and nongraphical characters from string that is read by read_text_file.
 * old_str: string that is read by read_text_file. We will remove all white spaces from the file and return a new string
 * with white spaces removed. */
static char * remove_space(char *old_str){
    int str_len = strlen(old_str);
    char *temp ; 
    char *new_str;
    int i, j, k;
    if((temp = malloc(str_len+1)) == NULL){
        fprintf(stderr,"Out of memory when removing white spaces from the string.\n");
        exit(EXIT_FAILURE);
    }
    j = 0;
    for(i = 0; i < str_len; i++){
        if (!(isspace(old_str[i])) && isgraph(old_str[i])){
            temp[j] = old_str[i];
            j++;
        }
    }
    temp[j] = '\0';
    if((new_str = realloc(temp,strlen(temp)+1)) == NULL) {
        fprintf(stderr,"Out of memory when removing white spaces from the string.\n");
        free(temp);
        exit(EXIT_FAILURE);
    }
    return new_str;
}

/* function that reads the content of a text file. The file cannot contain the null character
 * it will return a pointer to the string that contains the string.
 * char *file_name: string containing the file path.*/
static char * read_text_file(char *file_name){

    char str[BUFF_SIZE]={'\0'};
    int size;
    FILE * fp;
    int  i,j,n1,n2;
    char *final= malloc(1);
    char *temp;
    final[0] = '\0';

    if ((fp = fopen(file_name,"r")) == NULL){
        fprintf(stderr,"Cannot open file %s.\n",file_name);
        free(final);
        exit(EXIT_FAILURE);
    }
    /*    fprintf(stdout,"file %s was openned.\n",file_name);*/

    while (fgets(str, BUFF_SIZE ,fp) != NULL){
        n1 = strlen(final);
        n2 = strlen(str);
        size = n1 + n2 +1;
        if ((temp = realloc(final, size)) == NULL){
            fprintf(stderr,"Out of memory when reading file %s.\n",file_name);
            fclose(fp);
            free(final);
            exit(EXIT_FAILURE);
        }

        final = temp;
        j = 0;
        for(i = n1; i < n1+n2; i++){
            final[i] = str[j];
            j++;
        }
        final[n1+n2] = '\0';
    }
    if(ferror(fp)){
        fclose(fp);
        free(final);
        fprintf(stderr,"There was an error when reading file %s.\n",file_name);
        exit(EXIT_FAILURE);
    }
    if(feof(fp)){
        fclose(fp);
        return final;
    }
    fprintf(stderr,"The program should not get here. There is an error in reading the file.\n");
    exit(EXIT_FAILURE);
    /*fprintf(stderr,"We should not get here. There is an error in the logic of file.\n");*/
}

char * get_newick_from_file(char *file_name){
    char *temp;
    char *newick_str;

    temp = read_text_file(file_name);
    newick_str = remove_space(temp);
    /* let's ensure the string ends with ; if not then not valid newick str. */
    if (newick_str[strlen(newick_str)-1] != ';'){
        fprintf(stderr, "Not a valid newick string. A newick string should finish with \";\".\n");
        exit(EXIT_FAILURE);
    }
    newick_str[strlen(newick_str)-1] = '\0'; /* remove the ";" character from the end of the string.*/
    return newick_str;
}

/* This function trims the white spaces from the begining and end of the string and returns it back. 
 * str: input string which will be trimed and returned back.*/
char *trim(char *str)
{
    size_t len = 0;
    char *frontp = str;
    char *endp = NULL;

    if( str == NULL ) { return NULL; }
    if( str[0] == '\0' ) { return str; }

    len = strlen(str);
    endp = str + len;

    /* Move the front and back pointers to address the first non-whitespace
     * characters from each end.
     */
    while( isspace(*frontp) ) { ++frontp; }
    if( endp != frontp )
    {
        while( isspace(*(--endp)) && endp != frontp ) {}
    }

    if( str + len - 1 != endp )
        *(endp + 1) = '\0';
    else if( frontp != str &&  endp == frontp )
        *str = '\0';

    /* Shift the string so that it starts at str so that if it's dynamically
     * allocated, we can still free it on the returned pointer.  Note the reuse
     * of endp to mean the front of the string buffer now.
     */
    endp = str;
    if( frontp != str )
    {
        while( *frontp ) { *endp++ = *frontp++; }
        *endp = '\0';
    }
    return str;
}

/* Function that reads the phenotypes file and returns the result in pheno array.
 * file_name: path to the file.
 * pheno: pointer to array that will be filled in by the function. 
 * number_leaves: number of leaves on the tree.
 * */
/*
   int get_pheno_from_file(char *file_name, char ***names, int **phenos, int number_leaves){
   char *temp_line;
   int str_length;
   char *str;
   int temp_pheno, num_args, line_counter,i,j;
   char *temp_name ;
   char **temp_names;
   int *temp_phenos, counter;

   printf("I got to here.");
   if ((temp_names = malloc(number_leaves *sizeof(*temp_names))) == NULL){
   fprintf(stderr,"Out of memory when trying to read the phenotypes file.\n");
   exit(EXIT_FAILURE);
   }
   if ((temp_phenos = malloc(number_leaves *sizeof(int))) == NULL){
   fprintf(stderr,"Out of memory when trying to read the phenotypes file.\n");
   exit(EXIT_FAILURE);
   }
   for (i = 0; i<number_leaves; i++){
   temp_phenos[i] = -1;
   }

   str = read_text_file(file_name);
   trim(str);
   str_length = strlen(str);
   if ((str=realloc(str, str_length + 2)) == NULL){
   fprintf(stderr,"Out of memory when reading file %s.\n",file_name);
   free(str);
   exit(EXIT_FAILURE);
   }
   str[str_length]= '\n';
   str[str_length+1] = '\0';
   str_length = strlen(str);

   if ((temp_line = malloc(str_length+1)) == NULL){ 
   fprintf(stderr,"Out of memory when reading file %s.\n",file_name);
   free(str);
   exit(EXIT_FAILURE);
   }

   line_counter = 0;
   counter = 0;
   for(i = 0; i<str_length; i++){
   if (str[i] != '\n'){
   temp_line[counter] = str[i];
   counter++;
   }else if (str[i]=='\n'){
   temp_line[counter] = '\0';
   counter = 0;
   num_args= sscanf(temp_line,"%s\t%d",temp_name,&temp_pheno);
   if (num_args != 2){
   fprintf(stderr,"The phenotype file is not in the correct format.\n");
   exit(EXIT_FAILURE);
   }
   if (line_counter > number_leaves){
   fprintf(stderr,"There are more lines than the number of leaves on the tree in the phenotype file.\n");
   exit(EXIT_FAILURE);
   }
   if ((temp_names[line_counter] = malloc(strlen(temp_name)+1)) == NULL){
   fprintf(stderr,"Out of memory when trying to read the phenotypes file.\n");
   exit(EXIT_FAILURE);
   }
   strcpy(temp_names[line_counter], temp_name);
   temp_phenos[line_counter] = temp_pheno;
   line_counter++;
   }

   }
   return 0;
   }*/

/* An earlier version of the above function.*/

int get_pheno_from_file(char *file_name, char ***names, int **phenos, const int number_leaves){
    char str[BUFF_SIZE]={'\0'};
    int temp_pheno, num_args, line_counter,i;
    char temp_name[BUFF_SIZE] = {'\0'};
    FILE * fp;
    char **temp_names;
    int *temp_phenos;
    int maximum;

    if ((temp_names = malloc(number_leaves *sizeof(*temp_names))) == NULL){
        fprintf(stderr,"Out of memory when trying to read the phenotypes file.\n");
        exit(EXIT_FAILURE);
    }
    if ((temp_phenos = malloc(number_leaves *sizeof(int))) == NULL){
        fprintf(stderr,"Out of memory when trying to read the phenotypes file.\n");
        exit(EXIT_FAILURE);
    }

    for (i = 0; i<number_leaves; i++){
        temp_phenos[i] = -1;
    }
    if ((fp = fopen(file_name,"r")) == NULL){
        fprintf(stderr,"Cannot open file %s.\n",file_name);
        exit(EXIT_FAILURE);
    }

    line_counter = 0;
    while (fgets(str, BUFF_SIZE ,fp) != NULL){ /* assuming that the lines are less than 100 chars long. */
        trim(str);
        if (str[0]){

            num_args= sscanf(str,"%s\t%d",temp_name,&temp_pheno);
            if (num_args != 2){
                fprintf(stderr,"The phenotype file is not in the correct format.\n");
                exit(EXIT_FAILURE);
            }
            if (line_counter > number_leaves){
                fprintf(stderr,"There are more lines than the number of leaves on the tree.\n");
                exit(EXIT_FAILURE);
            }
            if ((temp_names[line_counter] = malloc(strlen(temp_name)+1)) == NULL){
                fprintf(stderr,"Out of memory when trying to read the phenotypes file.\n");
                exit(EXIT_FAILURE);
            }
            strcpy(temp_names[line_counter], temp_name);
            temp_phenos[line_counter] = temp_pheno;
            line_counter++;
        }
    } 

    if(ferror(fp)){
        fclose(fp);
        fprintf(stderr,"There was an error when reading file %s.\n",file_name);
        exit(EXIT_FAILURE);
    }
    if(feof(fp)){
        fclose(fp);
        *phenos = temp_phenos;
        *names = temp_names;
        /* let's check that the phenos start from 0 and go to K without any missing numbers in between.*/
        if (min(*phenos,number_leaves) != 0){
            fprintf(stderr,"Phenotypes have to be integer starting from 0. If there are \
                    K distinct phenotypes then all numbers between 0 and K-1 have to be used.\n");
            exit(EXIT_FAILURE);
        }else{
            maximum = max(*phenos,number_leaves);
            for(i = 1; i<= maximum; i++){
                if (find_first(i,*phenos,number_leaves) == -1){
                    fprintf(stderr,"Phenotypes have to be integers. If there are \
                            K distinct phenotypes then all numbers between 0 and K-1 have to be used.\n");
                    exit(EXIT_FAILURE);
                }
            }
        }
        return maximum+1;
    }
    fprintf(stderr,"The program should not get here. There is an error in reading the file.\n");
    exit(EXIT_FAILURE);
}

void get_tree_data(knhx1_t *tree, const int n_nodes, const int number_leaves, int **pars, double **dists, int **phenos){
    int i;
    int * temp_pars;
    double * temp_dists;
    int *temp_phenos;
    knhx1_t *p, *p_temp;

    if ((temp_pars = malloc(n_nodes * sizeof(int))) ==NULL){
        fprintf(stderr,"Out of memory when processing the newick tree.\n");
        exit(EXIT_FAILURE);
    }
    if ((temp_dists = malloc(n_nodes * sizeof(double))) ==NULL){
        fprintf(stderr,"Out of memory when processing the newick tree.\n");
        exit(EXIT_FAILURE);
    }
    if ((temp_phenos = malloc(number_leaves * sizeof(int))) ==NULL){
        fprintf(stderr,"Out of memory when processing the newick tree and its phenotypes.\n");
        exit(EXIT_FAILURE);
    }

    for(i = 0; i<n_nodes-1; i++){ /* while not at the root */
        p = tree + i;
        p_temp = tree + p->parent;
        temp_pars[p->index] = p_temp->index;
        temp_dists[p->index] = p->d;
        if (isleaf(p)){
            temp_phenos[p->index] = p->pheno;
        }
    }
    temp_pars[n_nodes - 1] = -1; /* for the root set to -1. */
    temp_dists[n_nodes - 1] = -1;
    *pars = temp_pars;
    *dists = temp_dists;
    *phenos = temp_phenos;
}

void get_leaves_under(int *par, int ***leaves_under, int **n_leaves_under, const int number_leaves, const int n_nodes){
    int p, leaf, counter,i;
    int **tmp_leaves_under = NULL;
    int *tmp_n_leaves_under = NULL;

    if((tmp_leaves_under = malloc( n_nodes* sizeof(*tmp_leaves_under))) == NULL){
        fprintf(stderr, "Out of memory when parsing the newick tree.\n");
        exit(EXIT_FAILURE);
    }
    if ((tmp_n_leaves_under = calloc(n_nodes, sizeof(int))) == NULL){
        fprintf(stderr, "Out of memory when parsing the newick tree.\n");
        exit(EXIT_FAILURE);
    }

    for(i = 0; i < n_nodes; i++){
        if ((tmp_leaves_under[i] = calloc(number_leaves, sizeof(int))) == NULL){
            fprintf(stderr, "Out of memory when parsing the newick tree.\n");
            exit(EXIT_FAILURE);
        }
    }

    for (leaf = 0; leaf < number_leaves; leaf++){
        p = leaf;
        while(p != -1){
            tmp_leaves_under[p][tmp_n_leaves_under[p]] = leaf;
            tmp_n_leaves_under[p]++;
            p = par[p];
        }
    }
    *leaves_under = tmp_leaves_under;
    *n_leaves_under = tmp_n_leaves_under;
}

void set_pheno_in_tree(knhx1_t *tree, const int n_nodes, const int number_leaves, char **names, int *phenos){
    int i,j;
    int flag;
    knhx1_t *p;
    for (i = 0; i < n_nodes; i++) {
        p = tree + i;
        if (isleaf(p)){
            flag = 1;
            for(j = 0; j<number_leaves;j++){
                if (strcmp(p->name,names[j]) == 0){
                    p->pheno = phenos[j];
                    flag =0;
                    break;
                }
            }
            if (flag){
                fprintf(stderr,"%s node had no phenotype information in the phenotype file. Check the name spelling and file format.\n",p->name);
                exit(EXIT_FAILURE);
            }
        }else{
            p->pheno= -1;
        }
    }

}


/*#define KNHX_MAIN */
#ifdef KNHX_MAIN
int main(int argc, char *argv[])
{
    char *s = "((a[abc],d1)x:0.5,((b[&&NHX:S=MOUSE],h2)[&&NHX:S=HUMAN:B=99][blabla][&&NHX:K=foo],c))";
    int *par = NULL;
    double *dists = NULL;
    knhx1_t *node;
    int i, j, n, error;
    kstring_t str;
    int **leaves = NULL;
    int *n_leaves = NULL;
    int num_leaves;

    node = kn_parse(s, &n, &error);
    get_tree_data(node, n,&par, &dists);
    num_leaves = get_number_leaves(node, n);


    for (i = 0; i < n; ++i) {
        knhx1_t *p = node + i;
        printf("[%2d] [%2d] %s\t%d\t%d\t%g", i,p->index, p->name, p->parent, p->n, p->d);
        for (j = 0; j < p->n; ++j)
            printf("\t%d", p->child[j]);
        putchar('\n');
    }

    for(i = 0; i<n;i++){
        printf("[%4d]\t %d\t %e\n",i,par[i],dists[i]);
    }

    get_leaves_under(par, &leaves, &n_leaves, num_leaves, n);
    printf("After getting the leaves we have:\n");
    for(i = 0; i < n; i++){
        printf("leaves under branch %d are: ",i);
        for(j = 0; j<n_leaves[i]; j++)
            printf("%4d ",leaves[i][j]);

        printf("\n");
    }


    str.l = str.m = 0; str.s = 0;
    kn_format(node, n-1, &str);
    puts(str.s);
    free(str.s);
    return 0;
}
#endif
