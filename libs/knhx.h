#ifndef KNHX_H_
#define KNHX_H_

#define KNERR_MISSING_LEFT   0x01
#define KNERR_MISSING_RGHT   0x02
#define KNERR_BRACKET        0x04
#define KNERR_COLON          0x08

/* tree structre. 
 * parent: is the index of parent node. For Root this is set to -1.
 * n: number of children nodes. Leaves have zero children node, all other nodes should have >= 2 children nodes.
 * child: array of ints with the index of children node. The number of entries is stored in n.
 * name: string containing the name of the node. If no name then empty string. If there is space, the name will break there.
 * index: This is sorted index that I use. For leaves it goes from 0 to n-1 and for nodes from n to 2n-2
 * pheno: the phenotype of the leaves. For internal nodes this will be set to -1. They have to be in the range 0 to K-1 if there are K distinct phenotypes.
 * d: branch length or distance to parent node. If no distance then it is set to -1. 
 * posterior: Posterior probability of having a change point on the branch above the node.*/
typedef struct {
    int parent, n,index, pheno;
    int *child;
    char *name;
    double d;
    double posterior;
} knhx1_t;

#ifndef KSTRING_T
#define KSTRING_T kstring_t
typedef struct __kstring_t {
    size_t l, m;
    char *s;
} kstring_t;
#endif

#ifdef __cplusplus
extern "C" {
#endif
    /* parse the newick string and return the tree structre that is stored in knhx1_t. The tree is an array of of _n nodes.
     * char *nhx: string containing the newick format. There should be no spaces and no newline characters in the string.
     * int *_n: It will contain the number of nodes on the tree including the root. They will be indexed from 0 to _n-1.
     * int *_error: not sure about this one.*/
    knhx1_t *kn_parse(const char *nhx, int *_n, int *_error);

    /* write the output tree with posterior probabilities and other stuff in comments [] for all nodes.
     * knhx1_t *node: this points to the start of the tree which is node[0].
     * int root: the index of the root node which should be from the kn_parse function output (_n-1).
     * kstring_t *s: a structure which returns the reformatted newick string in s.s */
    void kn_format(const knhx1_t *node, int root, kstring_t *s);

    /* return one if node is leaf, zero otherwise. */
    int isleaf(knhx1_t *node);

    /* return the number of leaves on the tree. 
     * knhx1_t *tree: pointer to the first node of the tree which is node 0 from the return value of kn_parse.
     * int n_nodes: number of nodes on the tree from the return value of kn_parse.*/
    int get_number_leaves(knhx1_t *tree, int n_nodes);

    /* functions to read newick string from file.*/
    
    /* read newick string from file and return a pointer to it. remove the ; from the end of the string.
     * char *file_name: string that contains the path to the file name. We also check to ensure the string finishes with
     * a ";". If it does not then exit. if it does then remove it and return the rest of the string.*/
    char * get_newick_from_file(char *file_name);
    
    /* get the leaves under each node (terminal and internal) and the number of leaves under the node in the correct array indexed. 
     * par: it should come from function get_tree_data. It is the parent array with the correct indices.
     * leaves_under: the function will return a pointer 2d array. Each array is the list of leaves under that node.
     * n_leaves_under: the function will return this 1d array. each number tells us how many leaves are under that node.
     * number_leaves: the number of leaves on the tree provided by us.
     * n_nodes: the number of nodes on the tree. provided by us.*/
    void get_leaves_under(int *par, int ***leaves_under, int **n_leaves_under, const int number_leaves, const int n_nodes);

    /* read the phenotype file and return the results. output should be the number of phenotypes
     * in the file. Assume that file is a tab seperated file with no header
     * The first column is the name of the leaf and the second column is the phenotype which must be an integer and starts at 0 to K.
     * The names must be the same as the one that is given by the newick string.
     * output: number of phenotypes present.
     * file_name: file path that contains the name of the phenotpyes file.
     * names: the names of the leaves of the tree in an array, will be filled in by function.
     * phenos: pointer to array that will be completed by the function.
     * number_leaves: */
    int get_pheno_from_file(char *file_name, char ***names, int **phenos, const int number_leaves);
    
    /* Given the names and phenos information go to the tree and set the pheno field for the leaves
     *  of the tree. For the internal nodes set the pheno to -1.
     *  tree: the tree structre.
     *  n_nodes: the nubmer of nodes on the tree
     *  number_leaves: the number of leaves on the tree.
     *  names: 2d array of names of the leaves of the tree.
     *  phenos: 1d array of the phenos which matches the names.*/
    void set_pheno_in_tree(knhx1_t *tree, const int n_nodes, const int number_leaves, char **names, int *phenos);

    /* get the arrays for the par and dists and phenos. This function has to be
     * called after set_pheno_in_tree_function and uses the output of that function.
     * par and dists should be n_nodes long. phenos should be number_leaves long.
     * tree: is the pointer to the tree structre that we got from kn_parse
     * n_nodes: is the number of nodes on the tree. It comes from kn_parse.
     * par: a pointer to the parents array. the parent of each node are given in the array. I use my own ->index. for root set to -1.
     * dist: the distance of each node to its parent. For root we set this to -1.
     * phenos: the phenos according to our new indexing. */
    void get_tree_data(knhx1_t *tree, const int n_nodes, const int number_leaves, int **par, double **dist, int **phenos);





#ifdef __cplusplus
}
#endif

#endif
