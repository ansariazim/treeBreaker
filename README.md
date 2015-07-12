# treeBreaker
Algorithm to divide a phylogenetic tree into segments based on phenotypes at the leaves of the tree

**Bulding the project**
- To build the project you need to install gsl.
- The command to build the project is as follows (assuming that you are in src directory):

    gcc -lgsl treeBreaker.c ../libs/knhx.c -o ../bin/treeBreaker
    
**Example run**

    ../bin/treeBreaker ../testData/testTree.newick ../testData/phenoTestFile.txt outfile
**Input arguments**

There are 3 mandatory input arguments which should have the following order:
  1. Input newick tree
  2. Input pheno file which has to be a tab seperated file. The first column is the name of the isolate and the second column is the phenotype.
    The isolate names in the pheno file have to match the isolate names in the newick file. The phenotypes have to be integers starting with 0.
  3. Name of output file
    
There are also some optional arguments, which must be indicated before the three mandatory input arguments:

    -x NUM      Sets the number of iterations after burn-in (default is 50000)
    -y NUM      Sets the number of burn-in iterations (default is 50000)
    -z NUM      Sets the number of iterations between samples (default is 100)
    -S NUM      Sets the seed for the random number generator to NUM
    -v          Verbose mode

**Output files**
- change_points_mcmc_state_file.txt: one line per change points array. If a branch above a node has a change point on it, then it is set to 1
    and to 0 otherwise. The root branch (the last column) is always set to 1. To know which column matches which node on the the tree, refer to
    the "index" attribute of each node in the output newick string in "outfile".
- outfile: contains the newick string of the tree with comments in []. Each node can have 3 attributes:
    - index: shows which column in "change_points_mcmc_state_file.txt" matches the node. The index starts at 0.
    - pheno: for leaf nodes, the phenotype that was read from the input phenotype file.
    - posterior: the posterior probability of having a change point on the branch above the node.
- lambda_mcmc_state_file.txt: Each line shows a value of the parameter lambda sampled in the MCMC.

