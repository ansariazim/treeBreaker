# treeBreaker
Algorithm to divide a phylogenetic tree into segments based on phenotypes at the leaves of the tree

**Bulding the project**
- To build the project you need to install gsl.
- The command to build the project is as follows (assuming that you are in src directory):

    gcc -lgsl treeBreaker.c ../libs/knhx.c -o ../bin/treeBreaker
    
**Example run**

    ../bin/treeBreaker ../testData/testTree.newick ../testData/phenoTestFile.txt 1000000
**Input arguments**

There are 3 input arguments which should have the following order:
  1. Input newick tree
  2. Input pheno file which has to be a tab seperated file. The first column is the name of the isolate and the second column is the phenotype.
    The isolate names in the pheno file have to match the isolate names in the newick file. The phenotypes have to be integers starting with 0.
  3. Number of steps in MCMC chain
    
    
