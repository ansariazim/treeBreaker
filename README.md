# TreeBreaker
TreeBreaker can infer the evolution of a phenotype distribution on a phylogenetic tree, and divide the tree into segments where this distribution is constant.

For more details about how TreeBreaker works, please see the manuscript "Bayesian Inference of the Evolution of a Phenotype Distribution on a Phylogenetic Tree" by M Azim Ansari and Xavier Didelot, http://biorxiv.org/content/early/2016/02/23/040980

**Bulding the project**
- To build the project you need to install gsl.
- The command to build the project is as follows (assuming that you are in src directory):
```bash
    gcc -lgsl treeBreaker.c ../libs/knhx.c -o treeBreaker
```
- If you need to install gsl locally you can do it as follows:  
```bash
    wget ftp://ftp.gnu.org/gnu/gsl/gsl-1.16.tar.gz  
    tar xfz gsl-1.16.tar.gz  
    rm gsl-1.16.tar.gz  
    cd gsl-1.16 && ./configure && make  
  
    cd ../src  
    gcc -I ../gsl-1.16  -lm treeBreaker.c ../libs/knhx.c -o treeBreaker ../gsl-1.16/.libs/libgsl.a   
```
**Example run**
```bash
    ./treeBreaker ../testData/testTree.newick ../testData/phenoTestFile.txt outfile
```
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

**Output file**

The output file contains: 
- One line for each MCMC sample, with columns indicating whether branches have a changepoint on them (1) or not (0). To know which column matches which node on the tree, refer to the "index" attribute described below. The last column contains the values of lambda.
- The last line contains the newick string of the tree with comments in square brackets. Each node can have 3 attributes:
    - index: shows which column in the first part of the output file matches the node. The index starts at 0.
    - pheno: for leaf nodes, the phenotype that was read from the input phenotype file.
    - posterior: the posterior probability of having a change point on the branch above the node.

