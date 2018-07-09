#' treeBreaker main function
#' @param inputfile_tree newick file to use as input tree
#' @param inputfile_phenotype tab file to use as input phenotype
#' @param outputfile output file
#' @param x the number of iterations after burn-in (default is 500000)
#' @param y the number of burn-in iterations (default is 500000)
#' @param z the number of iterations between samples (default is 1000)
#' @param seed the seed for the random number generator (default is a random seed)
#' @param verbose verbose mode
#' @return success status
#' @export
#' @importFrom Rcpp evalCpp
#' @useDynLib treeBreaker
treeBreaker = function(inputfile_tree,inputfile_phenotype,outputfile,x=500000,y=500000,z=1000,seed=NA,verbose=F) {
  args=c('-x',sprintf('%d',x),'-y',sprintf('%d',y),'-z',sprintf('%d',z),inputfile_tree,inputfile_phenotype,outputfile)
  if (verbose) args=c('-v',args)
  x=print('\n')
  if (!is.na(seed)) args=c('-S',sprintf('%d',seed),args)
  mainR(args)
}
