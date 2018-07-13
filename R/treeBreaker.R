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
treeBreaker = function(inputfile_tree,inputfile_phenotype,outputfile,x=500000,y=500000,z=1000,seed=NA,verbose=F) {
  args=c('-x',sprintf('%d',x),'-y',sprintf('%d',y),'-z',sprintf('%d',z),inputfile_tree,inputfile_phenotype,outputfile)
  if (verbose) args=c('-v',args)
  writeLines('     ')
  if (!is.na(seed)) args=c('-S',sprintf('%d',seed),args)
  status=mainR(args)
  if (status==1) return(NULL) else return(readOutFile(outputfile))
}

#' Reading an output file
#' @param outfile name of the output file
#' @return object containing output
#' @export
readOutFile = function(outfile) {
  #Read the tree and information attached to it
  res=list()
  res$tree=ape::read.tree(outfile)
  ntips=length(res$tree$tip.label)
  l=strsplit(res$tree$tip.label,'[{}|=]',perl=T)
  res$edge_index=rep(0,nrow(res$tree$edge))
  res$tip_pheno=rep(0,length(l))
  edge_posterior=rep(0,nrow(res$tree$edge))
  for (i in 1:length(l)) {
    res$tree$tip.label[i]=l[[i]][1]
    res$edge_index[which(res$tree$edge[,2]==i)]=as.numeric(l[[i]][3])
    res$tip_pheno[i]=as.numeric(l[[i]][5])
    res$edge_posterior[which(res$tree$edge[,2]==i)]=as.numeric(l[[i]][7])
  }
  l=strsplit(res$tree$node.label,'[{}|=]',perl=T)
  res$node_posterior=rep(0,length(l))
  for (i in 1:length(l)) {
    res$edge_index[which(res$tree$edge[,2]==i+ntips)]=as.numeric(l[[i]][3])
    res$edge_posterior[which(res$tree$edge[,2]==(i+ntips))]=as.numeric(l[[i]][5])
  }

  #Read rest of file
  t=read.table(outfile,comment.char='(')
  res$states =as.matrix(t[,2:(ncol(t)-1)])
  res$lambdas=as.vector(t[,ncol(t)])
  class(res)<-"resTreeBreaker"
  return(res)
}

#' Plotting function
#' @param x An object of class resTreeBreaker
#' @param type Type of plot to perform, can be cons (default), states, traces or correlation
#' @param ... Ignored
#' @return Invisible object
#' @export
plot.resTreeBreaker = function(x,type='cons',...) {
  if (type=='cons') {
    #Plot tree showing phenotype and posterior
    par(mfrow=c(1,1))
    ec=x$edge_posterior
    w=which(x$tree$edge[,1]==(ape::Ntip(x$tree)+1));if (length(w)==2) ec[w]=max(ec[w])
    ape::plot.phylo(x$tree,show.tip.label = F,edge.color=rgb(ec,0,0),edge.width=1+ec*10)
    ncols=length(unique(x$tip_pheno))
    ape::tiplabels(NULL,pch=16,col=rainbow(2*ncols)[ncols+x$tip_pheno])
  }

  if (type=='trace') {
    #Plot some MCMC traces
    par(mfrow=c(1,2))
    plot(x$lambdas,type = 'l',xlab='Sampled iteration',ylab='lambda')
    plot(rowSums(x$states)-1,type='l',xlab='Sampled iteration',ylab='Number of changepoints')
  }

  if (type=='states') {
    #Plot ten MCMC states
    par(mfrow=c(2,5))
    for (s in seq(nrow(x$states)/10,nrow(x$states),nrow(x$states)/10)) {
      ec=rep(0,nrow(x$tree$edge))
      for (i in 1:nrow(x$tree$edge)) {
        ec[i]=x$states[s,x$edge_index[i]+1]
      }
      w=which(x$tree$edge[,1]==(ape::Ntip(x$tree)+1));if (length(w)==2) ec[w]=max(ec[w])
      ape::plot.phylo(x$tree,show.tip.label = F,edge.color=rgb(ec,0,0),edge.width=1+ec*10)
      ncols=length(unique(x$tip_pheno))
      ape::tiplabels(NULL,pch=16,col=rainbow(2*ncols)[ncols+x$tip_pheno])
    }
  }

  if (type=='correlation') {
    #plot the correlation between the change points
    par(mfrow=c(1,1))
    image(t(x$states[,-ncol(x$states)]), axes=FALSE)
    axis(1, at=seq(0,1,length.out=10), labels= floor(seq(0,ncol(x$states)-1,length.out=10) ))
    axis(2, at=seq(0,1,length.out=10), labels= floor(rev( seq(0,nrow(x$states),length.out=10)) ))
  }

  return(invisible(x))
}


