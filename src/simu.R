library('ape')

#Simulate data under TreeBreaker model 
simu <- function(n=100,inputtreefile=NULL,prefixout='simu',lambda=NULL,nbChangingBranches=NULL,seed=NULL)  {
  if (!is.null(seed)) set.seed(seed)
  
  if (!is.null(inputtreefile)) {
    #Use tree given as input
    tr=read.tree(inputtreefile) 
    n=length(tr$tip.label)
    
  } else {  
    #Creates a random coalescent tree with n leaves
    tr<-list()
    tr$Nnode<-n-1
    tr$tip.label<-as.character(1:n)
    tr$edge<-matrix(0,n*2-2,2)
    tr$edge.length<-rep(0,n*2-2)
    toadd<-1:n
    time<-0
    tim<-rep(0,n)
    iedge<-1
    inode<-2*n-1
    while (length(toadd)>1) {
      l<-length(toadd)
      time<-time+rexp(1,l*(l-1)/2)
      w<-sample(1:l,2)
      tr$edge[iedge,]<-c(inode,toadd[w[1]])
      tr$edge.length[iedge]<-time-tim[w[1]]
      iedge<-iedge+1
      tr$edge[iedge,]<-c(inode,toadd[w[2]])
      tr$edge.length[iedge]<-time-tim[w[2]]
      iedge<-iedge+1
      toadd<-c(toadd[-w],inode)
      inode<-inode-1
      tim<-c(tim[-w],time)
    }
    class(tr)<-'phylo'
    write.tree(tr,file=sprintf('%s.nwk',prefixout))
  }
  
  #Changing branches
  tr=reorder.phylo(tr,'postorder')
  if (is.null(lambda)) lambda=rexp(1)/sum(tr$edge.length)
  probs=1-exp(-lambda*tr$edge.length)
  if (!is.null(nbChangingBranches)) {
    changing=rep(F,length(probs))
    changing[sample(1:length(probs),nbChangingBranches,replace=F,prob=probs)]=T
  } else {
    changing=runif(length(tr$edge.length))<probs
  }
  
  #Create hidden variables
  h=rep(NA,n+n-1)
  h[n+1]=runif(1)
  for (i in seq(nrow(tr$edge),1,-1)) {
    if (changing[i]) h[tr$edge[i,2]]=runif(1) else 
      h[tr$edge[i,2]]=h[tr$edge[i,1]]
  }

  #Create phenotypes and write to file
  pheno=as.numeric(runif(n)<h[1:n])
  write.table(cbind(tr$tip.label,pheno),sprintf('%s.pheno',prefixout),quote=F,row.names=F,col.names=F)
} 

testOnSimulation=function(nbChangingBranches=1) {
  #First simulate a dataset with given number of changepoints
  simu(inputtreefile = '../testData/tree1000.nwk',prefix='../testData/simu',nbChangingBranches=nbChangingBranches)
  #Run analysis
  system(sprintf('../bin/treeBreaker ../testData/tree1000.nwk ../testData/simu.pheno ../testData/simu.out'))
  #Read output file
  t=read.table('../testData/simu.out',comment.char='(')
  states =as.matrix(t[,2:(ncol(t)-1)])
  lambdas=as.vector(t[,ncol(t)])
  #Calculate number of inferred changepoints and Bayes Factor and return both
  inferredNbChangingBranches=mean(rowSums(states)-1)
  bf=length(which(lambdas>0))/length(which(lambdas==0))
  return(c(inferredNbChangingBranches,bf))
}