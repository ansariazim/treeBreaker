library('ape')

simu <- function(n=100,inputtreefile=NULL,prefixout='simu',lambda=NULL,seed=NULL)  {
  #Simulate data under TreeBreaker model 
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
  
  #Create hidden variables
  if (is.null(lambda)) lambda=rexp(1)/sum(tr$edge.length)
  tr=reorder.phylo(tr,'postorder')
  h=rep(NA,n+n-1)
  h[n+1]=runif(1)
  for (i in seq(nrow(tr$edge),1,-1)) {
    probcp=exp(-lambda*tr$edge.length[i])
    if (runif(1)<probcp) h[tr$edge[i,2]]=runif(1) else 
      h[tr$edge[i,2]]=h[tr$edge[i,1]]
  }
  
  #Create phenotypes and write to file
  pheno=as.numeric(runif(n)<h[1:n])
  sink(sprintf('%s.pheno',prefixout))
  for (i in 1:n) cat(i,pheno[i],'\n')
  sink()
} 

