rm(list=ls())
library('ape')

#Read the tree and information attached to it
tree=read.tree('outfile')
ntips=length(tree$tip.label)
l=strsplit(tree$tip.label,'[\\]\\[\\|=]',perl=T)
edge_index=rep(0,nrow(tree$edge))
tip_pheno=rep(0,length(l))
edge_posterior=rep(0,nrow(tree$edge))
for (i in 1:length(l)) {
  tree$tip.label[i]=l[[i]][1]
  edge_index[which(tree$edge[,2]==i)]=as.numeric(l[[i]][3])
  tip_pheno[i]=as.numeric(l[[i]][5])
  edge_posterior[which(tree$edge[,2]==i)]=as.numeric(l[[i]][7])
}
l=strsplit(tree$node.label,'[\\]\\[\\|=]',perl=T)
node_posterior=rep(0,length(l))
for (i in 1:length(l)) {
  edge_index[which(tree$edge[,2]==i+ntips)]=as.numeric(l[[i]][3])
  edge_posterior[which(tree$edge[,2]==(i+ntips))]=as.numeric(l[[i]][5])
}

#Read rest of file
t=read.table('outfile',comment.char='(')
states =as.matrix(t[,2:(ncol(t)-1)])
lambdas=as.vector(t[,ncol(t)])

#Plot tree showing phenotype and posterior
par(mfrow=c(1,1))
ec=edge_posterior
w=which(tree$edge[,1]==(ntips+1));if (length(w)==2) ec[w]=max(ec[w])
plot.phylo(tree,tip.color=rgb(tip_pheno,0,0),edge.color=rgb(ec,0,0),edge.width=1+ec*10)

#Plot some MCMC traces
par(mfrow=c(1,2))
plot(lambdas,type = 'l',xlab='Sampled iteration',ylab='lambda')
plot(rowSums(states)-1,type='l',xlab='Sampled iteration',ylab='Number of changepoints')

#Plot ten MCMC states
par(mfrow=c(2,5))
for (s in seq(nrow(states)/10,nrow(states),nrow(states)/10)) {
  ec=rep(0,nrow(tree$edge))
  for (i in 1:nrow(tree$edge)) {
    ec[i]=states[s,edge_index[i]+1]
  }
  w=which(tree$edge[,1]==(ntips+1));if (length(w)==2) ec[w]=max(ec[w])
  plot.phylo(tree,tip.color=rgb(tip_pheno,0,0),edge.color=rgb(ec,0,0),edge.width=1+ec*10)
}

#plot the correlation between the change points
image(t(states[,-39]), axes=FALSE)
axis(1, at=seq(0,1,length.out=10), labels= floor(seq(0,ncol(states)-1,length.out=10) ))
axis(2, at=seq(0,1,length.out=10), labels= floor(rev( seq(0,nrow(states),length.out=10)) ))