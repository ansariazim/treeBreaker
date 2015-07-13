rm(list=ls())
library('ape')

#Read the tree and information attached to it
tree=read.tree('outfile')
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
  edge_index[which(tree$edge[,2]==i+length(tree$tip.label))]=as.numeric(l[[i]][3])
  edge_posterior[which(tree$edge[,2]==(i+length(tree$tip.label)))]=as.numeric(l[[i]][5])
}

#Read rest of file
t=read.table('outfile',comment.char='(')
states =as.matrix(t[,2:(ncol(t)-1)])
lambdas=as.vector(t[,ncol(t)])

#Plot tree showing phenotype and posterior
par(mfrow=c(1,1))
plot.phylo(tree,tip.color=rgb(tip_pheno,0,0),edge.color=rgb(edge_posterior,0,0),edge.width=1+edge_posterior*10)

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
plot.phylo(tree,tip.color=rgb(tip_pheno,0,0),edge.color=rgb(ec,0,0),edge.width=1+ec*10)
}