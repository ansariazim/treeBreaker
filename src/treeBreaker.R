rm(list=ls())
library('ape')

#Read the tree and information attached to it
tree=read.tree('outfile')
l=strsplit(tree$tip.label,'[\\]\\[\\|=]',perl=T)
tip_index=rep(0,length(l))
tip_pheno=rep(0,length(l))
edge_posterior=rep(0,nrow(tree$edge))
for (i in 1:length(l)) {
  tree$tip.label[i]=l[[i]][1]
  tip_index[i]=as.numeric(l[[i]][3])
  tip_pheno[i]=as.numeric(l[[i]][5])
  edge_posterior[which(tree$edge[,2]==i)]=as.numeric(l[[i]][7])
}
l=strsplit(tree$node.label,'[\\]\\[\\|=]',perl=T)
node_index=rep(0,length(l))
node_posterior=rep(0,length(l))
for (i in 1:length(l)) {
  node_index[i]=as.numeric(l[[i]][3])
  edge_posterior[which(tree$edge[,2]==(i+length(tree$tip.label)))]=as.numeric(l[[i]][5])
}

#Plot tree showing phenotype and posterior
plot.phylo(tree,tip.color=rgb(tip_pheno,0,0),edge.color=rgb(edge_posterior,0,0),edge.width=1+edge_posterior*10)