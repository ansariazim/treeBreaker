rm(list=ls())
source('simu.R')

#Do the simulations if not already done
if (!file.exists('result.Rdata')) {
library(doParallel)
registerDoParallel(cores=6)

rep=100
codes=matrix(NA,11,rep)
res1=matrix(NA,11,rep)
res2=matrix(NA,11,rep)
for (i in 1:11) {
  out<-foreach(j=1:rep) %dopar% {
    return(testOnSimulation(i-1))
  }
  for (j in 1:length(out)) {res1[i,j]=out[[j]][1];res2[i,j]=out[[j]][2];codes[i,j]=out[[j]][3]}
}
save.image('result.Rdata')
}

#First plotting method
library(beanplot)
pdf('beanplot.pdf',width=7,height=8)
load('result.Rdata')
par(mfrow=c(2,1),mar=c(1.5,1.5,1.5,1.5),oma=c(4,4,0,2),xpd=NA)
beanplot(data.frame(t(res1)),ll=0.05,bw='nrd0',what=c(0,1,1,0),axes=F,ylab='Inferred changing branches')
axis(1,at=1:11,labels=0:10)
axis(2)
bf=res2
bf[bf>1000]=Inf
bf=log10(bf)
bf[is.infinite(bf)]=3#+(runif(length(which(is.infinite(bf))))-0.5)
beanplot(data.frame(t(bf)),ll=0.05,bw='nrd0',what=c(0,1,1,0),axes=F,ylim=c(-1,4),xlab='True number of changing branches',ylab='Bayes Factor (log scale)')
axis(1,at=1:11,labels=0:10)
axis(2,at=c(-1,0,1,2,3),labels=c('0.1','1','10','100','1000+'))
dev.off()

#Second plotting method
library(beeswarm)
load('result.Rdata')
pdf('beeswarm.pdf',width=7,height=8)
par(mfrow=c(2,1),mar=c(1.5,1.5,1.5,1.5),oma=c(4,4,0,2),xpd=NA)
plot(0,0,xlim=c(0,12),ylim=c(0,12),axes=F,xlab='',ylab='Inferred changing branches',col='white')
beeswarm(data.frame(t(res1)),col=rainbow(11),pch=21,corral='wrap',bg = "#00000050",method = 'hex',add=T,cex=0.5)
axis(1,at=1:11,labels=0:10)
axis(2)
bf=res2
bf[bf>1000]=Inf
bf=log10(bf)
bf[is.infinite(bf)]=3+(runif(length(which(is.infinite(bf))))-0.5)*0.5
plot(0,0,xlim=c(0,12),ylim=c(-1,4),axes=F,xlab='True number of changing branches',ylab='Bayes Factor (log scale)',col='white')
beeswarm(data.frame(t(bf)),col=rainbow(11),pch=21,corral='wrap',bg = "#00000050",method = 'hex',add=T,cex=0.5)
#,ll=0.05,bw='nrd0',what=c(0,1,1,0),axes=F,ylim=c(-1,4))
axis(1,at=1:11,labels=0:10)
axis(2,at=c(-1,0,1,2,3),labels=c('0.1','1','10','100','1000+'))
dev.off()
