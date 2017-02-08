library(devtools)
install_url("https://github.com/SooLee/Soo.plot.base/archive/0.9.0.zip")
library(Soo.plot.base)

x=read.table("out",skip=0,sep="\t",stringsAsFactors=F,header=T)
x=x[1:70,] # up to 10^8.

setwd("plots")

exp_axis<-function(x,axis_ind){
  at=pretty(x)
  at_exp = sapply(at, function(xx)as.expression(bquote(10^ .(xx))))
  axis(axis_ind, at=at, label=at_exp)
}


# proportion
pngpdf( function(){    
par(plt=c(0.2,0.8,0.2,0.8))
xlim=c(2,4)
matplot(x[,1],x[,12:15],pch=19,type='o',xlab="distance (log10)", ylab="Proportion",lwd=1,lty=1, axes=F, xlim=xlim)
exp_axis(x[,1],1)
axis(2)
box()
legend("topright",sub('proportion.','',colnames(x)[12:15]),col=1:4,bty="n",pch=19,lwd=1,lty=1)
},"proportion")

# log10 count
pngpdf( function(){
par(plt=c(0.2,0.8,0.2,0.8))
ylim=c(0,max(x[,7:10]))
xlim=c(2,4)
matplot(x[,1],x[,7:10],pch=19,type='o',xlab="distance (log10)", ylab="Contact frequency (log10)",lwd=1,lty=1,ylim=ylim, axes=F, xlim=xlim)
exp_axis(xlim,1)
exp_axis(ylim,2)
box()
legend("bottomright",sub('log10count.','',colnames(x)[7:10]),col=1:4,bty="n",pch=19,lwd=1,lty=1)
},"log10counts")

# log10 prob
pngpdf( function(){
par(plt=c(0.2,0.8,0.2,0.8))
ylim=range(x[,18])
xlim=c(3.2,8)
plot(x[,1],x[,18],pch=19,type='o',xlab="distance (log10)", ylab="Contact probability (log10)",lwd=1,lty=1,ylim=ylim,axes=F,xlim=xlim)
exp_axis(xlim,1)
exp_axis(ylim,2)
box()

# slope
tad=which(x$distance>=4 & x$distance<=5.5)
xx=x$distance[tad]
yy=x$log10prob[tad]
tad_slope=glm(yy~xx)$coefficients[2]
tad_intercept=glm(yy~xx)$coefficients[1]
spacing=abs(yy[1]*0.07)
text(xx[1], yy[1]+spacing, paste("slope=",round(tad_slope,2),sep=' '), pos=4)
abline(a=tad_intercept + spacing, b=tad_slope, lty=2)

}, "log10prob")

entropy<-function(xx)-sum(xx*log2(xx))

# entropy
pngpdf( function(){
par(plt=c(0.2,0.8,0.2,0.8))
entropies=apply(x[,12:15],1,entropy)
plot(entropies,type='o',pch=19,ylab="entropy",xlab="distance")
print(entropies)
},"entropy",height=4)


# sd with cut off
pngpdf( function(){
par(plt=c(0.2,0.8,0.2,0.8))
sds=apply(x[,12:15],1,sd)
xlim=c(2,4)
plot(x[,1], sds,type='o',pch=19,ylab="sd",xlab="distance",xlim=xlim)
cut = min(which(sds<0.02))
abline(v=x[cut ,1],col=2,lty=2)
print(sds)
},"sd_w_cut",height=4)


