library(devtools)
install_url("https://github.com/SooLee/Soo.plot.base/archive/0.9.0.zip")
library(Soo.plot.base)

exp_axis<-function(x, axis_ind){
  at=pretty(x)
  at_exp = sapply(at, function(xx)as.expression(bquote(10^ .(xx))))
  axis(axis_ind, at=at, label=at_exp)
}

pngpdf.nodate <- function(f,...) pngpdf(f,...,add.date=FALSE)

plot_orientation_proportion_vs_distance <- function(x, xlim=c(2,4), plt=c(0.2,0.95,0.2,0.95), no_xlabel=FALSE){
  par(plt=plt)
  matplot(x$distance,x[,12:15],pch=19,type='o',xlab="", ylab="Proportion",lwd=1,lty=1, axes=F, xlim=xlim)
  exp_axis(xlim,1)
  axis(2)
  box()
  if(no_xlabel==FALSE) mtext("distance",side=1,line=2.5)
  legend("topright",sub('proportion.','',colnames(x)[12:15]),col=1:4,bty="n",pch=19,lwd=1,lty=1)
} 


plot_orientation_log10count_vs_distance <- function(x, xlim=c(2,4), plt=c(0.2,0.95,0.2,0.95), no_xlabel=FALSE){
  par(plt=plt)
  ind = which(x$distance>=xlim[1] & x$distance<=xlim[2])
  ylim=range(x[ind,7:10])
  matplot(x$distance,x[,7:10],pch=19,type='o',xlab="", ylab="Contact frequency",lwd=1,lty=1,ylim=ylim, axes=F, xlim=xlim)
  exp_axis(xlim,1)
  exp_axis(ylim,2)
  box()
  if(no_xlabel==FALSE) mtext("distance",side=1,line=2.5)
  legend("bottomright",sub('log10count.','',colnames(x)[7:10]),col=1:4,bty="n",pch=19,lwd=1,lty=1)
}


plot_contact_probability_vs_distance <- function(x, xlim=c(3.2,8), plt=c(0.2,0.95,0.2,0.95), no_xlabel=FALSE) {
  par(plt=plt)
  ylim=range(x[,18])
  plot(x$distance,x[,18],pch=19,type='o',xlab="", ylab="Contact probability",lwd=1,lty=1,ylim=ylim,axes=F,xlim=xlim)
  exp_axis(xlim,1)
  exp_axis(ylim,2)
  box()
  if(no_xlabel==FALSE) mtext("distance",side=1,line=2.5)  

  # slope
  tad=which(x$distance>=4 & x$distance<=5.5)
  xx=x$distance[tad]
  yy=x$log10prob[tad]
  tad_slope=glm(yy~xx)$coefficients[2]
  tad_intercept=glm(yy~xx)$coefficients[1]
  spacing=abs(yy[1]*0.07)
  text(xx[1], yy[1]+spacing, paste("slope=",round(tad_slope,2),sep=' '), pos=4)
  abline(a=tad_intercept + spacing, b=tad_slope, lty=2)
}


# entropy
entropy<-function(xx)-sum(xx*log2(xx))
plot_entropy <- function(x, xlim=c(2,4), plt=c(0.2,0.95,0.2,0.95)){
  par(plt=plt)
  entropies=apply(x[,12:15],1,entropy)
  plot(entropies,type='o',pch=19,ylab="entropy",xlab="distance",xlim=xlim,axes=F)
  exp_axis(x$distance,1)
  axis(2)
}


plot_sd_with_cutoff <- function(x, xlim=c(2,4), plt=c(0.2,0.95,0.2,0.95)){
  par(plt=plt)
  sds=apply(x[,12:15],1,sd)
  plot(x$distance, sds,type='o',pch=19,ylab="sd",xlab="distance",xlim=xlim,axes=F)
  exp_axis(xlim,1)
  axis(2)
  box()
  cut = min(which(sds<0.02))
  abline(v=x[cut ,"distance"],col=2,lty=2)
}


