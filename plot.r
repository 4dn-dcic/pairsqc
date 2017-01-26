library(devtools)
install_url("https://github.com/SooLee/Soo.plot.base/archive/0.9.0.zip")
library(Soo.plot.base)

x=read.table("log",skip=3,sep="\t",stringsAsFactors=F,header=T)

# proportion
pngpdf( function(){    
par(plt=c(0.2,0.8,0.2,0.8))
matplot(x[,1],x[,12:15],pch=19,type='o',xlab="distance (log10)", ylab="Proportion",lwd=1,lty=1)
legend("topright",sub('proportion.','',colnames(x)[12:15]),col=1:4,bty="n",pch=19,lwd=1,lty=1)
},"proportion")

# log10 count
pngpdf( function(){
par(plt=c(0.2,0.8,0.2,0.8))
matplot(x[,1],x[,7:10],pch=19,type='o',xlab="distance (log10)", ylab="Contact frequency (log10)",lwd=1,lty=1,ylim=c(0,max(x[,7:10])))
legend("bottomright",sub('log10count.','',colnames(x)[7:10]),col=1:4,bty="n",pch=19,lwd=1,lty=1)
},"log10counts")

entropy<-function(xx)-sum(xx*log2(xx))

pngpdf( function(){
par(plt=c(0.2,0.8,0.2,0.8))
plot(apply(x[,12:15],1,entropy),type='o',pch=19,ylab="entropy",xlab="distance")
},"entropy",height=4)

#> apply(x[,12:15],1,entropy)
# [1] 1.1689760 1.2308699 1.1248270 1.0273758 0.5071788 0.3627795 0.7133532
# [8] 0.3290150 0.2805800 1.4788505 1.8796579 1.8370873 1.8808367 1.9624474
#[15] 1.9959724 1.9993521 1.9995802 1.9986929 1.9979038 1.9985432


pngpdf( function(){
par(plt=c(0.2,0.8,0.2,0.8))
sds=apply(x[,12:15],1,sd)
plot(x[,1], sds,type='o',pch=19,ylab="sd",xlab="distance")
},"sd",height=4)

pngpdf( function(){
par(plt=c(0.2,0.8,0.2,0.8))
sds=apply(x[,12:15],1,sd)
plot(x[,1], sds,type='o',pch=19,ylab="sd",xlab="distance")
cut = min(which(sds<0.02))
abline(v=x[cut ,1],col=2,lty=2)
},"sd_w_cut",height=4)


#sds
# [1] 0.297296709 0.325647048 0.338574832 0.362277058 0.448742688 0.467178321
# [7] 0.420129345 0.470673985 0.476000700 0.267620378 0.115675912 0.133108978
#[13] 0.114702514 0.065324829 0.021556128 0.011814539 0.006976150 0.009287088
#[19] 0.013375973 0.012961481

