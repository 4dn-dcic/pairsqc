
args = commandArgs(TRUE)
RE_len = args[1] # either 4 or 6
if(length(args)>1) { report_dir = args[2] } else { report_dir = './report' }

rainbow_w_offset <- function(L, offset = NA){
    if(is.na(offset)) offset = floor(L/2)
    color = rainbow(L+ offset)[1:L]
    return( color )
}

library(plotosaurus)

plot_orientation_proportion_vs_distance <- function(x, RE_len, xlim=c(2,4), no_xlabel=FALSE){
  matplot(x$distance,x[,c('proportion.Inner','proportion.Outer','proportion.Right','proportion.Left')],pch=19,type='o',xlab="", ylab="Proportion",lwd=1,lty=1, col=COLOR4(), axes=F, xlim=xlim)
  exp_axis(xlim,1)
  axis(2)
  box()
  if(no_xlabel==FALSE) mtext("distance",side=1,line=2.5)
  legend("topright",c('Inner','Outer','Right','Left'),col=COLOR4(),bty="n",pch=19,lwd=1,lty=1)
  if(RE_len==4) { abline(v=3.5, lty=2, col="grey") # 4-cutters at 3kb
  } else { abline(v=4.5, lty=2, col="grey") } # 6-cutters at 30kb
} 


plot_orientation_log10count_vs_distance <- function(x, RE_len, xlim=c(2,4), no_xlabel=FALSE){
  ind = which(x$distance>=xlim[1] & x$distance<=xlim[2])
  ylim=range(x[ind,c('log10count.Inner','log10count.Outer','log10count.Right','log10count.Left')])
  matplot(x$distance,x[,c('log10count.Inner','log10count.Outer','log10count.Right','log10count.Left')],pch=19,type='o',xlab="", ylab="Contact frequency",lwd=1,lty=1,ylim=ylim, col=COLOR4(), axes=F, xlim=xlim)
  exp_axis(xlim,1)
  exp_axis(ylim,2)
  box()
  if(no_xlabel==FALSE) mtext("distance",side=1,line=2.5)
  legend("bottomright", c('Inner','Outer','Right','Left'), col=COLOR4(), bty="n", pch=19, lwd=1, lty=1)
  if(RE_len==4) { abline(v=3.5, lty=2, col="grey") # 4-cutters at 3kb
  } else { abline(v=4.5, lty=2, col="grey") } # 6-cutters at 30kb
}


plot_contact_probability_vs_distance <- function(x, xlim=c(3.2,8), no_xlabel=FALSE) {
  ylim=range(x$log10prob)
  plot(x$distance,x$log10prob,pch=19,type='o',xlab="", ylab="Contact probability",lwd=1,lty=1,ylim=ylim,axes=F,xlim=xlim)
  exp_axis(xlim,1)
  exp_axis(ylim,2)
  box()
  if(no_xlabel==FALSE) mtext("distance",side=1,line=2.5)  
  
  abline(v=4, lty=2, col="grey80")
  abline(v=5.5, lty=2, col="grey80")

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
  par(las=1,lend=2)
  par(plt=plt)
  entropies=apply(x[,c('proportion.Inner','proportion.Outer','proportion.Right','proportion.Left')],1,entropy)
  plot(x$distance, entropies,type='o',pch=19,ylab="entropy",xlab="distance",xlim=xlim,axes=F)
  exp_axis(xlim,1)
  axis(2)
}


plot_sd_with_cutoff <- function(x, xlim=c(2,4)){
  sds=apply(x[,c('proportion.Inner','proportion.Outer','proportion.Right','proportion.Left')],1,sd)
  plot(x$distance, log10(sds), type='o',pch=19,ylab="sd",xlab="distance",xlim=xlim,axes=F, col=COLOR3B()[1])
  exp_axis(xlim,1)
  exp_axis(log10(sds), 2, 3)
  box()
  almost_converged = which(sds<0.02)
  well_converged = which(sds<0.005)
  points(x$distance[almost_converged], log10(sds[almost_converged]), col=COLOR3B()[2], pch=19)
  points(x$distance[well_converged], log10(sds[well_converged]), col=COLOR3B()[3], pch=19)
  abline(v=x[almost_converged[1] ,"distance"],col=COLOR3B()[2],lty=2)
  abline(v=x[well_converged[1] ,"distance"],col=COLOR3B()[3],lty=2)
  legend("topright",c("sd<0.02","sd<0.005"),col=COLOR3B()[2:3],pch=19,bty='n')
}


plot_contact_frequency_vs_genomic_separation_per_chr <- function(x) { 
  y=data.frame(x[,grep('log10prob_per_chr',colnames(x))])   
  valid = which(x$distance > 3.5 & x$distance < 7)
  y=y[valid,]
  ylim = range(apply(y,2,function(xx){ ind =which(xx>-90); if(length(ind)==0) return(NA) else return( range(xx[ind])) }), na.rm=T)
  L = ncol(y)
  chrnames = sub("log10prob_per_chr.","",colnames(y))
  
  plt0 = par('plt', no.readonly=FALSE)
  plt1 = plt0; plt2= plt0
  pltx.plotend = (plt0[2] - plt0[1]) * 7/8 + plt0[1]
  pltx.legendstart = (plt0[2] - plt0[1]) * 29/32 + plt0[1]
  plt1[2] = pltx.plotend
  plt2[1] = pltx.legendstart
  par(plt=plt1)
  matplot(x$distance[valid],y,ylim=ylim, pch=19, type='l', lty=1, lwd=1, col=rainbow_w_offset(L), xlab="genomic separation", ylab="Contact frequency", axes=F)
  exp_axis(x$distance[valid], 1)
  exp_axis(ylim, 2)
  box()
  abline(v=4, lty=2, col="grey80")
  abline(v=5.5, lty=2, col="grey80")
  par(new=T, plt=plt2)
  plot.new(); plot.window(c(0,1),c(0,1), xaxs='i', yaxs='i')
  legend("topleft", chrnames, col = rainbow_w_offset(L), lty=1, lwd=1, bty='n' ,xpd=T, cex=0.8)
}


########################


require( "Nozzle.R1" )

generate_pairsqc_report <- function ( sample_name = NA) {

  if(is.na(sample_name)) { 
    report_title = "PairsQC Report"
  } else {
    report_title = paste("PairsQC Report for sample", sample_name, sep=" ") 
  }

  # Phase 1: create report elements
  r <- newCustomReport( report_title ); 
  s1 <- newSection( "Cis-to-trans ratio" );
  s2 <- newSection( "Proportion of read orientation versus genomic separation");
  s3 <- newSection( "Number of reads versus genomic separation, stratified by read orientation" );
  s4 <- newSection( "Contact probability versus genomic separation" );
  s5 <- newSection( "Contact probability versus genomic separation, per chromosome" );
  references <- newSection( "References" );

  # References
  refRao = newCitation( authors= 'Rao et al.', title= 'A 3D Map of the Human Genome at Kilobase Resolution Reveals Principles of Chromatin Looping', publication= 'Cell', issue= '159:7', pages='1665-1680', year='2014' )
  refSanborn = newCitation( authors= 'Sanborn et al.', title= 'Chromatin extrusion explains key features of loop and domain formation in wild-type and engineered genome', publication= 'Proc. Natl. Acad. Sci. USA', issue= '112:47', pages='E6456-E6465', year='2015' )
  refImakaev = newCitation( authors= 'Imakaev et al.', title= 'Iterative correction of Hi-C data reveals hallmarks of chromosome organization', publication= 'Nature Methods', issue= '9:10', pages='999-1003', year='2012' )

  # section 1
  table_files = list.files('.', pattern='*cis_to_trans.out$')
  sample_names = sub('.?cis_to_trans.out','',table_files, perl=T, fixed=F)
  tableData1 = do.call(data.frame, sapply(table_files, function(table_file){
    tableData=read.table(table_file,stringsAsFactors=F, header=F, sep="\t")
    colnames(tableData) = c('QC field','value')
    return(tableData)
  }, simplify=F))[,c(1,seq(2,2*length(table_files),2))]
  colnames(tableData1) = c('QC field', sample_names)

  t1 <- newTable( tableData1, significantDigits=0, exportId="TABLE_1", "Cis-to-trans ratio" );
  p1 <- newParagraph( "Cis-to-trans ratio is computed as the ratio of long-range cis reads (>20kb) to trans reads plus long-range cis reads. Typically, cis-to-trans ratio higher than 40% is required. Percentage of long-range cis reads is the ratio of long-range cis reads to total number of reads. Minimum 15% is required and 40% or higher suggests a good library", asReference( refRao ), "." );
  
  # section 2
  # using sample_names from section 1 (same ordering)
  png2_files = paste('plots/', sample_names, '.proportion.png',sep='')
  pdf2_files = gsub('png$','pdf',png2_files)
  f2list <- mapply(function(png2, pdf2, sample_name) { 
                       newFigure( png2, fileHighRes=pdf2, exportId=paste("FIGURE_2", sample_name, sep="."), 
                                  paste("Proportion of read orientation versus genomic separation : ", asStrong(sample_name), sep=""));
                   }, png2_files, pdf2_files, sample_names, SIMPLIFY=FALSE)
  eq2a = '10^3.5'
  eq2b = '10^4.5'
  p2 <- newParagraph("Contact frequency (number of reads, left) and proportion of reads (right) are shown, stratified by read orientation. Good four-cutter and six-cutter samples would converge at ~3kb (", eq2a, ") and ~30kb (", eq2b, "), respectively", asReference( refRao ), ". Convergence is determined by the standard deviation of the proportions being less than 0.005.")
  
  # section 3
  #png3= 'plots/log10counts.png'
  #pdf3= 'plots/log10counts.pdf'
  #f3 <- newFigure( png3, fileHighRes=pdf3, exportId="FIGURE_3", "Number of reads versus genomic separation, stratified by read orientation"); 
  
  # section 4
  # using sample_names from section 1 (same ordering)
  png4_files = paste('plots/', sample_names, '.log10prob.png',sep='')
  pdf4_files = gsub('png$','pdf',png4_files)
  f4list <- mapply(function(png4, pdf4, sample_name) {
                       newFigure( png4, fileHighRes=pdf4, 
                                  paste("Contact probability versus genomic separation : ", asStrong(sample_name), sep=""));
                   }, png4_files, pdf4_files, sample_names, SIMPLIFY=FALSE)
  p4 <- newParagraph( "Contact probability (number of reads, normalized by number of bins and bin size) is shown with respect to genomic separation between mates", asReference( refImakaev), asReference( refSanborn ), ". The slope between distance 10kb ~ 300kb (10^4 ~ 10^5.5) representing a TAD is calculated. A good mitotic sample would have a slope close to ~ -0.76", asReference( refSanborn ), "." );
  
  # section 5
  # png5= 'plots/log10prob_chr.png'
  # pdf5= 'plots/log10prob_chr.pdf'
  # f5 <- newFigure( png5, fileHighRes=pdf5, exportId="FIGURE_5", "Contact probability versus genomic separation, separated by chromosome");

  # Phase 2: assemble report structure bottom-up
  s1 <- addTo( s1, p1, t1);
  s2 <- addTo( s2, p2);
  for(f2 in f2list) { s2 <- addTo( s2, f2); }
  #s3 <- addTo( s3, f3);
  s4 <- addTo( s4, p4 );
  for(f4 in f4list) { s4 <- addTo( s4, f4); }
  # s5 <- addTo( s5, f5);
  references <- addTo( references, refRao, refImakaev, refSanborn )
  #r <- addTo( r, s1, s2, s3, s4, references );
  #r <- addTo( r, s1, s2, s4, s5, references );
  r <- addTo( r, s1, s2, s4, references );
  
  # Phase 3: render report to file
  writeReport( r, filename="pairsqc_report" ); # w/o extension

}


##################


cwd = getwd()
sample_names = gsub('.?plot_table.out$', '', list.files(report_dir, pattern='*.plot_table.out$'), perl=T, fixed=F)

plot_for_sample<-function(sample_name, report_dir) {
  setwd(cwd)
  plot_dir = paste(report_dir,"plots",sep="/")
  plot_table_file = paste(report_dir, "/", sample_name, ".plot_table.out", sep="")
  
  x=read.table(plot_table_file,sep="\t",stringsAsFactors=F,header=T)
  dir.create(plot_dir, showWarnings = FALSE, recursive = TRUE)
  setwd(plot_dir)
  pngpdf_preset( list(
      function() plot_orientation_log10count_vs_distance(x, RE_len),
      function() plot_orientation_proportion_vs_distance(x, RE_len, no_xlabel=T),
      function() plot_sd_with_cutoff(x)
    ), paste(sample_name, "proportion", sep="."), stylefunc0, get_preset(3), add.date=FALSE
  )
  
  pngpdf_preset( list(
      function() plot_contact_probability_vs_distance(x),
      function() plot_contact_frequency_vs_genomic_separation_per_chr(x)
    ), paste(sample_name, "log10prob", sep="."), stylefunc0, get_preset(2,'h50'), add.date=FALSE
  )
  
  #pngpdf.nodate( function()plot_contact_probability_vs_distance(x) ,"log10prob")
  #pngpdf.nodate( function()plot_contact_frequency_vs_genomic_separation_per_chr(x), "log10prob_chr", height=5.5 )
  #pngpdf.nodate( function()plot_entropy(x), "entropy", height=4) 
}

sapply(sample_names, plot_for_sample, report_dir = report_dir)

setwd(cwd)
setwd(report_dir)
generate_pairsqc_report()

