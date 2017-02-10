#library(devtools)
#install_github("parklab/nozzle")

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
  
  tableData1=read.table('cis_to_trans.out',stringsAsFactors=F, header=F, sep="\t")
  colnames(tableData1) = c('QC field','value')
  t1 <- newTable( tableData1, significantDigits=0, exportId="TABLE_1", "Cis-to-trans ratio" );
  
  png2= 'plots/proportion.png'
  pdf2= 'plots/proportion.pdf'
  f2 <- newFigure( png2, fileHighRes=pdf2, exportId="FIGURE_2", "Proportion of read orientation versus genomic separation");
  
  png3= 'plots/log10counts.png'
  pdf3= 'plots/log10counts.pdf'
  f3 <- newFigure( png3, fileHighRes=pdf3, exportId="FIGURE_3", "Number of reads versus genomic separation, stratified by read orientation"); 
  
  png4= 'plots/log10prob.png'
  pdf4= 'plots/log10prob.pdf'
  f4 <- newFigure( png4, fileHighRes=pdf4, exportId="FIGURE_4", "Contact probability versus genomic separation");
  
  # Phase 2: assemble report structure bottom-up
  s1 <- addTo( s1, t1);
  s2 <- addTo( s2, f2);
  s3 <- addTo( s3, f3);
  s4 <- addTo( s4, f4);
  r <- addTo( r, s1, s2, s3, s4 );
  
  # Phase 3: render report to file
  writeReport( r, filename="pairsqc_report" ); # w/o extension

}

