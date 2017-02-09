source("plot_functions.r")
source("nozzle_report.r")

x=read.table("report/plot_table.out",skip=0,sep="\t",stringsAsFactors=F,header=T)

dir.create("report/plots", showWarnings = FALSE, recursive = TRUE)
setwd("report/plots")
pngpdf.nodate( function(){
  plot_orientation_proportion_vs_distance(x, plt=c(0.2,0.8,0.45,0.9), no_xlabel=T)
  par(new=T)
  plot_sd_with_cutoff(x, plt=c(0.2,0.8,0.2,0.35))
}, "proportion" )

pngpdf.nodate( function()plot_orientation_log10count_vs_distance(x) , "log10counts" )
pngpdf.nodate( function()plot_contact_probability_vs_distance(x) ,"log10prob")
pngpdf.nodate( function()plot_entropy(x), "entropy", height=4)
#pngpdf.nodate( function()plot_sd_with_cutoff(x),"sd_w_cut", height=4)

setwd("..")
generate_pairsqc_report()
setwd("..")

