#!/bin/bash
python pairsqc.py -p test_samples/merged_nodup.tab.chrblock_sorted.txt.gz -c test_samples/GRCh37.chrom.sizes.mainonly.female -t M
Rscript plot.r
open report/pairsqc_report.html
zip report.zip report

