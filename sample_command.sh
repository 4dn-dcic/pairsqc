#!/bin/bash
python pairsqc.py --pairs /n/data1/hms/dbmi/park/sl325/4dn/pairs/K562_insitu_merged_nodups.tab.chrblock_sorted.-4.txt.gz --chrsize test_samples/GRCh37.chrom.sizes.mainonly.female -t OM -O K562
#python pairsqc.py --pairs /home/sw229/Storage/HiC_SV/data/HiC_data/GM12878/pairFile/GM12878_insitu_DpnII_merged_nodups.tab.chrblock_sorted.txt.gz --chrsize test_samples/GRCh37.chrom.sizes.mainonly.female -t OM -O GM

