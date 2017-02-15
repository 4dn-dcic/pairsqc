#!/bin/bash
cp /home/sw229/Storage/HiC_SV/data/HiC_data/GM12878/pairFile/GM12878_insitu_DpnII_merged_nodups.tab.chrblock_sorted.txt $PAIRS/
sort -k3,3 -k7,7 -k4,4 -k8,8 /n/data1/hms/dbmi/park/sl325/4dn/pairs/GM12878_insitu_DpnII_merged_nodups.tab.chrblock_sorted.txt > sort -k3,3 -k7,7 -k4,4 -k8,8 /n/data1/hms/dbmi/park/sl325/4dn/pairs/GM12878_insitu_DpnII_merged_nodups.tab.bsorted.txt
bgzip -f /n/data1/hms/dbmi/park/sl325/4dn/pairs/GM12878_insitu_DpnII_merged_nodups.tab.bsorted.txt
pairix -f -s3 -d7 -b4 -e4 -u8 -T /n/data1/hms/dbmi/park/sl325/4dn/pairs/GM12878_insitu_DpnII_merged_nodups.tab.bsorted.txt.gz

