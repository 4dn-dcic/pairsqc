
# use pairsqc
python pairsqc --pairs ../pairix/samples/merged_nodup.tab.chrblock_sorted.txt.gz > testout1

# use awk & bash
cis=`gunzip -c ../pairix/samples/merged_nodup.tab.chrblock_sorted.txt.gz | awk '$2==$6' |wc -l`
trans=`gunzip -c ../pairix/samples/merged_nodup.tab.chrblock_sorted.txt.gz | awk '$2!=$6' |wc -l`
ratio=`expr $cis / $trans`

echo "Cis reads	$cis" > testout2
echo "Trans reads	$trans" >> testout2
echo -n "Cis/Trans ratio	" >> testout2
echo "scale=2; $cis/$trans" |bc -l >> testout2

diff testout1 testout2


