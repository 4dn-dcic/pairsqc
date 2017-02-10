import pypairix
import math
import os

SEPARATOR = '|'
OUTDIR = 'report'
CIS_TRANS_OUT_FILE = OUTDIR + '/cis_to_trans.out'
PLOT_TABLE_OUT_FILE = OUTDIR + '/plot_table.out'

## pairs ##
#POS1_COL = 2
#POS2_COL = 4
#STRAND1_COL = 5
#STRAND2_COL = 6

## merged_nodup ##
POS1_COL = 2
POS2_COL = 6
STRAND1_COL = 0
STRAND2_COL = 4

## old merged_nodup ##
#POS1_COL = 3
#POS2_COL = 7
#STRAND1_COL = 1
#STRAND2_COL = 5

## pairs
#orientation_list = ['+-','-+','++','--']

## merged_nodup and old merged_nodup
orientation_list = ['016','160','00','1616']

## common across input formats
orientation_names = ['Inner','Outer','Right','Left']


def get_chr_lens ( chromsize_file ):
    chrsize=dict()
    with open(chromsize_file,'r') as f:
        chr, size = f.readline().strip().split('\t')
        chrsize[chr] = int(size)
    return chrsize


def cis_trans_ratio ( pairs_file, DIST_THRES=20000, pos1_col=POS1_COL, pos2_col=POS2_COL, strand1_col=STRAND1_COL, strand2_col=STRAND2_COL):
    """measure cis/trans ratio for a given pairs file"""

    cis=0
    trans=0
     
    tb=pypairix.open( pairs_file )
    chrplist = tb.get_blocknames()
    for chrp in chrplist:
        it = tb.querys2D( chrp )
        chr1, chr2 = chrp.split( SEPARATOR )
        if chr1 == chr2:
            for x in it:
                distance = int( x[pos2_col] ) - int( x[pos1_col] )
    
                # distance will always be > 0 for upper triangle, but in case it is not true.
                if distance > 0:
                    orientation = str( x[strand1_col] ) + str( x[strand2_col] )
                else:
                    orientation = str( x[strand2_col] ) + str( x[strand1_col] )
                    distance = abs(distance)
                if distance > DIST_THRES:
                    cis += 1
        else:
            trans += sum(1 for x in it)
    
    with open(CIS_TRANS_OUT_FILE,'w') as f:
         f.write("Cis reads\t{:,}\n".format(cis))
         f.write("Trans reads\t{:,}\n".format(trans))
         f.write("Cis/Trans ratio\t{:.3f}\n".format(cis/(cis+trans)*100))


def distance_histogram ( pairs_file, chromsize_file, max_logdistance=math.log10(1E5), min_logdistance=math.log10(10), logdistance_binsize=0.1, pos1_col=POS1_COL, pos2_col=POS2_COL, strand1_col=STRAND1_COL, strand2_col=STRAND2_COL, pseudocount=1E-100 ):
    """create a log10-scale binned histogram table for read separation distance histogram
    The histogram is stratefied by read orientation (4 different orientations)
    The table includes raw counts, log10 counts (pseudocounts added), contact probability, log10 contact probability, and proportions for orientation (pseudocounts added)
    Bin is represented by the mid value at the log10 scale. 
    """

    chrsizes = get_chr_lens( chromsize_file )
    genomelen = sum(v for v in chrsizes.values())
    nChr = len(chrsizes)
    max_bin_number = int( ( max_logdistance - logdistance_binsize / 2 ) / logdistance_binsize )

    count = { a: [0] * (max_bin_number + 1) for a in orientation_list }
    log10count = { a: [0] * (max_bin_number + 1) for a in orientation_list }
    pcount = { a: [0] * (max_bin_number + 1) for a in orientation_list }
    sumcount = [0] * (max_bin_number + 1)
    log10sumcount = [0] * (max_bin_number + 1)
    prob = [0] * (max_bin_number + 1)
    log10prob = [0] * (max_bin_number + 1)
    allpossible_sumcount = [0] * (max_bin_number + 1)

    tb=pypairix.open( pairs_file )
    chrplist = tb.get_blocknames()

    # calculate histogram
    for chrp in chrplist:
        chr1, chr2 = chrp.split( SEPARATOR )
        if chr1 == chr2:
            it = tb.querys2D( chrp )
            for x in it:
                distance = int( x[pos2_col] ) - int( x[pos1_col] )

                # distance will always be > 0 for upper triangle, but in case it is not true.
                if distance > 0:
                    orientation = str( x[strand1_col] ) + str( x[strand2_col] )
                else:
                    orientation = str( x[strand2_col] ) + str( x[strand1_col] )
                    distance = abs(distance)
   
                if orientation not in orientation_list: # for some exceptional cases like '4' in merged_nodup
                    continue
 
                # remove zero distance
                if distance > 0:
                    log_distance = math.log10( distance ) 
                    bin_number = int( log_distance / logdistance_binsize ) 
                    if bin_number <= max_bin_number:
                        count[orientation][bin_number] += 1
                        sumcount[bin_number] += 1
 

    # calculate histogram in log10 counts and proportion
    for bin_number in range(0, max_bin_number + 1):
        sc = sumcount[bin_number] + pseudocount * 4
        log10sumcount[bin_number] = math.log10( sc )
        for orientation in orientation_list:
            c = count[orientation][bin_number] + pseudocount
            log10count[orientation][bin_number] = math.log10( c )  
            pcount[orientation][bin_number] = c / sc


    # calculate contact probability
    for bin_number in range(0, max_bin_number + 1):
        bin_x = bin_number * logdistance_binsize + logdistance_binsize/2
        bin_x_size = 10**( bin_x + logdistance_binsize/2 ) - 10**( bin_x - logdistance_binsize/2 )
        allpossible_sumcount[bin_number] = genomelen - nChr * ( 10**bin_x + 1 )  # same as sum of (chrlen - bin_x - 1).
        prob[bin_number] = sumcount[bin_number] / allpossible_sumcount[bin_number] / bin_x_size # normalize by bin size
        log10prob[bin_number] = math.log10( prob[bin_number] + pseudocount )


    # print histogram
    with open(PLOT_TABLE_OUT_FILE,'w') as f:
        header_str = "distance\t" \
            + '\t'.join('count.{}'.format(k) for k in orientation_names) \
            + '\tsum\t' \
            + '\t'.join('log10count.{}'.format(k) for k in orientation_names) \
            + '\tlog10sum\t' \
            + '\t'.join('proportion.{}'.format(k) for k in orientation_names) \
            + '\tallpossible_sumcount' \
            + '\tprob' \
            + '\tlog10prob\n'
        f.write(header_str)
        for bin_number in range(0, max_bin_number + 1):
            bin_x = bin_number * logdistance_binsize + logdistance_binsize/2
            if bin_x <= max_logdistance and bin_x >= min_logdistance:
                print_str = "{:.3f}\t".format(bin_x)
                print_str += '\t'.join('{}'.format(count[ori][bin_number]) for ori in orientation_list )
                print_str += "\t{}\t".format(sumcount[bin_number])
                print_str += '\t'.join('{:.3f}'.format(log10count[ori][bin_number]) for ori in orientation_list )
                print_str += "\t{:.3f}\t".format(log10sumcount[bin_number])
                print_str += '\t'.join('{:.3f}'.format(pcount[ori][bin_number]) for ori in orientation_list )
                print_str += "\t{:.3E}".format(allpossible_sumcount[bin_number])
                print_str += "\t{:.3E}".format(prob[bin_number])
                print_str += "\t{:.3f}\n".format(log10prob[bin_number])
                f.write( print_str )


if __name__ == '__main__':

   import argparse

   parser = argparse.ArgumentParser(description = 'QC for Pairs')
   parser.add_argument('--pairs', help = "input pairs file")
   parser.add_argument('--chrsize', help = "input chromsize file")
   args = parser.parse_args()

   if not os.path.exists(OUTDIR):
       os.mkdir(OUTDIR)

   cis_trans_ratio ( args.pairs )
   distance_histogram ( args.pairs, args.chrsize, max_logdistance = 8.4 , min_logdistance = 1 )


