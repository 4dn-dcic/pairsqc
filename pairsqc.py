import pypairix
import math
import os

SEPARATOR = '|'
OUTDIR = 'report'
CIS_TRANS_OUT_FILE = OUTDIR + '/cis_to_trans.out'
PLOT_TABLE_OUT_FILE = OUTDIR + '/plot_table.out'

# I'm testing this using old merged nodup but later we will switch to pairs.

class ColIndices(object):
    """Column indices for position1, position2, strand1 and strand2, 0-based"""

    def __init__(self, pos1, pos2, strand1, strand2):
        self.pos1 = pos1
        self.pos2 = pos2
        self.strand1 = strand1
        self.strand2 = strand2


## column indices per file type
cols_pairs = ColIndices(2, 4, 5, 6)
cols_merged_nodups = ColIndices(2, 6, 0, 4)
cols_old_merged_nodups = ColIndices(3, 7, 1, 5)

## orientation representation per file type
orientation_list_pairs = ['+-','-+','++','--']
orientation_list_merged_nodups = ['016','160','00','1616']

## common across input formats
orientation_names = ['Inner','Outer','Right','Left']


class GenomeSize(object):

    def __init__(self, chromsize_file):
        """return a dictionary of chromosome : size pairs from a chromsize file."""
        self.chrsize=dict()
        self.total_len=0
        with open(chromsize_file,'r') as f:
            chr, size = f.readline().strip().split('\t')
            self.chrsize[chr] = int(size)
            self.total_len += int(size)
        self.nChr = len(self.chrsize)


def get_distance_and_orientation (line, cols):
    """return distance and orientation 
    given a list representing a line from the pairs input file
    and a ColIndices object
    """
    distance = int(line[cols.pos2]) - int(line[cols.pos1])
    
    # distance will always be > 0 for upper triangle, but in case it is not true.
    if distance > 0:
        orientation = str(line[cols.strand1]) + str(line[cols.strand2])
    else:
        orientation = str(line[cols.strand2]) + str(line[cols.strand1])
        distance = abs(distance)

    return(distance, orientation)


def cis_trans_ratio (pairs_file, DIST_THRES=20000, cols= cols_pairs):
    """measure cis/trans ratio for a given pairs file"""

    cis=0
    trans=0
    cis_short=0
     
    tb=pypairix.open(pairs_file )
    chrplist = tb.get_blocknames()
    for chrp in chrplist:
        it = tb.querys2D(chrp )
        chr1, chr2 = chrp.split(SEPARATOR )
        if chr1 == chr2:
            for x in it:
                distance, orientation = get_distance_and_orientation (x, cols)
                if distance > DIST_THRES:
                    cis += 1
                else:
                    cis_short += 1
        else:
            trans += sum(1 for x in it)

    total = cis + cis_short + trans
    
    with open(CIS_TRANS_OUT_FILE,'w') as f:
         f.write("Total reads\t{:,}\n".format(total))
         f.write("Short cis reads (<20kb)\t{:,}\n".format(cis_short))
         f.write("Cis reads (>20kb)\t{:,}\n".format(cis))
         f.write("Trans reads\t{:,}\n".format(trans))
         f.write("Cis/Trans ratio\t{:.3f}\n".format(cis/(cis+trans)*100))


class SeparationStat(object):
    """Statistics to be calculated for each separation distance bin"""

    def __init__(self, orientation_list, pseudocount=1E-100):
        self.orientation_list = orientation_list
        self.count = { a: 0 for a in orientation_list }
        self.log10count = { a: 0 for a in orientation_list }
        self.pcount = { a: 0 for a in orientation_list }
        self.sumcount = 0
        self.log10sumcount = 0
        self.prob = 0
        self.log10prob = 0
        self.allpossible_sumcount = 0
        self.pseudocount = pseudocount

    def increment(self, orientation):
        self.count[orientation] += 1
        self.sumcount += 1

    def calculate_log10count(self):
        for orientation in self.orientation_list:
            c = self.count[orientation] + self.pseudocount
            self.log10count[orientation] = math.log10( c )  

    def calculate_log10sumcount(self):
        sc = self.sumcount + self.pseudocount * 4
        self.log10sumcount = math.log10( sc )

    def calculate_pcount(self):
        sc = self.sumcount + self.pseudocount * 4
        for orientation in self.orientation_list:
            c = self.count[orientation] + self.pseudocount
            self.pcount[orientation] = c / sc

    def calculate_contact_probability(self, s, bin_size, gs):
        """Calculate contact probability for a given separation distance and bin size
        s is the representative log10 separation distance.
        gs: GenomeSize object
        """
        self.allpossible_sumcount = gs.total_len - gs.nChr * ( 10**s + 1 )  # same as sum of (chrlen - s - 1).
        self.prob = self.sumcount / self.allpossible_sumcount / bin_size 
        self.log10prob = math.log10(self.prob + self.pseudocount)


    def print_content(self, fout, bin_mid):
        print_str = "{:.3f}\t".format(bin_mid)
        print_str += '\t'.join('{}'.format(self.count[ori]) for ori in self.orientation_list )
        print_str += "\t{}\t".format(self.sumcount)
        print_str += '\t'.join('{:.3f}'.format(self.log10count[ori]) for ori in self.orientation_list )
        print_str += "\t{:.3f}\t".format(self.log10sumcount)
        print_str += '\t'.join('{:.3f}'.format(self.pcount[ori]) for ori in self.orientation_list )
        print_str += "\t{:.3E}".format(self.allpossible_sumcount)
        print_str += "\t{:.3E}".format(self.prob)
        print_str += "\t{:.3f}\n".format(self.log10prob)
        fout.write(print_str)


    def print_header(fout):
        header_str = "distance\t" \
            + '\t'.join('count.{}'.format(k) for k in orientation_names) \
            + '\tsum\t' \
            + '\t'.join('log10count.{}'.format(k) for k in orientation_names) \
            + '\tlog10sum\t' \
            + '\t'.join('proportion.{}'.format(k) for k in orientation_names) \
            + '\tallpossible_sumcount' \
            + '\tprob' \
            + '\tlog10prob\n'
        fout.write(header_str)


def distance_histogram (pairs_file, chromsize_file, cols=cols_pairs, orientation_list = orientation_list_pairs, max_logdistance=math.log10(1E5), min_logdistance=math.log10(10), log_binsize=0.1):
    """create a log10-scale binned histogram table for read separation distance histogram
    The histogram is stratefied by read orientation (4 different orientations)
    The table includes raw counts, log10 counts (pseudocounts added), contact probability, log10 contact probability, and proportions for orientation (pseudocounts added)
    Bin is represented by the mid value at the log10 scale. 
    log_binsize: distance bin size in log10 scale.
    """

    gs = GenomeSize(chromsize_file)
    max_bin_number = int( ( max_logdistance - log_binsize / 2 ) / log_binsize )

    ss = []
    for a in range(0, max_bin_number+1):
        ss.append(SeparationStat(orientation_list))

    tb=pypairix.open( pairs_file )
    chrplist = tb.get_blocknames()

    # calculate histogram
    for chrp in chrplist:
        chr1, chr2 = chrp.split( SEPARATOR )
        if chr1 == chr2:
            it = tb.querys2D( chrp )
            for x in it:
                distance, orientation = get_distance_and_orientation (x, cols)
                if orientation not in orientation_list: # for some exceptional cases like '4' in merged_nodup
                    continue
 
                # remove zero distance, count.
                if distance > 0:
                    log_distance = math.log10( distance ) 
                    bin_number = int( log_distance / log_binsize ) 
                    if bin_number <= max_bin_number:
                        ss[bin_number].increment(orientation)

    # calculate histogram in log10 counts and proportion
    for bin_number in range(0, max_bin_number + 1):
        ss[bin_number].calculate_log10count()
        ss[bin_number].calculate_log10sumcount()
        ss[bin_number].calculate_pcount()

    # calculate contact probability
    for bin_number in range(0, max_bin_number + 1):
        bin_mid = bin_number * log_binsize + log_binsize/2
        bin_size = 10**( bin_mid + log_binsize/2 ) - 10**( bin_mid - log_binsize/2 )
        ss[bin_number].calculate_contact_probability(bin_mid, bin_size, gs)

    # print histogram
    with open(PLOT_TABLE_OUT_FILE,'w') as f:
        SeparationStat.print_header(f)
        for bin_number in range(0, max_bin_number + 1):
            bin_mid = bin_number * log_binsize + log_binsize/2
            if bin_mid <= max_logdistance and bin_mid >= min_logdistance:
                ss[bin_number].print_content(f, bin_mid)


if __name__ == '__main__':

   import argparse

   parser = argparse.ArgumentParser(description = 'QC for Pairs')
   parser.add_argument('--pairs', help = "input pairs file")
   parser.add_argument('--chrsize', help = "input chromsize file")
   args = parser.parse_args()

   if not os.path.exists(OUTDIR):
       os.mkdir(OUTDIR)

   cis_trans_ratio (args.pairs, cols = cols_merged_nodups )
   distance_histogram (args.pairs, args.chrsize, max_logdistance = 8.4 , min_logdistance = 1, cols = cols_merged_nodups, orientation_list = orientation_list_merged_nodups)


