import pypairix
import math
import os

SEPARATOR = '|'
CIS_TRANS_OUT_FILE_NAME =  '/cis_to_trans.out'
PLOT_TABLE_OUT_FILE_NAME = '/plot_table.out'


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
            for line in f:
                chr, size = line.strip().split('\t')
                self.chrsize[chr] = int(size)
                self.total_len += int(size)
        self.nChr = len(self.chrsize)


class CisTransStat(object):
    """Summary statistics including cis-trans ratio"""

    def __init__(self):
        self.cis = 0
        self.trans = 0
        self.cis_short = 0
        self.total = 0

    def calculate_total(self):
        self.total = self.cis + self.cis_short + self.trans

    def print_stat(self, fout):
        fout.write("Total reads\t{:,}\n".format(self.total))
        fout.write("Short cis reads (<20kb)\t{:,}\n".format(self.cis_short))
        fout.write("Cis reads (>20kb)\t{:,}\n".format(self.cis))
        fout.write("Trans reads\t{:,}\n".format(self.trans))
        fout.write("Cis/Trans ratio\t{:.3f}\n".format(float(self.cis)/float(self.cis+self.trans)*100))


class SeparationStat(object):
    """Statistics to be calculated for each separation distance bin"""

    def __init__(self, orientation_list, gs, pseudocount=1E-100):
        """gs: GenomeSize object"""
        self.orientation_list = orientation_list
        self.gs = gs
        self.chr_list = list(gs.chrsize.keys())
        self.chr_list.sort()
        self.pseudocount = pseudocount

        # per-orientation
        self.count_per_ori = { a: 0 for a in orientation_list }
        self.log10count_per_ori = { a: 0 for a in orientation_list }
        self.pcount_per_ori = { a: 0 for a in orientation_list }

        # per-chromosome
        self.count_per_chr = { a: 0 for a in gs.chrsize.keys() }
        self.allpossible_count_per_chr = { a: 0 for a in gs.chrsize.keys() }
        self.prob_per_chr = { a: 0 for a in gs.chrsize.keys() }
        self.log10prob_per_chr = { a: 0 for a in gs.chrsize.keys() }

        # total
        self.sumcount = 0
        self.log10sumcount = 0
        self.prob = 0
        self.log10prob = 0
        self.allpossible_sumcount = 0

    def increment(self, orientation, chr):
        """increment both count_per_ori and count_per_chr together, so that we don't count the read on a weird chromosome for orientation and vice versa"""
        if orientation in self.orientation_list: # skip if not included in orientation list
            if chr in self.chr_list: # skip if not included in chr list
                self.count_per_ori[orientation] += 1
                self.count_per_chr[chr] += 1

    def calculate_sumcount(self):
        self.sumcount = sum(self.count_per_ori.values())
        assert self.sumcount == sum(self.count_per_chr.values())

    def calculate_log10count_per_ori(self):
        for orientation in self.orientation_list:
            c = self.count_per_ori[orientation] + self.pseudocount
            self.log10count_per_ori[orientation] = math.log10( c )  

    def calculate_log10sumcount(self):
        sc = self.sumcount + self.pseudocount * 4
        self.log10sumcount = math.log10( sc )

    def calculate_pcount_per_ori(self):
        sc = self.sumcount + self.pseudocount * 4
        for orientation in self.orientation_list:
            c = self.count_per_ori[orientation] + self.pseudocount
            self.pcount_per_ori[orientation] = c / sc

    def calculate_contact_probability_per_chr(self, s, bin_size):
        """Calculate contact probability for a given separation distance and bin size
        s is the representative log10 separation distance.
        """
        for chr in self.chr_list:
            self.allpossible_count_per_chr[chr] = self.gs.chrsize[chr] - 10**s - 1
            if self.allpossible_count_per_chr[chr] <= 0: # the chromosome is smaller than s
                self.allpossible_count_per_chr[chr] = 0
                self.prob_per_chr[chr] = 0
            else:
                self.prob_per_chr[chr] = self.count_per_chr[chr] / self.allpossible_count_per_chr[chr] / bin_size 
                self.log10prob_per_chr[chr] = math.log10(self.prob_per_chr[chr] + self.pseudocount)

    def calculate_contact_probability(self, s, bin_size):
        """Calculate contact probability for a given separation distance and bin size
        s is the representative log10 separation distance.
        """
        self.allpossible_sumcount = sum(self.allpossible_count_per_chr.values()) 
        self.prob = self.sumcount / self.allpossible_sumcount / bin_size 
        self.log10prob = math.log10(self.prob + self.pseudocount)

    def print_content(self, fout, bin_mid, bin_range_string):
        print_str = "{:.3f}\t".format(bin_mid)
        print_str += "{}\t".format(bin_range_string)
        print_str += '\t'.join('{}'.format(self.count_per_ori[ori]) for ori in self.orientation_list )
        print_str += "\t{}\t".format(self.sumcount)
        print_str += '\t'.join('{:.3f}'.format(self.log10count_per_ori[ori]) for ori in self.orientation_list )
        print_str += "\t{:.3f}\t".format(self.log10sumcount)
        print_str += '\t'.join('{:.3f}'.format(self.pcount_per_ori[ori]) for ori in self.orientation_list )
        print_str += "\t{:.3E}".format(self.allpossible_sumcount)
        print_str += "\t{:.3E}".format(self.prob)
        print_str += "\t{:.3f}\t".format(self.log10prob)
        print_str += '\t'.join('{:.3E}'.format(self.count_per_chr[chr]) for chr in self.chr_list )
        print_str += '\t'
        print_str += '\t'.join('{:.3E}'.format(self.allpossible_count_per_chr[chr]) for chr in self.chr_list )
        print_str += '\t'
        print_str += '\t'.join('{:.3E}'.format(self.prob_per_chr[chr]) for chr in self.chr_list )
        print_str += '\t'
        print_str += '\t'.join('{:.3f}'.format(self.log10prob_per_chr[chr]) for chr in self.chr_list )
        print_str += '\n'
        fout.write(print_str)

    def print_header(self, fout):
        header_str = "distance\t" \
            + '\tdistance_range(bp)\t' \
            + '\t'.join('count.{}'.format(k) for k in orientation_names) \
            + '\tsum\t' \
            + '\t'.join('log10count.{}'.format(k) for k in orientation_names) \
            + '\tlog10sum\t' \
            + '\t'.join('proportion.{}'.format(k) for k in orientation_names) \
            + '\tallpossible_sumcount' \
            + '\tprob' \
            + '\tlog10prob\t' \
            + '\t'.join('count_per_chr.{}'.format(k) for k in self.chr_list) \
            + '\t' \
            + '\t'.join('allpossible_count_per_chr.{}'.format(k) for k in self.chr_list) \
            + '\t' \
            + '\t'.join('prob_per_chr.{}'.format(k) for k in self.chr_list) \
            + '\t' \
            + '\t'.join('log10prob_per_chr.{}'.format(k) for k in self.chr_list) \
            + '\n'
        fout.write(header_str)


class DistanceBin(object):
    """class related to conversion between distance, log distance, distance bin number, bin size, etc"""

    def __init__(self, min_logdistance, max_logdistance, log_binsize):
        self.min_logdistance = min_logdistance
        self.max_logdistance = max_logdistance
        self.log_binsize = log_binsize
        self.max_bin_number = int( ( max_logdistance - log_binsize / 2 ) / log_binsize )
        self.range = range(0, self.max_bin_number+1)

    def get_bin_size(self, bin_mid):
        return(10**( bin_mid + self.log_binsize/2 ) - 10**( bin_mid - self.log_binsize/2 ))

    def get_bin_mid(self, bin_number):
        """return midpoint of a bin at log scale"""
        return(bin_number * self.log_binsize + self.log_binsize/2)

    def get_bin_number(self, distance):
        log_distance = math.log10(distance) 
        bin_number = int(log_distance / self.log_binsize)
        return(bin_number) 

    def get_bin_range_string(self, bin_mid):
        minval = int(round(10**(bin_mid - self.log_binsize/2)))
        maxval = int(round(10**(bin_mid + self.log_binsize/2)))
        return("{:,}~{:,}".format(minval, maxval))


def get_distance_and_orientation (line, cols):
    """return distance and orientation 
    given a list representing a line from the pairs input file and a ColIndices object
    """
    distance = int(line[cols.pos2]) - int(line[cols.pos1])
    
    # distance will always be > 0 for upper triangle, but in case it is not true.
    if distance > 0:
        orientation = str(line[cols.strand1]) + str(line[cols.strand2])
    else:
        orientation = str(line[cols.strand2]) + str(line[cols.strand1])
        distance = abs(distance)

    return(distance, orientation)


def cis_trans_ratio (pairs_file, outdir='report', DIST_THRES=20000, cols= cols_pairs):
    """measure cis/trans ratio for a given pairs file"""

    cts = CisTransStat()
     
    tb=pypairix.open(pairs_file)
    chrplist = tb.get_blocknames()
    for chrp in chrplist:
        it = tb.querys2D(chrp)
        chr1, chr2 = chrp.split(SEPARATOR)
        if chr1 == chr2:
            for x in it:
                distance = get_distance_and_orientation(x, cols)[0]
                if distance > DIST_THRES:
                    cts.cis += 1
                else:
                    cts.cis_short += 1
        else:
            cts.trans += sum(1 for x in it)
    cts.calculate_total()
    
    # print stats
    with open(outdir + '/' + CIS_TRANS_OUT_FILE_NAME,'w') as f:
        cts.print_stat(f)


def distance_histogram (pairs_file, chromsize_file, outdir='report', cols=cols_pairs, orientation_list = orientation_list_pairs, max_logdistance=8.4, min_logdistance=1, log_binsize=0.1):
    """create a log10-scale binned histogram table for read separation distance histogram
    The histogram is stratefied by read orientation (4 different orientations)
    The table includes raw counts, log10 counts (pseudocounts added), contact probability, log10 contact probability, and proportions for orientation (pseudocounts added)
    Bin is represented by the mid value at the log10 scale. 
    log_binsize: distance bin size in log10 scale.
    """
    gs = GenomeSize(chromsize_file)
    bins = DistanceBin(min_logdistance, max_logdistance, log_binsize)

    ss = []
    for a in bins.range:
        ss.append(SeparationStat(orientation_list,gs))

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
                    bin_number = bins.get_bin_number(distance)
                    if bin_number <= bins.max_bin_number:
                        ss[bin_number].increment(orientation, chr1)

    # calculate total
    for bin_number in bins.range:
        ss[bin_number].calculate_sumcount()

    # calculate histogram in log10 counts and proportion
    for bin_number in bins.range:
        ss[bin_number].calculate_log10count_per_ori()
        ss[bin_number].calculate_log10sumcount()
        ss[bin_number].calculate_pcount_per_ori()

    # calculate contact probability
    for bin_number in bins.range:
        bin_mid = bins.get_bin_mid(bin_number)
        bin_size = bins.get_bin_size(bin_mid)
        ss[bin_number].calculate_contact_probability_per_chr(bin_mid, bin_size)
        ss[bin_number].calculate_contact_probability(bin_mid, bin_size)

    # print histogram
    with open(outdir + '/' + PLOT_TABLE_OUT_FILE_NAME,'w') as f:
        ss[0].print_header(f)
        for bin_number in bins.range:
            bin_mid = bins.get_bin_mid(bin_number)
            if bin_mid <= bins.max_logdistance and bin_mid >= bins.min_logdistance:
                ss[bin_number].print_content(f, bin_mid, bins.get_bin_range_string(bin_mid))


if __name__ == '__main__':

    import argparse
 
    parser = argparse.ArgumentParser(description = 'QC for Pairs')
    parser.add_argument('-p','--pairs', help = "input pairs file")
    parser.add_argument('-c','--chrsize', help = "input chromsize file")
    parser.add_argument('-t','--input_type', help = "input file type (P:pairs, M:merged_nodups, OM:old_merged_nodups)")
    parser.add_argument('-O','--output_prefix', help = "prefix of output directory (output directory name will be <output_prefix>_report")
    args = parser.parse_args()
 
    if args.output_prefix:
        outdir = args.output_prefix + '_report'
    else:
        outdir = 'report'
    if not os.path.exists(outdir):
        os.mkdir(outdir)
 
    # input type selection
    if args.input_type == 'P':
        cols = cols_pairs
        orientation_list = orientation_list_pairs
    elif args.input_type == 'M':
        cols = cols_merged_nodups
        orientation_list = orientation_list_merged_nodups
    elif args.input_type == 'OM':
        cols = cols_old_merged_nodups 
        orientation_list = orientation_list_merged_nodups
    else:
        print("Unknown input type"); exit(1)
 
    # get the stats
    cis_trans_ratio (args.pairs, outdir = outdir, cols = cols )
    distance_histogram (args.pairs, args.chrsize, outdir = outdir, cols = cols, orientation_list = orientation_list)
 
 
