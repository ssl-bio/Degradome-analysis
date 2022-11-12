# PyDegradome, a program to find cuts site in degradome-seq/PARE experiments
# with control and test datasets
#
# Authors   : Marta Gaglia (Tufts Univ.), Chris H. Rycroft (Harvard/LBL)
# Email     : marta.gaglia@tufts.edu

from sys import *
from math import *
from bisect import bisect
from scipy import interpolate
import numpy as np
from thr_class import *
import re
import argparse

# Defining external arguments
parser = argparse.ArgumentParser(
    prog='PyDegradome.py',
    description='''
Identifies coordinates in a test sample where the number of reads
are significatively higher than those of the corresponding control sample.
Library requirements: NumPy, SciPy
Authors: Marta Gaglia (Tufts Univ.), Chris H. Rycroft (Harvard/LBL)
''')

parser.add_argument(
    "-gtf",
    "--gtf_genome",
    dest="gtf",
    required=True,
    help="Genome annotation file",
)

parser.add_argument(
    "-ctrl",
    "--sam_ctrl",
    dest="sam_ctrl",
    required=True,
    help="Mapped file (sam) from control sample"
)

parser.add_argument(
    "-test",
    "--sam_test",
    dest="sam_test",
    required=True,
    help="Mapped file (sam) from test sample"
)

parser.add_argument(
    "-iconf",
    "--conf_level",
    dest="iconf",
    required=True,
    help="Confidence level (e.g. 0.95)"
)

parser.add_argument(
    "-w",
    "--window",
    default=4,
    dest="window",
    help="window size (Integer, default 4) for merging reads"
)

parser.add_argument(
    "-imf",
    "--mult_factor",
    dest="imf",
    required=True,
    help="Multiplicative factor (integer) to calculate significance threshold"
)

parser.add_argument(
    "-t",
    "--tot_sites",
    dest="tot_sites",
    default=32348300,
    help="Estimate for the total number of sites that could possibly be cut\
(1% of human genome)"
)

parser.add_argument(
    "-o",
    "--out_file",
    dest="out_file",
    required=True,
    help="Output filename"
)

args = parser.parse_args()


# Check if help requested
# if len(argv) > 1 and (argv[1] == "-h" or argv[1] == "--help"):
#     help_message()

# # Check for enough command-line arguments
# if len(argv) < 8:
#     error_message()

# # Parse any options
# while len(argv) > 8:
#     cmd=argv.pop(1)
#     if cmd == "-h" or cmd == "--help":
#         help_message()
#     elif cmd == "-t":
#         if argv == 8:
#             error_message()
#             tot_sites=float(argv.pop(1))
#     else:
#         error_message()

# Default value for the total number of sites that could possibly be cut, based
# on 1% of the human genome
tot_sites = float(args.tot_sites)

# Define class to store annotation, peaks indexed by chromosome
class chromosomes:
    def __init__(self, ID):
        self.ID = ID
        self.exonstart = [[], []]
        self.exonstop = [[], []]

        # key = exon indentifier; value = list of reads (5') that map to that exon;
        # position in list: 0 = + strand, 1 = - strand
        self.exoncov = [{}, {}]

        # key = exon indentifier, value = list of lists, one for each sample,
        # each as long as the exon length;
        # position in list: 0 = + strand, 1 = - strand
        self.exon_l = [{}, {}]

        # key = exon indentifier, value = list of hit windows [start,stop];
        # position in list: 0 = + strand, 1 = - strand
        self.hit = [{}, {}]

        # key = exon indentifier, value = list of hit windows [start,stop];
        # position in list: 0 = + strand, 1 = - strand
        self.hitrefined = [{}, {}]
        self.order = [[], []]

    # Collapses alternative exons into a big one
    def add_exons(self, strand, start, stop):
        if strand == "+":
            ind = 0
        elif strand == "-":
            ind = 1
        if len(self.exonstart[ind]) > 0:
            if start < self.exonstart[ind][-1]:
                print start, self.exonstart[ind][-1]
                exit()
            if start < self.exonstop[ind][-1]:
                if stop > self.exonstop[ind][-1]:
                    self.exonstop[ind][-1] = stop
            else:
                self.exonstart[ind].append(start)
                self.exonstop[ind].append(stop)
        else:
            self.exonstart[ind].append(start)
            self.exonstop[ind].append(stop)

    def add_hit(self, s, ex_id, p1, p2):
        if ex_id not in self.hit[s]:
            self.hit[s][ex_id] = []
            self.hit[s][ex_id].append([p1, p2])

           
# Function to find which exon the aligned read comes from
def find_exon(ex_start_list, ex_stop_list, n):
    poss_ex = []
    # support for maintaining a list in sorted order
    # without having to sort the list after each insertion
    a = bisect(ex_start_list, n) - 1
    if a >= 0 and n <= ex_stop_list[a]:
        poss_ex.append(ex_start_list[a])
        poss_ex.append(ex_stop_list[a])
    return poss_ex


# Adds to read count per exon and per position
def add_cov(l, d1, d2, s, p):
    k = str(l[0])+"-"+str(l[1])
    d1[k][s] += 1
    p2 = p-l[0]
    d2[k][s][p2] += 1

def parsecigar(cigarstring, pos_ref):
        '''
        Modifies a sequence according to its CIGAR string, modified for speed according to frequency of CIGAR tags

        :param cigarstring: cigar string, v3 or v4 are allowed
        :param pos_ref: position in the reference of the leftmost aligned nucleotide (0-based)
        :return: position
        '''

        matches = re.findall(r"(\d+)([MIDNSHPX=]{1})", cigarstring)
        cigar = [{"type": m[1], "length": int(m[0])} for m in matches]
        start = 0
        start_ref = int(pos_ref)
        for c in range(0, len(cigar)):
                l = cigar[c]["length"]
                if cigar[c]["type"] in ["=", "X", "M", "N"]:
                        # seqout += seq[start:start + l]
                        # start += l
                        start_ref += l
                elif cigar[c]["type"] in ["N", "P"]:
                        start_ref += l
                elif cigar[c]["type"] == "S":
                        start += l
        return start_ref
    
# Set up the gd dictionary based on the gtf file and record position of exons
# in each chromosome
stderr.write("processing annotation:\n")

gtf = open(args.gtf, "r")
gd = {}

for l in gtf:
    l = l.strip()
    if not l.startswith("#"):
        sl = l.split("\t")
        if sl[2] == "exon":
            chrm = sl[0]
            if chrm not in gd:
                gd[chrm] = chromosomes(chrm)
            strand = sl[6]
            start = int(sl[3])
            stop = int(sl[4])
            gd[chrm].add_exons(strand, start, stop)

# Create a list for each exon to record the total number of reads per exon and
# a list to record the counts at each position
for chrm in gd:
    stderr.write(chrm+"\n")
    for ind in [0, 1]:
        temp = gd[chrm].exonstart[ind]
        temp2 = gd[chrm].exonstop[ind]
        for x in range(len(temp)):
            k = str(temp[x])+"-"+str(temp2[x])
            gd[chrm].exoncov[ind][k] = [0, 0]
            ex_len = temp2[x]-temp[x]+1
            gd[chrm].exon_l[ind][k] = []
            gd[chrm].exon_l[ind][k].append([0]*ex_len)
            gd[chrm].exon_l[ind][k].append([0]*ex_len)

# Go through the read file, identify the exon the read matches to, generate a
# list for that exon if needed and count where the 5' most end of the read
# falls
stderr.write("\nprocessing control and test read files:\n")

index_list = [args.sam_ctrl, args.sam_test]  # sets input file list

tcr = [0, 0]  # records total # of aligned reads in input
cr = [0, 0]  # records total # of reads used in analysis

p = re.compile("^\d*M$")

for s in range(len(index_list)):
    filename = index_list[s]
    stderr.write(filename+"\n")
    f = open(filename, "r")
    for l in f:
        if l[0] != "@":
            l = l.strip()
            sl = l.split("\t")
            pos = "not set"
            tcr[s] += 1
            icigar = sl[5]
            pos_ref = sl[3]
            for fl in range(11, len(sl)):
                if sl[fl][0:2] == "XS":
                    XS = sl[fl]

            # For simplicity, reads with insertions or deletions relative to
            # the reference are excluded
            if icigar in ["I", "D", "H"]:
                continue

            if int(sl[1]) & 16:
               
                # Removes reads that mapped to strand complementary to the
                # reference
                if XS == "XS:A:-":
                    x = 1
                    if (p.search(icigar)):
                        pos = int(pos_ref) + len(sl[9]) - 1
                    else:
                        # int1 = sl[5].split("M")
                        # int2 = int1[1].split("N")
                        # pos = int(pos_ref) + int(int1[0]) + int(int2[0]) + int(int2[1]) - 1
                        pos = parsecigar(icigar, pos_ref)
            else:
                # Removes reads that mapped to strand complementary to the
                # reference
                if XS == "XS:A:+":
                    x = 0
                    pos = int(pos_ref)
            chrm = sl[2]
            if pos != "not set":
                poss_ex = find_exon(gd[chrm].exonstart[x],
                                    gd[chrm].exonstop[x], pos)
                if len(poss_ex) != 0:
                    add_cov(poss_ex, gd[chrm].exoncov[x],
                            gd[chrm].exon_l[x], s, pos)
                    cr[s] += 1

tcr_s = "total reads control: "+str(tcr[0])+"\ntotal reads test: "+str(tcr[1])+"\n"
cr_s = "reads used control: "+str(cr[0])+"\nreads used test: "+str(cr[1])+"\n"
stderr.write(tcr_s)
stderr.write(cr_s)
stderr.write("\nparameters of prior distributions:\n")

# Generate a table of number of positions with a certain read count to it
#
# First generate a list where pos = read count, value = freq of that read count
# 0 = control, 1 = test
chd = [{},{}]

# Scanning window size to integrate over
w = int(args.window)

# Record the number of windows with specific read count as a dictionary with
# key = read count, value = # of instances (one for test, and one for control
# samples)
for chrm in gd:
    for x in [0, 1]:
        for ex in gd[chrm].exon_l[x]:
            for z in [0, 1]:
                temp_hl = gd[chrm].exon_l[x][ex][z]
                for nt in range(0, len(temp_hl)-w+1):
                    w_hl = sum(temp_hl[nt:nt+w])
                    ch = l_scale(w_hl)
                    if ch != 0:
                        if ch not in chd[z]:
                            chd[z][ch] = 0
                        chd[z][ch] += 1

# Convert the dictionary to a list
chd_l = [[], []]
for x in [0, 1]:
    max_val = 1+max(chd[x].keys())
    chd_l[x] = [0]*max_val
    for k in chd[x]:
        chd_l[x][k] = chd[x][k]

# Print a table of the list
fn_l = ["ctl", "test"]
for x in [0, 1]:
    fn = args.out_file + "_" + fn_l[x] + "_heightfreq.txt"
    fh = open(fn, "w")
    for y in range(len(chd_l[x])):
        nl = [str(y), str(chd_l[x][y])]
        nlj = "\t".join(nl)
        fh.write(nlj+"\n")
    fh.close()

# Generate the threshold table using the lists above with the given input
# confidence level (stored in argv[4])
th = thresh_tab()
th.fit_prior(chd_l, tot_sites)
th.make_table(float(args.iconf))
th.save_table(args.out_file+"_test.tab")

# Go through counts in exons, identify count in controls, test whether count in
# test exceed the threshold by looking up the threshold in the generated table,
# and record windows that passed
count = 0
co = float(args.imf)

stderr.write("finding peaks:\n")

for chrm in gd:
    for x in [0,1]:
        for ex in gd[chrm].exon_l[x]:
            count += 1
            if count % 100000 == 0:
                stderr.write("processed "+str(count)+" exons\n")
            flag = 0
            for ds in [0,1]:
                if gd[chrm].exoncov[x][ex][ds] > 10:
                    flag += 1
            if flag == 2:
                ex_cov_ratio = float(gd[chrm].exoncov[x][ex][1])/float(gd[chrm].
                                                                       exoncov[x][ex][0])
                factor = co*ex_cov_ratio

                # If the factor is larger 100, it will be artificially set to
                # 100 and a warning message will be printed
                if factor > 100:
                    err_mess = "Ratio of " + str(factor) + " found, chrm "+ chrm+", exon "+ex
                    stderr.write(err_mess+"\n")
                    factor = 100
                temp_t = gd[chrm].exon_l[x][ex][1]
                temp_c = gd[chrm].exon_l[x][ex][0]
                for nt in range(0, len(temp_c) - w + 1):
                    w_t = sum(temp_t[nt:nt+w])
                    w_c = sum(temp_c[nt:nt+w])
                    if w_t != 0 or w_c != 0:
                        stat = th.thr(w_c, factor)
                        if w_t >= stat and max(temp_t[nt:nt+w]) > 1:
                            gd[chrm].add_hit(x, ex, nt, nt+w)

# Trim windows so that 0 or 1 count positions are eliminated. For example,
# [0,1,2,4,5] becomes [2,4,5].
for c in gd:
    for x in [0, 1]:
        for ex in gd[c].hit[x]:
            temp = gd[c].hit[x][ex]
            cov = gd[c].exon_l[x][ex][1]
            if len(temp) != 0:
                for i in range(0, len(temp)):
                    a = cov[temp[i][0]:temp[i][1]]
                    c1 = 0
                    c2 = 0
                    for j in range(0, len(a)):
                        if a[j] < 2:
                            c1 += 1
                        else:
                            break
                    for j in range(len(a) - 1, -1, -1):
                        if a[j] < 2:
                            c2 += 1
                        else:
                            break
                    gd[c].hit[x][ex][i] = [temp[i][0]+c1, temp[i][1]-c2]

# Merge windows into bigger peaks
for c in gd:
    for x in [0, 1]:
        for ex in gd[c].hit[x]:
            temp_hit = gd[c].hit[x][ex]
            if len(temp_hit) != 0:
                prev_h = temp_hit[0][0]
                temp_start = [prev_h]
                temp_end = []
                for i in range(1, len(temp_hit)):
                    next_h = temp_hit[i][0]
                    if next_h > prev_h+10:
                        temp_start.append(next_h)
                        temp_end.append(temp_hit[i-1][1])
                    prev_h = next_h
                temp_end.append(temp_hit[-1][1])
                for j in range(0, len(temp_start)):
                    if ex not in gd[c].hitrefined[x]:
                        gd[c].hitrefined[x][ex] = []
                    gd[c].hitrefined[x][ex].append([temp_start[j],
                                                    temp_end[j]])

# Construct main output file
#
# Format is 3 lines per peak:
#
# The first line includes chr, strand, position of peak start, position of peak
# end, the position with the highest read count, the ratio between all reads in
# the exon in test vs. control, the factor used for thresholding (user-defined
# ratio * ratio of exon reads), the read count at the position listed as
# "position with the highest read count" in test and control sample, the total
# read count in the peak in the test and control samples.
#
# The second and third lines print out the list of counts for positions
# starting 10 nt 5' of peak start and ending 10 nt 3' of peak end (fewer if the
# peak is less than 10 nt away from the end of an annotated exon) for each of
# control and test dataset.

# Open file and print header line
outf = open(args.out_file, "w")
outf.write("#chr\tstrand\tpeak_start\tpeak_stop\tmax_peak_test\tratio\tfactor\tmax_count_test\tmax_count_ctl\ttot_count_test\ttot_count_ctl\n")

chrm_list = gd.keys()
chrm_list.sort()

strand_l = ["+", "-"]

for c in chrm_list:
    for x in [0, 1]:
        order = gd[c].hitrefined[x].keys()
        order.sort()
        for ex in order:
            temp_fin = gd[c].hitrefined[x][ex]
            cov_t = gd[c].exon_l[x][ex][1]
            cov_c = gd[c].exon_l[x][ex][0]
            info = ex.split("-")
            ex_start = int(info[0])
            strand = strand_l[x]
            if len(temp_fin) != 0:
                for h in temp_fin:
                    test = cov_t[h[0]:h[1]]
                    control = cov_c[h[0]:h[1]]
                    maxtest = max(test)
                    maxcontrol = max(control)
                    peak_rel = test.index(maxtest)
                    start = ex_start + h[0]
                    end = ex_start + h[1]
                    peak = start + peak_rel
                    tot_t = sum(test)
                    tot_c = sum(control)
                    ratio = float(gd[c].exoncov[x][ex][1])/float(gd[c].exoncov[x][ex][0])
                    factor = ratio*co
                    new_line = [c, strand, str(start), str(end),
                                str(peak), str(ratio), str(factor),
                                str(maxtest), str(maxcontrol),
                                str(tot_t), str(tot_c)]
                    nlj2 = "\t".join(new_line)
                    outf.write(nlj2+"\n")
                    if h[0]-10 > 0:
                        outf.write(str(cov_c[h[0]-10:h[1]+10])+"\n")
                        outf.write(str(cov_t[h[0]-10:h[1]+10])+"\n")
                    else:
                        outf.write(str(cov_c[0:h[1]+10])+"\n")
                        outf.write(str(cov_t[0:h[1]+10])+"\n")

outf.close()
