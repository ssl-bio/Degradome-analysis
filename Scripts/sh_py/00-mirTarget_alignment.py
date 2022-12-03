# Test script to align known miRNA sequences and their targets
# pairs of sequences were downloaded from TarDB (www.biosequencing.cn/TarDB/)
# Liu et al. 2021 BMC Genomics 22: 348
# Two fasta files were obtained from the downloaded files

# Type of alignment is global, with -10, -2 open, extend penalty score

# Execution:
# python Scripts/sh_py/00-mirTarget_alignment.py -rd $(pwd) -bn $ibase
# Where $(pwd) is the input root directory and
# $ibase refers to the project base name

# Output: summary of the alignment scores is saved as txt file
# The lower value was used as threshold for the alignment of peak
# sequences and miRNAs (05-Peak_miRNA_alignment.py)


import argparse  # Parse arguments
from collections import defaultdict
from Bio import SeqIO  # parse fasta
# from Bio import pairwise2  # align
from Bio import Align  # align
import os  # access os directories and files
from Bio.Align import substitution_matrices as submat
import multiprocessing as mp  # parallelization
import pandas as pd  # Get summary stats
from time import process_time


# Get the alignment score
def getAlignScore(seq1, seq2):
    """Align pairs of sequences indexed from two equal multifasta objects and
      returns the  alingment score. Sequence two is complemented"""
    atX = seq1
    mirXc = seq2.complement()
    # i_align = pairwise2.align.globalds(atX, mirXc, matrix, open=-10,
    #                                    extend=-2, score_only=True)
    i_score = aligner.score(atX, mirXc)
    return (i_score)


# Parse input argument
parser = argparse.ArgumentParser()

parser.add_argument(
    "-rd",
    "--rootdir",
    dest="root_dir",
    nargs="?",
    help="Root directory",
)

parser.add_argument(
    "-bn",
    "--basename",
    dest="base_name",
    nargs="?",
    help="Name of the project",
)

args = parser.parse_args()

os.chdir(args.root_dir)

# Import variables
ifile = open("Env_variables/Degradome_"+args.base_name+".txt", "r")
ivars = defaultdict(str)
for line in ifile:
    ivar = line.strip().split('=')
    ivars[ivar[0].strip()] = ivar[1].strip()

# Define directories
SeqDir = os.path.join(ivars["basedir"], "Genetic_data/Fasta")
out_dir = os.path.join(ivars["basedir"], "Genetic_data/Others")

# Define output files
out_file = os.path.join(
    out_dir, "Alignment_score_summary_"+ivars["sp"]+"_miRNA.txt")

# Parse fasta files
mirSeq = SeqIO.parse(os.path.join(SeqDir, "ath_mir.fa"), "fasta")
athSeq = SeqIO.parse(os.path.join(SeqDir, "ath_tar.fa"), "fasta")

nuc44 = submat.load("NUC.4.4")  # substitution matrix

# Aligner settings
aligner = Align.PairwiseAligner()
aligner.substitution_matrix = nuc44
aligner.extend_gap_score = -2
aligner.open_gap_score = -10

mirSequences = []
for iseq in mirSeq:
    mirSequences.append(iseq)

AtSequences = []
for iseq in athSeq:
    AtSequences.append(iseq)

t1_start = process_time()
pool = mp.Pool(mp.cpu_count()-1)
irange = range(len(mirSequences))
l_score = [pool.apply(getAlignScore,
                      args=(AtSequences[i].seq,
                            mirSequences[i].seq)
                      ) for i in irange]
pool.close()
t1_stop = process_time()
print("Elapsed time:", t1_stop, t1_start)


# Get summary statistics of alignments
s = pd.Series(l_score)

isumm = s.describe()
fileObj = open(out_file, 'w')  # wb if python 2
fileObj.writelines(str(isumm))
fileObj.close()
