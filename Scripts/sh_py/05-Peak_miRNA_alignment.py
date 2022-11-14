import os  # access os directories and files
import argparse
from collections import defaultdict
import numpy as np
import itertools
from Bio import SeqIO  # parse fasta
from Bio import pairwise2  # align
# from Bio.pairwise2 import format_alignment  # pretty print alignment
from Bio.Align import substitution_matrices as submat  # Subs. matrices
import multiprocessing as mp  # parallelization
import pickle  # Save objects
from time import process_time


# Get the index of those alignments whose score is higher than a threshold
def getAlignIndx(seq1, seq2, index, threshold, matrix):
    """Align pairs of sequences indexed from two equal multifasta objects and
      returns the  alingment score. Sequence two is complemented"""
    atX = seq1
    mirXc = seq2[index].seq.reverse_complement()  # set2
    len_diff = len(atX)-len(mirXc)
    if (len_diff % 2 == 1):
        offset = ((len_diff-1)/2) - 1
    else:
        offset = ((len_diff)/2) - 1

    offset = int(offset)
    i_score = pairwise2.align.globalds(atX[offset:], mirXc, matrix, open=-10, extend=-2, score_only=True)
    if i_score > threshold:
        return ((index, i_score))


# Parse input argument
parser = argparse.ArgumentParser()

parser.add_argument(
    "-rd",
    "--rootdir",
    # default=Droot_dir,
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

peakSeqDir = os.path.join(ivars["supp_data_dir"], "Peak_sequences")
mirSeqDir = os.path.join(ivars["supp_data_dir"], "miRNA_seq/input")
out_dir = os.path.join(ivars["output_dirR"], "03-Report/Summary")

pydeg_settings = np.array(ivars["pydeg_script_settings"].
                          replace('"', '').
                          strip(')(').split(' '))
isel = list(
    itertools.chain.from_iterable(
        itertools.repeat([True, False],
                         int(len(pydeg_settings)/2))))
isel_not = [not i for i in isel]
iMF_list = [int(x) for x in pydeg_settings[np.array(isel_not)]]
iConf_list = [float(x) for x in pydeg_settings[np.array(isel)]]


nuc44 = submat.load("NUC.4.4")  # substitution matrix

# Parse fasta files
mirSeq = SeqIO.parse(os.path.join(mirSeqDir, "miRNA_sequences.fa"),
                     "fasta")

mirSequences = []
for iseq in mirSeq:
    mirSequences.append(iseq)

for i in range(len(iMF_list)):
    iMF = iMF_list[i]
    iConf = iConf_list[i]
    iConf2 = str(iConf).replace(".", "_")

    # Output file names
    out_obj = os.path.join(mirSeqDir,
                           "Alignment_indices_MF-" + str(iMF) +
                           "_iConf-" + iConf2 + ".obj")
    out_file = os.path.join(out_dir,
                            "miRNA_targets_MF-" + str(iMF) +
                            "_iConf-" + iConf2 + ".txt")

    if os.path.exists(out_obj) and os.path.exists(out_file):
        print(f"""Files:
    - {os.path.basename(out_obj)}
    - {os.path.basename(out_file)}     
Already exists""")
    else:
        ifasta = os.path.join(peakSeqDir, "PeakRegioncDNA_category_1_" +
                              iConf2 + "_4_" + str(iMF) + ".fa")
        if os.path.exists(ifasta):
            peakSequences = SeqIO.parse(ifasta, "fasta")

            # Create lists of sequences
            peakRegions = []
            for iseq in peakSequences:
                peakRegions.append(iseq)

            # Align
            t1_start = process_time()
            pool = mp.Pool(mp.cpu_count()-1)
            irange = range(len(mirSequences))
            l_index = []
            for i in range(len(peakRegions)):
                ipeak = peakRegions[i].seq
                aln_ind = [pool.apply(getAlignIndx,
                                      args=(ipeak,
                                            mirSequences,
                                            j, 59, nuc44))
                           for j in irange]
                aln_ind_filtered = [i for i in aln_ind if i is not None]
                results = [aln_ind_filtered, i]
                l_index.append(results)
            pool.close()
            t1_stop = process_time()
            print("Mapping time:", t1_stop, t1_start)

            # Check which sequences may be a miRNA target
            l_index_filtered = []
            for item in l_index:
                if len(item[0]) > 0:
                    l_index_filtered.append(item)

            # Save indices
            fileObj = open(out_obj, 'wb')
            pickle.dump(l_index_filtered, fileObj)
            fileObj.close()

            # Output summary
            print('# List of pydegradome identified targets whose peak region align with a miRNA', file=open(out_file, 'w'))

            for i in l_index_filtered:
                list_miRNAs = i[0]
                peak_indx = i[1]
                for j in list_miRNAs:
                    mir_indx = j[0]
                    print(f'Transcript: {peakRegions[peak_indx]. name}\nComparison: {peakRegions[peak_indx].description.split(" ")[1]}\nmiRNA: {mirSequences[mir_indx].name}', file=open(out_file, 'a'))            
