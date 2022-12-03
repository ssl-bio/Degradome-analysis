# Aligns RNA sequences around peaks identified from PyDegradome analysis
# against known miRNA sequences downloaded from www.mirbase.org

# Type of alignment is global, with -10, -2 open, extend penalty score
# Alignments with a score above a threshold are reported.

# Peak sequences are trimmed from both ends to match the length of
# miRNA sequences tested plus 2nt

# Execution: python Scripts/sh_py/05-Peak_miRNA_alignment.py -rd $dir -bn $ibase
# Where $dir is the input root directory and
# $ibase refers to the project base name


# Output: Aligments are saved as images (*.png)
# list of putative targets are saved as txt files
# tuples of peak indices, peak sequences ID, miRNA indices and
# alignment scores are saved as python objects (pickle)

import os  # access os directories and files
import argparse
from collections import defaultdict
import numpy as np
import itertools
from Bio import SeqIO  # parse fasta
# from Bio import pairwise2  # align
from Bio import Align  # align
from Bio.pairwise2 import format_alignment  # pretty print alignment
from Bio.Align import substitution_matrices as submat  # Subs. matrices
import multiprocessing as mp  # parallelization
import pickle  # Save objects
from PIL import ImageFont, Image, ImageDraw
import tkinter as tk
from tkinter import font as tkFont
from time import process_time


# Checks if directory exists if not it creates it
def dir_exist(ipath):
    if(not os.path.exists(ipath)):
        os.mkdir(ipath)


# Get the index of those alignments whose score is higher than a threshold
def getAlignIndx(peak, mirna, mir_index, threshold):
    """Align pairs of sequences indexed from two equal multifasta objects and
      returns the  alingment score. Sequence two is complemented"""
    mirXc = mirna[mir_index].seq.reverse_complement()  # set2

    # Set the length of the peak seq 2nt longer than miRNA
    len_diff = len(peak)-len(mirXc)
    if (len_diff % 2 == 1):
        offset = ((len_diff-1)/2)
        l_indx = int(offset)
        u_indx = int(len(peak)-offset+1)
        atX = peak[l_indx:u_indx]
    else:
        offset = ((len_diff)/2)
        l_indx = int(offset-1)
        u_indx = int((len(peak)-offset+1))
        atX = peak[l_indx:u_indx]

    # i_score = pairwise2.align.globalds(  # globalds
    #     atX, mirXc, matrix,
    #     open=-10, extend=-2, score_only=True)
    i_score = aligner.score(atX, mirXc)
    if i_score > threshold:
        return ((mir_index, i_score))


# Get the index of those alignments whose score is higher than a threshold
def getAlign(peak, mirna):
    """Align pairs of sequences indexed from two equal multifasta objects and
      returns the  alingment score. Sequence two is complemented"""
    mirXc = mirna.reverse_complement()

    # Set the length of the peak seq 2nt longer than miRNA
    len_diff = len(peak)-len(mirXc)
    if (len_diff % 2 == 1):
        offset = ((len_diff-1)/2)
        l_indx = int(offset)
        u_indx = int(len(peak)-offset+1)
        atX = peak[l_indx:u_indx]
    else:
        offset = ((len_diff)/2)
        l_indx = int(offset-1)
        u_indx = int((len(peak)-offset+1))
        atX = peak[l_indx:u_indx]

    offset = int(offset)
    # i_align = pairwise2.align.globalds(
    #     atX, mirXc, matrix,
    #     one_alignment_only=True,
    #     open=-10, extend=-2)
    i_align = aligner.align(atX, mirXc)
    return (i_align)


def drawAlign(text, outfile, font, size):
    tk.Frame().destroy()
    txt = tkFont.Font(family=font, size=size)
    ilist = text.split("\n")
    list_width = []
    for i in ilist:
        iwidth = int(txt.measure(i))
        list_width.append(iwidth)

    width = int(max(list_width)*0.88)
    height = int(width / 4)
    font = ImageFont.truetype(font.replace(" ", "-") + ".ttf",
                              size)
    img = Image.new(
        "RGBA", (width, height), (255, 255, 255, 255)  # Dimensions  # Bg color
    )
    draw = ImageDraw.Draw(img)
    draw.text((0, 0), text, (0, 0, 0), font=font)
    # draw = ImageDraw.Draw(img)
    img.save(outfile)


def irev(input_strand):
    reversed_strand = ""
    length = len(input_strand)
    for i in range(length):
        character = input_strand[length - 1 - i]
        if character == "A":
            reversed_strand = reversed_strand + "U"
        elif character == "U":
            reversed_strand = reversed_strand + "A"
        elif character == "G":
            reversed_strand = reversed_strand + "C"
        elif character == "C":
            reversed_strand = reversed_strand + "G"
        else:
            reversed_strand = reversed_strand + character
    reversed_strand = reversed_strand[::-1]
    return (reversed_strand)


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

# Image settings
img_font = "Hack Regular"
font_size = 20

# Directory definition
peakSeqDir = os.path.join(ivars["supp_data_dir"], "Peak_sequences")
mirSeqDir = os.path.join(ivars["supp_data_dir"], "miRNA_seq")
out_dir = os.path.join(ivars["output_dirR"], "03-Report/Summary")
out_dir_img_root = os.path.join(ivars["output_dirR"], "03-Report/Alignment")
dir_exist(out_dir_img_root)
out_dir_img = os.path.join(ivars["output_dirR"], "03-Report/Alignment/global")
dir_exist(out_dir_img)

# PyDegradome settings
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

# Aligner settings
aligner = Align.PairwiseAligner()
aligner.substitution_matrix = nuc44
aligner.extend_gap_score = -2
aligner.open_gap_score = -10

# Parse fasta files
mirSeq = SeqIO.parse(os.path.join(mirSeqDir, "input/miRNA_sequences.fa"),
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
                           "output/Alignment_indices_MF-" + str(iMF) +
                           "_iConf-" + iConf2 + ".obj")
    out_file = os.path.join(mirSeqDir,
                            "output/miRNA_targets_MF-" + str(iMF) +
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
            peakHeader = []
            for iseq in peakSequences:
                peakRegions.append(iseq)
                peakHeader.append(iseq.description)

            if os.path.exists(out_obj):
                fileobj = open(out_obj, "rb")
                l_index_filtered = pickle.load(fileobj)
                fileobj.close()
            else:
                peakIndx = []
                for ihead in peakHeader:
                    tmp = ihead.split()
                    itx = tmp[0].replace(".", "_")
                    icomp = tmp[1]
                    ipeak = tmp[2].split(":")[1]
                    peakIndx.append(icomp + itx + ipeak)

                # Align
                t1_start = process_time()
                irange = range(len(mirSequences))
                l_index = []
                for i in range(len(peakRegions)):
                    ipeak = peakRegions[i].seq
                    pool = mp.Pool(mp.cpu_count()-1)
                    aln_ind = [pool.apply(getAlignIndx,
                                          args=(ipeak,
                                                mirSequences,
                                                j, 33)) for j in irange]
                    pool.close()
                    aln_ind_filtered = [i for i in aln_ind if i is not None]
                    results = [aln_ind_filtered, i, peakIndx[i]]
                    l_index.append(results)
                t1_stop = process_time()
                print("Mapping time:", t1_stop, t1_start)

                # Check which sequences may be a miRNA target
                l_index_filtered = []
                for item in l_index:
                    if len(item[0]) > 0:
                        isumm = {"peak_indx": item[1],
                                 "peak_indxFull": item[2],
                                 "miRNA_indx": item[0][0][0],
                                 "algn_score": item[0][0][1]}
                        l_index_filtered.append(isumm)

                # Save indices
                fileObj = open(out_obj, 'wb')
                pickle.dump(l_index_filtered, fileObj)
                fileObj.close()

            if len(l_index_filtered) > 0:
                # Get full alignment
                l_algn = []
                for idict in l_index_filtered:
                    ialign = getAlign(peakRegions[idict["peak_indx"]].seq,
                                      mirSequences[idict["miRNA_indx"]].seq)
                    l_algn.append((ialign,
                                   idict["miRNA_indx"],
                                   idict["peak_indx"],
                                   idict["peak_indxFull"]))

                # Output summary
                print('# List of pydegradome identified targets whose peak region align with a miRNA', file=open(out_file, 'w'))

                for item in l_algn:
                    ialign = item[0]
                    iscore = str(ialign[0].score)
                    mir_indx = item[1]
                    peak_indx = item[2]
                    peak_indxFull = item[3]
                    icomp = peakRegions[peak_indx].description.split(" ")[1]
                    itx = peakRegions[peak_indx].name
                    imir = mirSequences[mir_indx].name
                    # ireport = format_alignment(*ialign[0])
                    # ireport = str(ialign[0]).replace("T","U") + str(ialign[0].score)
                    ialign0 = str(ialign[0]).replace("T", "U").split("\n")
                    imir_seq = irev(ialign0[2])
                    if len(itx) <= len(imir):
                        isep = " " * (len(imir) + 4)
                        ispc = " " * (len(imir) - len(itx))
                        pre1 = itx + ispc + " 5' "
                        pre2 = imir + " 3' "
                    else:
                        isep = " " * (len(itx) + 4)
                        ispc = " " * (len(itx) - len(imir))
                        pre1 = itx + " 5' "
                        pre2 = imir + ispc + " 3' "

                    ireport = pre1+ialign0[0]+" 3' "+"\n" + \
                        isep+ialign0[1] + "\n"+pre2 + imir_seq + " 5' " + \
                        "\n\n"+"Score: " + iscore

                    # Save alignment as image
                    out_img = "_".join(["Aln_global",
                                        itx.replace(".", "_"),
                                        imir, iConf2, str(iMF)]) + ".png"
                    comp_dir = os.path.join(out_dir_img, icomp)
                    dir_exist(comp_dir)
                    img_path = os.path.join(comp_dir, out_img)
                    drawAlign(ireport, img_path, img_font, font_size)

                    # Write report
                    print(f'Comparison: {icomp}\nTranscript: {itx}\nmiRNA: {imir}\nScore: {iscore}\nIndex: {peak_indxFull}',
                          file=open(out_file, 'a'))
