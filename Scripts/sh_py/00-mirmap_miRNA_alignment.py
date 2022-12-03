# Test script to align known miRNA sequences and their targets
# pairs of sequences were downloaded from TarDB (www.biosequencing.cn/TarDB/)
# Liu et al. 2021 BMC Genomics 22: 348
# Two fasta files were obtained from the downloaded files

# The script uses a library called 'mirmap' (https://mirmap.ezlab.org/)
# Vejnar and Zdobnov, 2012 Nucleic Acids Research 40: 11673-11683

# The library requires a python 2 enviroment
# Other dependencies include: 'viennarna' 'unittest2' 'dendropy'

# Execution:
# python Scripts/sh_py/00-mirmap_miRNA_alignment.py -rd $(pwd) -bn $ibase
# Where $(pwd) is the input root directory and
# $ibase refers to the project base name

# Output: Aligments are saved as python object (pickle)
# Summary of 'DeltaG open' values are saved as txt file
# These are intended to be used as reference for the alignment between
# peak sequences and miRNAs (05-Peak_miRNA_mirmap.py)

import os
import sys
import argparse
from Bio import SeqIO  # parse fasta
from collections import defaultdict
import multiprocessing as mp  # parallelization
import pickle  # Save objects
import pandas as pd  # Get summary stats


# def f_mirmap(itarget, imirna):
#     # Append extra bases
#     seq_target = "GC" + itarget
#     seq_mirna = imirna[::-1]
#     # Create the mm object
#     mim = mirmap.mm(seq_target, seq_mirna)

#     # Search for seeds with defined parameters
#     mim.find_potential_targets_with_seed(
#         allowed_lengths=[5, 6],
#         allowed_gu_wobbles={5: 1, 6: 1},
#         allowed_mismatches={5: 1, 6: 1},
#         take_best=True,
#     )
#     mim.eval_tgs_au(with_correction=False)  # AU content
#     # mim.eval_tgs_position(with_correction=False)  # UTR possition
#     # mim.eval_tgs_pairing3p(with_correction=False)  # 3'pairing
#     mim.eval_dg_open()  # Delta G open
#     mim.eval_prob_binomial()  # P binomial
#     return mim


def f_mirmapindx(itarget, imirna, indx):
    # Append extra bases
    if len(itarget) <= len(imirna)-1:
        seq_target = "GC" + itarget
    else:
        seq_target = itarget

    seq_mirna = imirna[::-1]
    # Create the mm object
    mim = mirmap.mm(seq_target, seq_mirna)

    # Search for seeds with defined parameters
    mim.find_potential_targets_with_seed(
        allowed_lengths=[11, 12],
        allowed_gu_wobbles={11: 1, 12: 1},
        allowed_mismatches={11: 1, 12: 1},
        take_best=True,
    )
    if len(mim.report()) > 0:
        return (indx)


# Get the mirmap object
def f_mirmap(itarget, imirna):
    # Append extra bases
    if len(itarget) <= len(imirna)-1:
        seq_target = "GC" + itarget
    else:
        seq_target = itarget

    seq_mirna = imirna[::-1]
    # Create the mm object
    mim = mirmap.mm(seq_target, seq_mirna)

    # Search for seeds with defined parameters
    mim.find_potential_targets_with_seed(
        allowed_lengths=[5, 6],
        allowed_gu_wobbles={5: 1, 6: 1},
        allowed_mismatches={5: 1, 6: 1},
        take_best=True)
    if len(mim.report()) > 0:
        mim.libs = mirmap.library_link.LibraryLink(
            os.path.join(ivars['mirmap_script_dir'],
                         'libs/lib-archlinux-x86_64'))
        mim.exe_path = os.path.join(ivars['mirmap_script_dir'],
                                    'libs/exe-archlinux-x86_64')
        mim.eval_dg_open()  # Delta G open
        mim.eval_dg_duplex()
        mim.eval_prob_binomial()  # P binomial
        return (mim)

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

# Add mirmap path to sys
sys.path.append(ivars["mirmap_script_dir"])
import mirmap
import mirmap.library_link

# Define directories
SeqDir = os.path.join(ivars["basedir"], "Genetic_data/Fasta")
out_dir = os.path.join(ivars["basedir"], "Genetic_data/Others")

# Define output files
out_obj = os.path.join(out_dir, "mirmap_indices_" + ivars["sp"] + "_miRNA.obj")
out_file = os.path.join(out_dir,
                        "DeltaG_duplex_summary_" + ivars["sp"] + "_miRNA.txt")

# Parse fasta files
mirSeq = SeqIO.parse(os.path.join(SeqDir, "ath_mir.fa"), "fasta")
mirSequences = []
for iseq in mirSeq:
    mirSequences.append(iseq)

atSeq = SeqIO.parse(os.path.join(SeqDir, "ath_tar.fa"), "fasta")
atSequences = []
for iseq in atSeq:
    atSequences.append(iseq)

pool = mp.Pool(mp.cpu_count() - 1)
mm_indices = [
    pool.apply(
        f_mirmapindx,
        args=(
            str(atSequences[indx].seq).replace("T", "U"),
            str(mirSequences[indx].seq).replace("T", "U"),
            indx
        ),
    )
    for indx in range(len(atSequences))
]
pool.close()
mm_indices_filtered = [
    i for i in mm_indices if i is not None]

mm_list = []
for i in mm_indices_filtered:
    itarget = str(atSequences[i].seq).replace("T", "U")
    imirna = str(mirSequences[i].seq).replace("T", "U")
    imm = f_mirmap(itarget, imirna)
    mm_list.append(imm)

mm_list_filtered = [
    i for i in mm_list if i is not None]
# ----------------------------------
dg_duplex_list = []
for i in range(len(mm_list_filtered)):
    try:
        dg_duplex_list.append(mm_list_filtered[i].dg_duplex)
    except ValueError:
        print("No alignment was found for item {}".format(i))

# Get summary statistics of alignments
s = pd.Series(dg_duplex_list)
isumm = s.describe()
fileObj = open(out_file, "wb")
fileObj.writelines(str(isumm))
fileObj.close()


# Save indices
fileObj = open(out_obj, "wb")
pickle.dump(mm_indices_filtered, fileObj)
fileObj.close()
