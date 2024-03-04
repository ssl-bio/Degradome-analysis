# 05-Peak_miRNA_mirmap.py
'''
Aligns RNA sequences around peaks identified from PyDegradome analysis
against known miRNA sequences downloaded from www.mirbase.org

The script uses a library called 'mirmap' (https://mirmap.ezlab.org/)
Vejnar and Zdobnov, 2012 Nucleic Acids Research 40: 11673-11683

The library requires a python 2 enviroment
Other dependencies include: 'viennarna' 'unittest2' 'dendropy'

Execution: python Scripts/sh_py/05-Peak_miRNA_mirmap.py -rd $dir -bn $ibase
Where $dir is the input root directory and
$ibase refers to the project base name


Output: Aligments are saved as images (*.png)
list of putative targets are saved as txt files
tuples of alignments, miRNA indices, peak indices and ID are saved as
python objects (pickle)
'''

import os  # access os directories and files
import sys
import argparse
from collections import defaultdict
import numpy as np
import itertools
from Bio import SeqIO  # parse fasta
import multiprocessing as mp  # parallelization
import pickle  # Save objects
from PIL import ImageFont, Image, ImageDraw
import Tkinter as tk
import tkFont


# Checks if directory exists if not it creates it
def dir_exist(ipath):
    if not os.path.exists(ipath):
        os.makedirs(ipath)


# Get the index of those alignments whose score is higher than a threshold
def f_mirmapindx(itarget, imirna, mirindx, peakindx):
    # Append extra bases
    if len(itarget) <= len(imirna) - 1:
        seq_target = "GC" + itarget
    else:
        seq_target = itarget

    seq_mirna = imirna  # [::-1]
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
        return (mirindx, peakindx)


# Get the mirmap object
def f_mirmap(itarget, imirna, threshold):
    # Append extra bases
    if len(itarget) <= len(imirna) - 1:
        seq_target = "GC" + itarget
    else:
        seq_target = itarget

    seq_mirna = imirna  # [::-1]
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
        mim.libs = mirmap.library_link.LibraryLink(
            os.path.join(ivars["mirmap_script_dir"],
                         "libs/lib-archlinux-x86_64")
        )
        mim.exe_path = os.path.join(
            ivars["mirmap_script_dir"], "libs/exe-archlinux-x86_64"
        )
        # mim.eval_dg_open()  # Delta G open
        if mim.dg_duplex < threshold:
            # mim.eval_prob_binomial()  # P binomial
            return mim
        else:
            return None


def drawAlign(text, outfile, font, size):
    tk.Frame().destroy()
    txt = tkFont.Font(family=font, size=size)
    ilist = text.split("\n")
    list_width = []
    for i in ilist:
        iwidth = int(txt.measure(i))
        list_width.append(iwidth)

    width = int(max(list_width) * 0.9)
    height = int(width / 2.4)
    font = ImageFont.truetype(font.replace(" ", "-") + ".ttf",
                              size, encoding="utf-8")
    img = Image.new(
        "RGBA", (width, height), (255, 255, 255, 255)  # Dimensions  # Bg color
    )
    draw = ImageDraw.Draw(img)
    draw.text((20, 20), text, (0, 0, 0), font=font)
    # draw = ImageDraw.Draw(img)
    img.save(outfile)


# Alignment threshold
threshold = -13.8

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
ifile = open("Env_variables/Degradome_" + args.base_name + ".txt", "r")
ivars = defaultdict(str)
for line in ifile:
    ivar = line.strip().split("=")
    ivars[ivar[0].strip()] = ivar[1].strip()

# Add mirmap path to sys
sys.path.append(ivars["mirmap_script_dir"])
import mirmap
import mirmap.library_link

# Image settings
img_font = "FreeMono"  # "Hack Regular"
font_size = 20

# Directory definition
peakSeqDir = os.path.join(ivars["supp_data_dir"], "Peak_sequences")
mirSeqDir = os.path.join(ivars["supp_data_dir"], "miRNA_seq")
out_dir = os.path.join(ivars["output_dirR"], "03-Report/Summary")
out_dir_img_root = os.path.join(ivars["output_dirR"],
                                "03-Report/Alignment")

dir_exist(out_dir_img_root)
out_dir_img = os.path.join(out_dir_img_root, "mirmap")
dir_exist(out_dir_img)

# Define PyDegradome settings from variable definition file
pydeg_settings = np.array(
    ivars["pydeg_script_settings"].replace('"', "").strip(")(").split(" ")
)
isel = list(
    itertools.chain.from_iterable(
        itertools.repeat([True, False], int(len(pydeg_settings) / 2))
    )
)
isel_not = [not i for i in isel]
iMF_list = [int(x) for x in pydeg_settings[np.array(isel_not)]]
iConf_list = [float(x) for x in pydeg_settings[np.array(isel)]]


# Parse fasta files
mirSeq = SeqIO.parse(os.path.join(mirSeqDir,
                                  "input/miRNA_sequences.fa"),
                     "fasta")

mirSequences = []
for iseq in mirSeq:
    mirSequences.append(iseq)

for indx in range(len(iMF_list)):
    iMF = iMF_list[indx]
    iConf = iConf_list[indx]
    iConf2 = str(iConf).replace(".", "_")

    # Output file names
    out_obj = os.path.join(
        mirSeqDir, "output/mirmap_indices-" + str(iMF) +
        "_iConf-" + iConf2 + ".obj"
    )
    out_file_list = os.path.join(
        mirSeqDir,
        "output/mirmap_miRNA_targets_MF-" + str(iMF) +
        "_iConf-" + iConf2 + ".txt",
    )

    # Check if output files exist
    if os.path.exists(out_obj) and os.path.exists(out_file_list):
        print(
            "Files,\n- {os.path.basename(out_obj)}\
\n- {os.path.basename(out_file)}\nAlready exists"
        )
    else:
        ifasta = os.path.join(
            peakSeqDir, "PeakRegioncDNA_" +
            iConf2 + "_4_" + str(iMF) + ".fa"
        )
        if os.path.exists(ifasta):
            peakSequences = SeqIO.parse(ifasta, "fasta")

            peakRegions = []
            peakHeader = []
            for iseq in peakSequences:
                peakRegions.append(iseq)
                peakHeader.append(iseq.description)

            peakIndxFull = []
            for ihead in peakHeader:
                tmp = ihead.split()
                itx = tmp[0].replace(".", "_")
                icomp = tmp[1]
                ipeak = tmp[2].split(":")[1]
                peakIndxFull.append(icomp + itx + ipeak)

            if os.path.exists(out_obj):
                fileobj = open(out_obj, "rb")
                mm_indx_all = pickle.load(fileobj)
                fileobj.close()
            else:
                # Align
                # pool = mp.Pool(mp.cpu_count()-1)
                irange = range(len(mirSequences))
                mm_indx_all = []
                for peak_indx in range(len(peakRegions)):
                    ipeak = str(peakRegions[peak_indx].seq).replace("T", "U")

                    pool = mp.Pool(mp.cpu_count() - 1)
                    mm_indices = [
                        pool.apply(
                            f_mirmapindx,
                            args=(
                                ipeak,
                                str(mirSequences[mirindx].seq).replace("T",
                                                                       "U"),
                                mirindx,
                                peak_indx,
                            ),
                        )
                        for mirindx in range(len(mirSequences))
                    ]
                    pool.close()
                    mm_indices_filtered = [i for i in mm_indices if
                                           i is not None]
                    if len(mm_indices_filtered) > 0:
                        mm_indx_all.append(mm_indices_filtered)

                # Save indices
                if len(mm_indx_all) > 0:
                    fileObj = open(out_obj, "wb")
                    pickle.dump(mm_indx_all, fileObj)
                    fileObj.close()

            # If there are aligned sequences
            mirmap_final = []
            for indx in range(len(mm_indx_all)):
                ilist = mm_indx_all[indx]
                for item in ilist:
                    mir_indx = item[0]
                    peak_indx = item[1]
                    peak_indxFull = peakIndxFull[peak_indx]
                    ipeak = str(peakRegions[peak_indx].seq).replace("T", "U")
                    imirna = str(mirSequences[mir_indx].seq).replace("T", "U")
                    imirmap = f_mirmap(ipeak, imirna, threshold)
                    if imirmap is not None:
                        mirmap_final.append(
                            {
                                "imm": imirmap,
                                "miRNA_indx": mir_indx,
                                "peak_indx": peak_indx,
                                "peak_indxFull": peak_indxFull,
                            }
                        )

            if len(mirmap_final) > 0:
                # Output summary and aligment images
                file = open(out_file_list, "w")
                file.write(
                    "# List of pydegradome identified targets whose peak region align with a miRNA"
                )
                file.close()
                file = open(out_file_list, "a")
                for idic in mirmap_final:
                    imm = idic["imm"]
                    idg_duplex = imm.dg_duplex
                    mir_indx = idic["miRNA_indx"]
                    peak_indx = idic["peak_indx"]
                    peak_indxFull = idic["peak_indxFull"]
                    icomp = peakRegions[peak_indx].description.split(" ")[1]
                    itx = peakRegions[peak_indx].name
                    imir = mirSequences[mir_indx].name
                    report_split = imm.report().split("\n")

                    # Format delta symbol
                    for i in range(len(report_split)):
                        item = report_split[i]
                        if "\xce\x94" in item:
                            temp = item.replace("\xce\x94", "\u0394").decode(
                                "unicode-escape"
                            )
                            report_split[i] = temp

                    # Add orientation and align
                    if len(itx) <= len(imir):
                        isep = " " * (len(imir) + 4)
                        ispc = " " * (len(imir) - len(itx))
                        pre1 = itx + ispc + " 5' "
                        pre2 = imir + " 3' "
                    else:
                        isep = " " * (len(imir) + 4)
                        ispc = " " * (len(itx) - len(imir))
                        pre1 = itx + " 5' "
                        pre2 = imir + ispc + " 3' "

                    ireport = (
                        isep
                        + report_split[0]
                        + "\n"
                        + isep
                        + report_split[1]
                        + "\n"
                        + pre1
                        + report_split[2]
                        + " 3' "
                        + "\n"
                        + isep
                        + report_split[3]
                        + "\n"
                        + pre2
                        + report_split[4]
                        + " 5' "
                        + "\n\n"
                        + report_split[5]
                        + "\n"
                        + report_split[6]
                    )

                    # Save alignment as image
                    out_img = (
                        "_".join(
                            [
                                "Aln_mirmap",
                                itx.replace(".", "_"),
                                imir,
                                iConf2,
                                str(iMF),
                            ]
                        )
                        + ".png"
                    )
                    comp_dir = os.path.join(out_dir_img, icomp)
                    dir_exist(comp_dir)
                    img_path = os.path.join(comp_dir, out_img)
                    drawAlign(ireport, img_path, img_font, font_size)

                    # Write report
                    print >> file, (
                        "\nComparison: {0}\nTranscript: {1}\nmiRNA: {2}\nScore: {3}\nIndex: {4}".format(
                            icomp, itx, imir, idg_duplex, peak_indxFull
                        )
                    )
                file.close()
