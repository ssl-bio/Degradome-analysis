# 06-Process_data.py

"""
Description: Script to summarize the tables generated during post pyDegradome
analysis and export them to be used in a dash app.
Tables of candidate peaks obtained for all pyDegradome settings tested
(e.g. Pooled_0_95_4_2) are combined into a single data frame.
The same is done for the summary of peaks obtained at each filtering step.
Similarly, for the list of sequences (peak and around) that could be
potentially targeted by miRNA, these are combined into a single dataframe

Execution:
python Scripts/sh_py/06-process_data.py -rd $(pwd) -bn $ibase
Where $(pwd) is the input root directory and
$ibase refers to the project base name (e.g. Zhang-2021)

Output: Three text files, located in:
Degradome-$ibase/output_02/04-Dash_app/data
- Candidate_peaks_degradome_$ibase
- Peak_classification_summary_$ibase
- miRNA_alignment_global_mirmap_$ibase
"""
import argparse  # Parse arguments
from collections import defaultdict
import os  # access os directories and files
import pandas as pd  # Get summary stats
import numpy as np
import multiprocessing as mp
from pdf2image import convert_from_path
import re
import json
import itertools
import shutil


# Create a list to store results from all tested settings
def process_miRNAData(files_list):
    df_list = []
    for ifile in files_list:
        # print(f'Working on {ifile}')
        df = pd.read_csv(
            os.path.join(miRNA_dir, ifile),
            comment="#", sep=":", header=None
        )

        # Remove white spaces
        df[1] = df[1].str.strip()

        # Create an index column to group per comparison
        df["ID"] = df.groupby(0).cumcount() + 1

        # Pivot table focusing on comparison (ID) and column 0
        df_wide = df.pivot(index="ID", columns=0)

        # Edit column names (rm top level)
        df_wide.columns = df_wide.columns.levels[1]

        # Remove old index column
        df_wide.drop("Index", axis=1, inplace=True)

        # Add pydeg settings column
        iMF = ifile.split("-")[1].split("_")[0]
        iMF = int(iMF)
        iconf = ifile.split("-")[2].split(".")[0]
        iconf = float(iconf.replace("_", "."))
        pydeg_settings_label = f"MF: {iMF} - conf: {iconf}"

        df_wide["pydeg_settings_label"] = pydeg_settings_label
        df_wide.index = range(0, len(df_wide))
        df_list.append(df_wide)

    df_merged = pd.concat(df_list)

    return df_merged


# Drop redundant columns from second dataframe
def merge_columns(df, col_a, col_b, new_col):
    merged_values = np.where(
        df[col_a].notnull(),
        df[col_a],
        np.where(df[col_b].notnull(), df[col_b], df[col_a]),
    )
    df[new_col] = merged_values
    df.drop(col_a, axis=1, inplace=True)
    df.drop(col_b, axis=1, inplace=True)
    return df


def set_plotLink(df, base_name):
    """
    Process the columns of a data frame to
    produce a series with links to plots
    """
    tx_list = [itx.replace(".", "_") for itx in df['tx_name']]
    classification = df['category_1'].astype(str).str.cat(df['category_2'],
                                                          sep='-')
    path_list = [df['img_rank'], ['Peak']*len(df), classification,
                 tx_list, [iconf.replace(".", "_")]*len(df),
                 ['4']*len(df), [iMF+".jpg"]*len(df)]
    path_list_t = list(map(list, zip(*path_list)))

    # Convert transposed list to DataFrame
    path_df = pd.DataFrame(path_list_t)
    end_path = path_df.astype(str).apply('_'.join, axis=1)
    base_path = f'assets/Dplots/{base_name}/Peak_' + df['comparison'] + "/"
    plot_path = pd.DataFrame([base_path, end_path]).apply("".join)
    return plot_path


def set_plotLink_mirna(df, align_type, base_name):
    """
    Process the columns of a data frame to
    produce a series with links to plots
    """
    df.index = range(0, len(df))
    tx_list = [itx.replace(".", "_") for itx in df['Transcript']]
    path_list = [[f"Aln_{align_type}"]*len(df),
                 tx_list, df['miRNA'],
                 [iconf.replace(".", "_")]*len(df),
                 [iMF+".png"]*len(df)]
    path_list_t = list(map(list, zip(*path_list)))

    # Convert transposed list to DataFrame
    path_df = pd.DataFrame(path_list_t)
    end_path = path_df.astype(str).apply('_'.join, axis=1)
    base_path = f'assets/Alignment/{base_name}/{align_type}/' +\
        df['Comparison'] + "/"
    plot_path = pd.DataFrame([base_path, end_path]).apply("".join)
    return plot_path


def plot2jpg(plot, indir, outdir, dpi=150):
    """
    Converts a pdf into a jpg file if the latter doesn't
    exist in the 'outdir'
    """
    in_file = os.path.join(indir, plot)
    plot_jpg = plot.replace(".pdf", ".jpg")
    out_file = os.path.join(outdir, plot_jpg)
    if ~os.path.isfile(out_file):
        img = convert_from_path(in_file, dpi=dpi)
        img[0].save(out_file, 'JPEG')


def samples2dict(names_str, values_str):
    names_list = names_str.replace('("', "").replace('")', "").split('" "')
    values_list = values_str.replace('("', "").replace('")', "").split('" "')
    idict = {k: v for k, v in zip(names_list, values_list)}
    return idict


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

# Define directories
pooledDir = os.path.join(ivars["output_dirR"], "02-PyDegradome_pooled")
summaryDir = os.path.join(ivars["output_dirR"], "03-Report", "Summary")
miRNA_dir = os.path.join(
    ivars["output_dir_base"], "Supporting_data", "miRNA_seq", "output"
)
output_dir = os.path.join(
    ivars["output_dirR"], '04-Dash_app'
)
dplot_path = os.path.join(ivars['output_dirR'], "03-Report", "Dplots")
dplot_path_out = os.path.join(output_dir, "assets",
                              "Dplots", ivars['ibase'])
dplot_dirs = os.listdir(dplot_path)
idirs = [f"assets/Dplots/{ivars['ibase']}/" + idir for idir in dplot_dirs]

# Create dash file structure
idirs += ['data', 'assets/Images']
for idir in idirs:
    ipath = os.path.join(output_dir, idir)
    if (not os.path.exists(ipath)):
        os.makedirs(ipath)

# Copy files to the dash app folder
shutil.copytree(os.path.join('Scripts/sh_py/dash_app'),
                output_dir, dirs_exist_ok=True)

shutil.copy(os.path.join('./Env_variables',
                         'PostPydeg_factor_description.tsv'),
            os.path.join(output_dir, 'data'))

shutil.copy(os.path.join('Env_variables', 'custom.css'),
            os.path.join(output_dir, 'assets'))

shutil.copytree(os.path.join('Env_variables', 'Images'),
                os.path.join(output_dir, 'assets/Images'),
                dirs_exist_ok=True)

# Convert plots from pdf to jpg
for idir in dplot_dirs:
    dplot_dir = os.path.join(dplot_path, idir)
    dplot_dir_out = os.path.join(dplot_path_out, idir)
    dplot_files = os.listdir(dplot_dir)
    dplot_files_out = os.listdir(dplot_dir_out)
    dplot_files_in = [iplot for iplot in dplot_files if
                      iplot.replace(".pdf", ".jpg") not in dplot_files_out]
    if dplot_files_in:
        irange = range(len(dplot_files_in))
        pool = mp.Pool(mp.cpu_count()-2)
        l_score = [pool.apply(plot2jpg,
                              args=(dplot_files_in[i],
                                    dplot_dir,
                                    dplot_dir_out)
                              ) for i in irange]
        pool.close()

# Process pooled files of candidate peaks
pooledFiles = [item for item in os.listdir(pooledDir)
               if re.search(r"^Pooled.*", item)]
pooled_list = []
pydeg_settings = {}
for i in range(0, len(pooledFiles)):
    in_file = os.path.join(pooledDir, pooledFiles[i])
    pooled_df = pd.read_csv(in_file, sep="\t")
    comparisons = pooled_df["comparison"].unique()

    # Add a column for plot number
    comparisons_list = []
    for icomp in comparisons:
        ibool = pooled_df["comparison"] == icomp
        df_comparison = pooled_df[ibool]
        # Treat category_1 as factor
        custom_order1 = [1, 2, 3, 4, 0]
        custom_order2 = ['A', 'B', 'C']
        df_comparison['category_1'] = df_comparison['category_1'].astype(
            pd.CategoricalDtype(categories=custom_order1,
                                ordered=True))
        # Separate each unique combination of categories 1 and 2
        # and rank plots based on the ratio of peak:transcript
        cat1 = pd.Series(df_comparison['category_1'].unique())
        cat2 = pd.Series(df_comparison['category_2'].unique())
        cat1 = cat1.sort_values(
            key=lambda x: x.map({k: i for i, k in enumerate(
                custom_order1)}))
        cat2 = cat2.sort_values(
            key=lambda x: x.map({k: i for i, k in enumerate(
                custom_order2)}))
        categories = list(itertools.product(cat1, cat2))

        for j in categories:
            cat1 = j[0]
            cat2 = j[1]
            ibool = (df_comparison['category_1'] == cat1) &\
                (df_comparison['category_2'] == cat2)
            df_cat = df_comparison[ibool]

            df_cat = df_cat.sort_values(by=["category_1",
                                            "category_2",
                                            "ratioPTx"],
                                        ascending=[True,
                                                   True,
                                                   False])
            df_cat["img_rank"] = [
                "{:02d}".format(k) for k in range(1, len(df_cat) + 1)
            ]
            comparisons_list.append(df_cat)
    pooled_df2 = pd.concat(comparisons_list)

    # Add data
    file_data_list = in_file.split("/")[-1].split("_")
    iconf = ".".join(file_data_list[1:3])
    iMF = file_data_list[-1]
    pydeg_settings_label = f"MF: {iMF} - conf: {iconf}"
    pooled_df2["pydeg_settings_label"] = pydeg_settings_label
    pooled_df2["pydeg_settings"] = i

    pydeg_settings[i] = pydeg_settings_label
    # Add link to plots
    pooled_df2.index = range(0, len(pooled_df2))

    # peak plot
    link_plot = set_plotLink(pooled_df2, ivars['ibase'])
    link_test = [os.path.isfile(os.path.join(output_dir, ilink)) for
                 ilink in link_plot]
    pooled_df2['peak_plot_link'] = np.where(link_test, link_plot, np.nan)

    # gene plot
    link_plot2 = [link.replace('Peak', 'Gene') for link in link_plot]
    link_test = [os.path.isfile(os.path.join(output_dir, ilink)) for
                 ilink in link_plot2]
    pooled_df2['gene_plot_link'] = np.where(link_test, link_plot2, np.nan)
    pooled_list.append(pooled_df2)

pydeg_df = pd.concat(pooled_list)
pydeg_df.index = range(0, len(pydeg_df))

# Define redundant columns
redundant_cols = [
    "ratio",
    "factor",
    "max_count_test_2",
    "max_count_ctl",
    "tot_count_test",
    "tot_count_ctl",
    "peak_width",
    "ID",
    "rep_gene_note",
    "max_count_test_1"
]

# Remove redundant columns
selected_cols = [col for col in pydeg_df.columns.to_list()
                 if col not in redundant_cols]
pydeg_df = pydeg_df[selected_cols]

# Change numeric values to strings or factors
# Dictionary for yes/no replacement
yes_no_mapping = {0: 'No', 1: 'Yes'}
pydeg_df['rep_gene'] = pydeg_df['rep_gene'].replace(yes_no_mapping)
pydeg_df['MorePeaks'] = pydeg_df['MorePeaks'].replace(yes_no_mapping)

# Dictionary for shared peaks btw replicates
shared_mapping = {1: 'Only Rep. 1', 2: 'Only Rep. 2', 3: 'Both Reps.'}
pydeg_df['shared'] = pydeg_df['shared'].replace(shared_mapping)

# Summary
ifile = "Peak_counts_BeforeAfter_Filtering"
sum_df = pd.read_csv(os.path.join(summaryDir, ifile), sep="\t")
isettings = sum_df["Settings"].tolist()
iMF_list, iconf_list = [], []
for i in isettings:
    isetting = i.split("-")
    iMF = isetting[0]
    iconf = isetting[1]
    iMF_list.append(iMF)
    iconf_list.append(iconf)

sum_df["MF"] = iMF_list
sum_df["conf"] = iconf_list
sum_df.drop("Settings", axis=1, inplace=True)

# Write data frame to file
out_file = os.path.join(
    output_dir, "data", "Peak_classification_summary_" +
    ivars["ibase"] + ".tsv")
sum_df.to_csv(out_file, sep='\t')

# miRNA
# List files from global alignment
miRNAFiles_global = [
    item for item in os.listdir(miRNA_dir)
    if re.search(r"^miRNA_targets.*", item)
]

miRNAFiles_mirmap = [
    item for item in os.listdir(miRNA_dir)
    if re.search(r"^mirmap.*.txt", item)
]

# Add link to existing plots
df_miRNAalignment_global = process_miRNAData(miRNAFiles_global)
link_plot = set_plotLink_mirna(df_miRNAalignment_global,
                               'global', ivars['ibase'])
link_test = [os.path.isfile(os.path.join(output_dir, ilink)) for
             ilink in link_plot]
df_miRNAalignment_global['global_link'] = np.where(link_test,
                                                   link_plot, np.nan)

df_miRNAalignment_mirmap = process_miRNAData(miRNAFiles_mirmap)
link_plot = set_plotLink_mirna(df_miRNAalignment_mirmap,
                               'mirmap', ivars['ibase'])
link_test = [os.path.isfile(os.path.join(output_dir, ilink)) for
             ilink in link_plot]
df_miRNAalignment_mirmap['mirmap_link'] = np.where(link_test,
                                                   link_plot, np.nan)

# Create an index column for merging
icols = ["Comparison", "pydeg_settings_label", "Transcript", "miRNA"]
df_miRNAalignment_global["index"] = (
    df_miRNAalignment_global[icols].astype(str).apply("_".join, axis=1)
)
df_miRNAalignment_mirmap["index"] = (
    df_miRNAalignment_mirmap[icols].astype(str).apply("_".join, axis=1)
)

df_miRNAalignment_all = pd.merge(
    df_miRNAalignment_global, df_miRNAalignment_mirmap, on="index", how="outer"
)

# Remove duplicated columns (ending in _x & _y)
for icol in icols:
    iregex = icol + ".*"
    df_col = df_miRNAalignment_all.filter(regex=iregex).columns.to_list()
    df_miRNAalignment_all = merge_columns(
        df_miRNAalignment_all, df_col[0], df_col[1],
        df_col[0].rsplit('_', 1)[0]  # split from right
    )

# Add a pydeg_settings (int) column
pydeg_settings_inv = {v: k for k, v in pydeg_settings.items()}
df_miRNAalignment_all['pydeg_settings'] = df_miRNAalignment_all[
    'pydeg_settings_label'].replace(pydeg_settings_inv)
df_miRNAalignment_all.drop("index", axis=1, inplace=True)

# Remove row with no links for both alignments
ibool = df_miRNAalignment_all['global_link'].isna() &\
    df_miRNAalignment_all['mirmap_link'].isna()
df_miRNAalignment_all = df_miRNAalignment_all[~ibool]

# Write data frame to file
out_file = os.path.join(
    output_dir, "data", "miRNA_alignment_global_mirmap_" +
    ivars["ibase"] + ".tsv")
df_miRNAalignment_all.to_csv(out_file, sep='\t', index=False)


# Create dummy columns for
# transcripts with links to plots and alignment
mirna_df = df_miRNAalignment_all[['Transcript',
                                  'pydeg_settings',
                                  'global_link',
                                  'mirmap_link']]
# remove rows with no links
no_align = mirna_df['global_link'].isna() & mirna_df['mirmap_link'].isna()
mirna_df = mirna_df[~no_align]

# create index column for merging
mirna_df['indx'] = mirna_df[['Transcript',
                             'pydeg_settings']].astype(str).apply(
                                 "_".join, axis=1)

# remove duplicated index
indx_dup = mirna_df["indx"].duplicated()
mirna_df = mirna_df[~indx_dup]

# create dummy column for links and remove unnecessary columns
mirna_df.drop(['Transcript', 'pydeg_settings',
               'global_link', 'mirmap_link'], axis=1, inplace=True)
mirna_df['miRNA_link'] = 'Yes'


# Create indx column for merging and dummy column for plots
pydeg_df['indx'] = pydeg_df[
    ['tx_name',
     'pydeg_settings']].astype(str).apply("_".join, axis=1)
plot_link_bool = pydeg_df['peak_plot_link'].isna() &\
    pydeg_df['peak_plot_link'].isna()
pydeg_df['plot_link'] = plot_link_bool.map({True: 'No', False: 'Yes'})

# merge and replace NaN values with 'No'
pydeg_df = pd.merge(pydeg_df, mirna_df, on='indx', how='left')
pydeg_df["miRNA_link"] = pydeg_df[
    "miRNA_link"].map({'Yes': 'Yes', np.nan: 'No'})

# Write list of peaks to file
out_file = os.path.join(
    output_dir, "data", "Candidate_peaks_degradome_" +
    ivars["ibase"] + ".tsv")
pydeg_df.to_csv(out_file, sep='\t')

# Add additional variables and export
ivars['pydeg_settings'] = pydeg_settings
ivars['pydeg_settings_inv'] = pydeg_settings_inv

# Dictionary of sample code and sample name (per replicate)
test_samples = samples2dict(ivars['test_samples'],
                            ivars['test_samples_name'])
control_samples = samples2dict(ivars['control_samples'],
                               ivars['control_samples_name'])
# Merge dictionaries
sample_dict = dict(test_samples.items() | control_samples.items())
ivars['sample_dict'] = sample_dict

# Dictionary of comparison code and name
comparison_dict = {}
icomp_list = pydeg_df.comparison.unique()
for icomp in icomp_list:
    isamples = icomp.replace("_and_", "-").split("-")
    isamples_names = (
        pd.Series(isamples)
        .replace(sample_dict)
        .str.replace(" \[[1-2]\]$", "", regex=True)
        .unique()
    )
    comparison = f"{isamples_names[0]} vs. {isamples_names[1]}"
    comparison_dict[icomp] = comparison
ivars['comparison_dict'] = comparison_dict

# Dictionary for classification
cat1_dict = {0: "Discarded",
             1: "Category 1",
             2: "Category 2",
             3: "Category 3",
             4: "Category 4"}

cat2_dict = {"A": "Category A",
             "B": "Category B",
             "C": "Category C"}

ivars['cat1_dict'] = cat1_dict
ivars['cat2_dict'] = cat2_dict

out_file = os.path.join(
    output_dir, "data", f"{ivars['ibase']}_local_vars.json")
with open(out_file, 'w') as out:
    json.dump(ivars, out)
