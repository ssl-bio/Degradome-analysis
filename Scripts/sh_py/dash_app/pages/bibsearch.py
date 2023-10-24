import re
import plotly.graph_objects as go
from dash import dcc, html
import dash_dangerously_set_inner_html
import dash_bootstrap_components as dbc
import numpy as np
import pandas as pd
from Bio import Entrez
import nbib
import json

# Custom colors
mybg = '#ddeff3'
my_page_bg = 'whitesmoke'

# Icons
elink = '<i class="fa fa-arrow-up-right-from-square text-primary"></i>'

# pattern for transcript matching
tx_pattern = r'\bAT\w*\.[0-9]'

# Dictionary for renaming columns
new_columns = {
    "category_1": "Classification 1",
    "category_2": "Classification 2",
    "tx_name": "Transcript",
    "rep_gene": "Is representative",
    "feature_type": "Feature type",
    "gene_name": "Gene name",
    "ratioPTx": "Peak(T):Max(C)",
    "MorePeaks": "Peaks 2>",
    "Description": "Description",
    "comparison": "Comparison",
    "max_peak_1": "Max Peak (rep. 1)",
    "max_np_gene_1": "Max NonPeak (rep. 1)",
    "max_non_peak_ratio_1": "Peak:NonPeak (rep. 1)",
    "max_peak_2": "Max Peak (rep. 2)",
    "max_np_gene_2": "Max NonPeak (rep. 2)",
    "max_non_peak_ratio_2": "Peak:NonPeak (rep. 2)",
    "shared": "Shared (PyDegradome)",
    "max_peak_scaled": "Mean Peak Test (Scaled)",
    "max_read_tx_scaled": "Max Tx Ctrl (Mean Scaled)",
    "chr": "Chr",
    "strand": "Strand",
    "max_peak_test": "Max peak coor.",
    "peak_start": "Peak start",
    "peak_stop": "Peak stop",
    "gene_region_start": "Gene region start",
    "gene_region_end": "Gene region end",
    "feature_start": "Feature start",
    "feature_end": "Feature end",
    "tx_width": "Transcript len.",
    "cds_len": "CDS length",
    "width": "Feature width",
    "img_rank": "img_rank",
    "pydeg_settings": "pyDegradome settings",
}

selected_cols = [
    'tx_name',
    'chr',  # Level 2
    'feature_type',
    'comparison',
    'rep_gene',  # Level 3
    'strand',
    'shared',
    'category_1',
    'category_2',
    'MorePeaks',
    'max_peak_scaled',  # Quantitative
    'max_read_tx_scaled',
    'ratioPTx',
    'img_rank',  # Not used
    'gene_name',
    'Description',
    'pydeg_settings',  # Level 1
    'pydeg_settings_label',
    'peak_plot_link',
    'gene_plot_link',
    'miRNA_link',
    'plot_link',
    'id'
]

pytable_dropdown_default = dbc.Row([
    dbc.Col([
        dcc.Dropdown(
            id="class1_drop",
            placeholder='Classification 1',
            options=[0, 1, 2, 3, 4],
            className="dropdownFont"
        )], className='five.columns'),
    dbc.Col([
        dcc.Dropdown(
            id="class2_drop",
            placeholder='Classification 2',
            options=['A', 'B', 'C'],
            className="dropdownFont"
        )], className='five.columns'),
    dbc.Col([
        dcc.Dropdown(
            id="feat_drop",
            placeholder='Feature',
            options=['3UTR', '5UTR', 'CDS'],
            className="dropdownFont"
        )], className='five.columns'),
    dbc.Col([
        dcc.Dropdown(
            id="plot_drop",
            placeholder='Has peak plot',
            options=['No', 'Yes'],
            className="dropdownFont"
        )], className='five.columns'),
    dbc.Col([
        dcc.Dropdown(
            id="mirna_drop",
            placeholder='Has miRNA alignment',
            options=['No', 'Yes'],
            className="dropdownFont"
        )], className='five.columns')
])

pydeg_group_description = pd.read_csv(
    './data/PostPydeg_factor_description.tsv', sep='\t'
)


def import_vars(base):
    # Import database variables
    filename = f'./data/{base}_local_vars.json'
    with open(filename) as f:
        ivars = json.load(f)

    return ivars


def import_data(base, ivars, py_settings_str):
    base = ivars['ibase']
    py_settings=int(py_settings_str)
    
    # Import miRNA alignment dataframe
    miRNA_df = pd.read_csv(
        f"./data/miRNA_alignment_global_mirmap_{base}.tsv",
        sep='\t',
        low_memory=False,
    )
    miRNA_df = miRNA_df.loc[(miRNA_df['pydeg_settings'].eq(py_settings))]
    miRNA_df['id'] = range(0, len(miRNA_df))
    miRNA_df.loc[:, 'Score_y'] = miRNA_df['Score_y'].\
        round(2)
    miRNA_df.loc[:, 'Comparison'] = miRNA_df['Comparison'].\
        replace(ivars['comparison_dict'])

    # Import and process pydegradome data
    pydeg_main = pd.read_csv(
        f"./data/Candidate_peaks_degradome_{base}.tsv",
        sep='\t', low_memory=False
    )
    pydeg_main = pydeg_main.loc[(pydeg_main['pydeg_settings'].eq(py_settings))]
    pydeg_main['id'] = range(0, len(pydeg_main))
    pydeg_df = pydeg_main[selected_cols]
    pydeg_df.loc[:, 'ratioPTx'] = pydeg_df['ratioPTx'].round(2)
    pydeg_df.loc[:, 'comparison'] = pydeg_df['comparison'].\
        replace(ivars['comparison_dict'])

    return {"pydeg_df": pydeg_df, "miRNA_df": miRNA_df}


def toggle_show(n, is_open):
    # Class definition for headers
    icon_hide = "fas fa-angle-up me-4 align-self-center icon"
    icon_show = "fas fa-angle-down me-4 align-self-center icon"
    hdr_hide = 'd-flex bg-opacity-25 bg-primary'
    hdr_show = 'd-flex bg-opacity-10 bg-primary'
    # hdr_div = 'd-flex bg-opacity-25 bg-primary'
    if n:
        hide_show = not is_open
        if hide_show:
            btn = icon_hide
            hdr = hdr_hide
        else:
            btn = icon_show
            hdr = hdr_show
    else:
        hide_show = is_open
        btn = icon_hide

    return hide_show, btn, hdr


def calcFontSize(width):
    Wmin = 360
    Wmax = 2048
    Fmax = 30
    Fmin = 13
    fontSize = ((width - Wmin) / (Wmax - Wmin)) * (Fmax - Fmin) + Fmin
    return fontSize


def make_empty_fig():
    fig = go.Figure()
    fig.layout.paper_bgcolor = mybg
    fig.layout.plot_bgcolor = mybg
    return fig


def draw_miRNAplot(mirmap_link, globalAln_link):
    if isinstance(mirmap_link, str):
        miRNA_alignment_plot = html.Div(
            [html.P('Peak alignment (mirmap)'),
             html.Img(src=mirmap_link,
                      style={'width': '100%',
                             'height': 'auto'})],
            className='miRNA_alignment')
    elif isinstance(globalAln_link, str):
        miRNA_alignment_plot = html.Div(
            [html.P('Peak region alignment (global)'),
             html.Img(src=globalAln_link,
                      style={'width': '100%',
                             'height': 'auto'})],
            className='miRNA_alignment')
    return miRNA_alignment_plot


def get_description(group, category=None, html=True):
    if category:
        text_description = pydeg_group_description[
            (pydeg_group_description["group"] == group)
            & (pydeg_group_description["category"] == category)
        ]["description"]
    else:
        text_description = pydeg_group_description[
            pydeg_group_description["group"] == group
        ]["description"]
    if html:
        description = dash_dangerously_set_inner_html.DangerouslySetInnerHTML(
            text_description
        )
    else:
        description = text_description
    return description


def setEmail(email):
    Entrez.email = email
    return None


def concat_keys(dict, key_1, sep, key_2=None, last_sep=None):
    keys = ""
    ilen = len(dict[key_1])
    for j in range(ilen):
        if key_2:
            key = dict[key_1][j][key_2]
        else:
            key = dict[key_1][j]
        if j == 0:
            keys += key
        elif j < ilen - 1:
            keys += sep + " " + key
        else:
            if last_sep:
                keys += last_sep + " " + key
            else:
                keys += sep + " " + key
    return keys


def author_etal(auth_str):
    auth_list = re.split(r"[;&]", auth_str)
    if len(auth_list) > 2:
        authors = f"{auth_list[0]}; et. al"
        return authors
    else:
        return auth_str


def esearch(db, search_term, n_items=5):
    handle = Entrez.esearch(db=db, term=search_term, retmax=n_items)
    search_results = Entrez.read(handle)
    handle.close()
    if len(search_results) > 1:
        return search_results["IdList"]
    else:
        return None


def getPubmedId(IdList):
    handle = Entrez.efetch(db="pmc", id=IdList,
                           retmode="text", rettype="medline")
    biblio_data = handle.read()
    handle.close()
    records = nbib.read(biblio_data)

    # list ids of references
    ref_ids = [str(item["pubmed_id"]) for item in records]
    return ref_ids


def fetchRefs(db, idList):
    handle = Entrez.efetch(db=db, id=idList,
                           retmode="text", rettype="medline")
    biblio_data = handle.read()
    handle.close()

    return biblio_data


def getRefRecords(db, idList,
                  ret_mode="text", ret_type="medline"):
    # Fetch data
    handle = Entrez.efetch(db=db, id=idList,
                           retmode=ret_mode, rettype=ret_type)
    biblio_data = handle.read()
    handle.close()
    records = nbib.read(biblio_data)

    return records


def getBibDF(entries, transcript):
    keys = ["title", "authors", "journal",
            "publication_date", "doi", "abstract"]
    authors_list = []
    title_list = []
    publication_date_list = []
    journal_list = []
    doi_list = []
    transcript_list = [transcript] * len(entries)
    abstract_list = []

    counter = 1
    for entrie in entries:
        for key in keys:
            if key in entrie.keys():
                if key == "authors":
                    authors = concat_keys(
                        dict=entrie,
                        key_1="authors",
                        key_2="author",
                        sep=";",
                        last_sep=" &",
                    )
                    icolumn = author_etal(authors)
                elif key == "publication_date":
                    icolumn = entrie[key][:4]
                elif key == "doi":
                    icolumn = f'<a href="https://doi.org/{entrie[key]}"\
                    target=" blank"> {elink} </a>'
                else:
                    icolumn = entrie[key]
            else:
                icolumn = np.nan

            list_name = key + "_list"
            ilist = locals()[list_name]
            ilist.append(icolumn)

    data = {
        "Author(s)": authors_list,
        "Title": title_list,
        "Year": publication_date_list,
        "Journal": journal_list,
        "Transcript": transcript_list,
        "Doi": doi_list,
        "Abstract": abstract_list
    }

    return pd.DataFrame(data)


def printBibSection(iter_query, iter_ref, sep):
    bib_items = []
    for key in iter_ref:
        if key in iter_query.keys():
            if key == "authors":
                authors = concat_keys(
                    dict=iter_query,
                    key_1="authors",
                    key_2="author",
                    sep=";",
                    last_sep=" &",
                )
                bib_items.append(authors)
            elif key == "keywords":
                keywords = concat_keys(dict=iter_query,
                                       key_1="keywords", sep=";")
                bib_items.append(f"keywords: {keywords}")
            elif key == "doi":
                bib_items.append(f"doi: {iter_query[key]}")
            elif key == "abstract":
                bib_items.append(f"Abstract:\n{iter_query[key]}")
            else:
                bib_items.append(iter_query[key])
    return sep.join(bib_items)
