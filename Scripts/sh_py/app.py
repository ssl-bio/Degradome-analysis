import re
import dash
from dash import html, dcc, dash_table, Input, Output, callback, no_update
import dash_bootstrap_components as dbc
import dash_dangerously_set_inner_html
import plotly.graph_objects as go
import plotly.express as px
import pandas as pd
import json


app = dash.Dash(__name__,
                # external_scripts=[jquery_lib],
                external_stylesheets=[dbc.themes.COSMO])

# Import variables
filename = './data/local_vars.json'
with open(filename) as f:
    ivars = json.load(f)

ivars['cat1_dict'] = {int(k): v for k, v in ivars['cat1_dict'].items()}

# Code for formating data
# Data construction
pydeg_main = pd.read_csv(
    f"./data/Candidate_peaks_degradome_{ivars['ibase']}.tsv",
    sep='\t', low_memory=False
)
summary = pd.read_csv(
    f"./data/Peak_classification_summary_{ivars['ibase']}.tsv", sep='\t'
)
miRNA_alignment_df = pd.read_csv(
    f"./data/miRNA_alignment_global_mirmap_{ivars['ibase']}.tsv",
    sep='\t',
    low_memory=False,
)
pydeg_group_description = pd.read_csv(
    './data/PostPydeg_factor_description.tsv', sep='\t'
)

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
    'plot_link'
]

# pattern for transcript matching
tx_pattern = r'\bAT\w*\.[0-9]'

# filter columns
pydeg_df = pydeg_main[selected_cols]

# Create id column
pydeg_df.loc[:, 'id'] = pydeg_df.index
miRNA_alignment_df.loc[:, 'id'] = miRNA_alignment_df.index

# round values
pydeg_df.loc[:, 'ratioPTx'] = pydeg_df['ratioPTx'].round(2)
miRNA_alignment_df.loc[:, 'Score_y'] = miRNA_alignment_df['Score_y'].round(2)
ratioPTx_range = [int(pydeg_df['ratioPTx'].min()),
                  int(pydeg_df['ratioPTx'].max()) + 1]

# Treat numeric column as categorical
pydeg_df.loc[:, 'category_1'] = pydeg_df['category_1'].astype('category')

# Replace codes for values
pydeg_df.loc[:, 'comparison'] = pydeg_df['comparison'].\
    replace(ivars['comparison_dict'])
pydeg_df.loc[:, 'category_1'] = pydeg_df['category_1'].\
    replace(ivars['cat1_dict'])
pydeg_df.loc[:, 'category_2'] = pydeg_df['category_2'].\
    replace(ivars['cat2_dict'])

miRNA_alignment_df.loc[:, 'Comparison'] = miRNA_alignment_df['Comparison'].\
    replace(ivars['comparison_dict'])

# Define active cell in peak data frame.
# Depends on order of columns of pydeg_table
initial_active_cell = {'row': 0, 'column': 3, 'row_id': 0}

# header
row0 = initial_active_cell['row_id']
tx_name_default = pydeg_df.loc[row0,
                               'tx_name']
comparison_default = pydeg_df.loc[row0,
                                  'comparison']
category_default = pydeg_df.at[row0, 'category_1'].replace("Category ", "") +\
    "-" + pydeg_df.at[row0, 'category_2'].replace("Category ", "")
plot_header = f'Decay plot for {tx_name_default} ({category_default}) from \
comparison {comparison_default}'
row_default = 0


# Define active cell in peak data frame.
initial_active_cell_miRNA = {'row': 0, 'column': 0, 'row_id': 0}

mybg = '#ddeff3'
my_page_bg = 'whitesmoke'


def make_empty_fig():
    fig = go.Figure()
    fig.layout.paper_bgcolor = mybg
    fig.layout.plot_bgcolor = mybg
    return fig


def multiline_indicator(indicator):
    split = indicator.split()
    return '<br>'.join(split)


def get_description(group, author=None):
    if author:
        text_description = pydeg_group_description[
            (pydeg_group_description[
                'group'] == group) & (pydeg_group_description[
                    'category'] == author)]['description']
    else:
        text_description = pydeg_group_description[
            pydeg_group_description[
                'group'] == group]['description']
    description = dash_dangerously_set_inner_html.DangerouslySetInnerHTML(
        text_description
    )
    return description


app.layout = html.Div([
    #  [S] Title
    dbc.Row([
        dbc.Col(lg=1),
        dbc.Col([
            html.H1(ivars['doc_title']),
        ]),
        dbc.Col(lg=1),
    ]),
    html.Div(className='vspace_title'),
    # PyDegradome explanantion
    dbc.Row([
        dbc.Col(lg=1),
        dbc.Col([
            html.H2('Introduction'),
            html.H3('PyDegradome description'),
            html.Div(
                get_description('pydeg_description'),
                className='description_h3', id='pydeg_description'),
            html.Div(className='vspace_title'),
            html.H3('Test case description'),
            html.Div(
                get_description('test_case_description', ivars['ibase']),
                className='description_h3', id='study_description')
        ]),
        dbc.Col(lg=1)
    ]),
    html.Div(className='vspace_section'),
    #  [S] Classification post PyDegradome
    dbc.Row([
        dbc.Col(lg=1),
        dbc.Col([
            html.H2('Classification post PyDegradome'),
            html.Div(
                get_description('postpydeg_description'),
                className='description_h3', id='postpydeg_description')
        ]),
        dbc.Col(lg=1)
    ]),
    html.Div(className='vspace_section'),
    #  [S] PyDegradome settings
    dbc.Row([
        dbc.Col(lg=1),
        dbc.Col([
            html.H2('PyDegradome settings'),
        ]),
        dbc.Col(lg=1)
    ]),

    # pydeg settings slider
    dbc.Row([
        dbc.Col(lg=1),
        dbc.Col([
            dcc.RadioItems(
                id='pydeg_settings_item',
                options=[
                    {'label': [
                        html.Span(ivars['pydeg_settings']['0'])
                    ],
                     'value': 0},
                    {'label': [
                        html.Span(ivars['pydeg_settings']['1'])
                    ], 'value': 1},
                    {'label': [
                        html.Span(ivars['pydeg_settings']['2'])
                    ], 'value': 2},
                ],
                value=0, className='radio_items'
            )
        ], lg=2),
        dbc.Col([
            html.Div(
                get_description('pydeg_settings'),
                className='description_h3', id='pydeg_settings')
        ], lg=8),
    ]),
    html.Div(className='vspace_section'),
    #  [S] Summary of peak classification
    dbc.Row([
        dbc.Col(lg=1),
        dbc.Col([
            html.H2('Summary of peak classification'),
        ]),
        dbc.Col(lg=1)
    ]),
    # Plot group settings
    dbc.Row([
        dbc.Col(lg=1),
        dbc.Col([
            dbc.Label('Group value for X-axis',
                      className='description_h3'),
            dcc.Dropdown(
                id='x_axis_value',
                value='feature_type',
                placeholder='Value for x axis',
                options=[
                    {'label': new_columns[x_axis],
                     'value': x_axis} for x_axis in
                    selected_cols[1:4]
                ],
                className='description_h4'
            ),
        ]),
        dbc.Col([
            dbc.Label('Factor to further divide data',
                      className='description_h3'),
            dcc.Dropdown(
                id='group_value',
                value='shared',
                placeholder='Select a value divide bars',
                options=[
                    {'label': new_columns[group],
                     'value': group} for group in
                    selected_cols[4:10]
                ],
                className='description_h4'
            ),
        ]),
        dbc.Col(lg=1)
    ]),
    # summary plot
    dbc.Row([
        dbc.Col(lg=1),
        dbc.Col(dcc.Graph(id='peak_count_barplot',
                          figure=make_empty_fig())),
        dbc.Col(lg=1),
    ]),
    # Factors description
    dbc.Row([
        dbc.Col(lg=1),
        dbc.Col([
            html.Div(
                id='factor_description',
                className='description_h3')
        ]),
        dbc.Col(lg=1)
    ]),
    html.Div(className='vspace_subsection'),
    dbc.Row([
        dbc.Col(lg=1),
        dbc.Col([
            html.H3('List of candidates'),
            html.Div(
                get_description('peak_list'),
                className='description_h3', id='peak_list')
        ])
    ]),
    # Table's columns dropdown menu
    dbc.Row([
        dbc.Col(lg=1),
        dbc.Col([
            class1_drop := dcc.Dropdown(
                placeholder='Classification 1',
                options=[x for x in
                         sorted(pydeg_df.category_1.unique())]
            )],
                width=2
                ),
        dbc.Col([
            class2_drop := dcc.Dropdown(
                placeholder='Classification 2',
                options=[x for x in
                         sorted(pydeg_df.category_2.unique())]
            )],
                width=2,
                ),
        dbc.Col([
            feat_drop := dcc.Dropdown(
                placeholder='Feature',
                options=[x for x in
                         sorted(pydeg_df.feature_type.unique())]
            )],
                width=2,
                ),
        dbc.Col([
            plot_drop := dcc.Dropdown(
                placeholder='Has peak plot',
                options=[x for x in
                         sorted(pydeg_df.plot_link.unique())]
            )],
                width=2,
                ),
        dbc.Col([
            mirna_drop := dcc.Dropdown(
                placeholder='Has miRNA alignment',
                options=[x for x in
                         sorted(pydeg_df.miRNA_link.unique())]
            )],
                width=2,
                ),
        dbc.Col(lg=1),
    ], justify='start'),
    # Table peak candidates
    dbc.Row([
        dbc.Col(lg=1),
        dbc.Col([
            pydeg_table := dash_table.DataTable(
                id='py_table',
                columns=[
                    {'name': 'Rank',
                     'id': 'img_rank'},
                    {'name': 'Class. 1',
                     'id': 'category_1'},
                    {'name': 'Class. 2',
                     'id': 'category_2'},
                    {'name': 'Transcript',
                     'id': 'tx_name'},
                    {'name': 'Repr.?',
                     'id': 'rep_gene'},
                    {'name': 'Feat. type',
                     'id': 'feature_type'},
                    {'name': 'Gene name',
                     'id': 'gene_name'},
                    {'name': 'Peak(T):Tx Max(C)',
                     'id': 'ratioPTx'},
                    {'name': 'Peaks 2>',
                     'id': 'MorePeaks'},
                    {'name': 'Comparison',
                     'id': 'comparison'},
                    {'name': 'Plot link',
                     'id': 'plot_link'},
                    {'name': 'Alignment link',
                     'id': 'miRNA_link'},
                ],
                data=pydeg_df.to_dict('records'),
                filter_action='native',
                page_size=5,
                hidden_columns=['comparison', 'plot_link', 'miRNA_link'],
                style_table={'overflowX': 'auto'},
                active_cell=initial_active_cell,
                style_cell_conditional=[
                    {'if': {'column_id': 'MorePeaks'},
                     'width': '10%'},
                    {'if': {'column_id': 'category_1'},
                     'width': '9%'},
                    {'if': {'column_id': 'category_2'},
                     'width': '9%'},
                    {'if': {'column_id': 'tx_name'},
                     'width': '15%'},
                    {'if': {'column_id': 'rep_gene'},
                     'width': '9%'},
                    {'if': {'column_id': 'feature_type'},
                     'width': '12%'},
                    {'if': {'column_id': 'ratioPTx'},
                     'width': '15%'},
                    {'if': {'column_id': 'MorePeaks'},
                     'width': '9%'}
                ],
                style_data_conditional=([
                    {'if': {
                        'filter_query': '{plot_link} = Yes',
                        'column_id': 'tx_name'
                    },
                     'backgroundColor': '#48a4ba',
                     'fontWeight': 'bold'},
                    {'if': {
                        'filter_query': '{miRNA_link} = Yes && \
{category_1} = "Category 1" && {category_2} = "Category A"',
                        'column_id': 'tx_name'
                    },
                     'color': 'crimson',
                     'fontWeight': 'bold'}
                ]),
                style_cell={'fontSize': '19px'},
            )]),
        dbc.Col(lg=1),
    ]),
    html.Div(className='vspace_subsection'),
    # [S] Decay plots
    dbc.Row([
        dbc.Col(lg=1),
        dbc.Col([
            html.H3(id='tx_name', children=plot_header)]),
        dbc.Col(lg=1)
    ]),
    # Plots
    dbc.Row([
        dbc.Col(lg=1, md=0),
        dbc.Col([
            # html.H4('Transcript level'),
            html.Div(id='gene_plot')], sm=12, md=5),
        dbc.Col([
            # html.H4('Peak region'),
            html.Div(id='peak_plot'),
        ], sm=12, md=5),
    ]),
    html.Div(className='vspace_section'),
    # [T] miRNA alignment
    dbc.Row([
        dbc.Col(lg=1),
        dbc.Col([
            html.H2('Alignment of sequences around peaks with known miRNA'),
            html.Div(
                get_description('mirna_alignment'),
                className='description_h3', id='mirna_alignment')
        ]),
        dbc.Col(lg=1)
    ]),
    # Table miRNA alignment and plot
    dbc.Row([
        dbc.Col(lg=1),
        dbc.Col([
            html.Div(
                id='miRNA_div'
            ),
            html.Div([
                miRNA_alignment_results := dash_table.DataTable(
                    id='miRNA_datatable',
                    columns=[
                        {'name': 'Transcript',
                         'id': 'Transcript'},
                        {'name': 'miRNA',
                         'id': 'miRNA'},
                        {'name': 'Peak region alignment',
                         'id': 'Score_x'},
                        {'name': 'Peak alignment',
                         'id': 'Score_y'},
                        {'name': 'Comparison',
                         'id': 'Comparison'}],
                    data=miRNA_alignment_df.to_dict('records'),
                    filter_action='native',
                    style_cell_conditional=[
                        {'if': {'column_id': 'Transcript'},
                         'width': '20%'},
                        {'if': {'column_id': 'miRNA'},
                         'width': '20%'},
                        {'if': {'column_id': 'Score_x'},
                         'width': '9%'},
                        {'if': {'column_id': 'Score_y'},
                         'width': '15%'},
                        {'if': {'column_id': 'Comparison'},
                         'width': '30%'}
                    ],
                    style_cell={'fontSize': '19px'},
                    active_cell=initial_active_cell_miRNA,
                    page_size=3

                )
            ], id='miRNA_table', style={'display': 'block'})
        ], lg=6),
        dbc.Col([
            html.Div(id='miRNA_alignment_plot')
        ], sm=12, md=4, lg=4),
        dbc.Col(lg=1)
    ]),
    html.Div(className='vspace_section'),
    # References
    dbc.Row([
        dbc.Col(lg=1),
        dbc.Col([
            html.H2('References'),
            html.Div(
                get_description('references'),
                className='csl-entry', id='references')
        ]),
        dbc.Col(lg=1)
    ]),
    dbc.Col(lg=1)
], style={'backgroundColor': my_page_bg})


# Call back functions
@app.callback(
    Output('peak_count_barplot', 'figure'),
    Input('pydeg_settings_item', 'value'),
    Input('x_axis_value', 'value'),
    Input('group_value', 'value'),
)
def peak_count_barplot(settings, x1, x2):
    df_setting = pydeg_df[pydeg_df['pydeg_settings'] == settings]
    if x2:
        counts = df_setting.groupby([x1, x2]).size()
        df_counts = counts.reset_index().rename(
            columns={'level_0': x1, 'level_1': x2, 0: 'Peak number'}
        )
        df_counts = df_counts.rename(columns=new_columns)
        x1_plot = new_columns[x1]
        x2_plot = new_columns[x2]
        fig = px.bar(df_counts, x=x1_plot, y='Peak number',
                     color=x2_plot, barmode='group',
                     color_discrete_sequence=px.colors.qualitative.Plotly,
                     text_auto=True)
    else:
        counts = df_setting.groupby([x1]).size()
        df_counts = counts.reset_index().rename(
            columns={'level_0': x1, 0: 'Peak number'}
        )
        df_counts = df_counts.rename(columns=new_columns)
        x1_plot = new_columns[x1]
        fig = px.bar(df_counts, x=x1_plot, y='Peak number', text_auto=True)

    fig.layout.paper_bgcolor = mybg
    fig.update_layout(
        xaxis=dict(title_font=dict(size=22), tickfont=dict(size=19)),
        yaxis=dict(title_font=dict(size=22), tickfont=dict(size=19)),
        legend=dict(title_font=dict(size=17), font=dict(size=17)),
        legend_title=multiline_indicator(x2_plot)
    )

    return fig


@callback(
    Output(pydeg_table, 'data'),
    Input('pydeg_settings_item', 'value'),
    Input(class1_drop, 'value'),
    Input(class2_drop, 'value'),
    Input(feat_drop, 'value'),
    Input(plot_drop, 'value'),
    Input(mirna_drop, 'value')
)
def update_dropdown_options(pydeg_settings, class_1,
                            class_2, feat, plot, mirna):
    df = pydeg_df.copy()
    dff = df[df['pydeg_settings'] == pydeg_settings]

    if class_1:
        dff = dff[dff.category_1 == class_1]
    if class_2:
        dff = dff[dff.category_2 == class_2]
    if feat:
        dff = dff[dff.feature_type == feat]
    if plot:
        dff = dff[dff.plot_link == plot]
    if mirna:
        dff = dff[dff.miRNA_link == mirna]

    return dff.to_dict('records')


@app.callback(
    Output('tx_name', 'children'),
    Output('gene_plot', 'children'),
    Output('peak_plot', 'children'),
    Input('py_table', 'active_cell'),
)
def cell_clicked(active_cell):
    if active_cell is None:
        return no_update

    row = active_cell['row_id']
    transcript = pydeg_df.at[row, 'tx_name']
    comparison = pydeg_df.at[row, 'comparison']
    category = pydeg_df.at[row, 'category_1'].replace("Category ", "") + "-" +\
        pydeg_df.at[row, 'category_2'].replace("Category ", "")
    header = f'Decay plots for {transcript} ({category}) \
from comparison {comparison}'
    gene_plot_link = pydeg_df.at[row, 'gene_plot_link']

    # Check if link to plot exists
    if isinstance(gene_plot_link, str):
        gene_plot = html.Img(src=gene_plot_link,
                             style={'width': '100%',
                                    'height': 'auto'},
                             **{'data-table': transcript})
    else:
        gene_plot = html.P(f'No plot was produced for transcript {transcript}',
                           className='alt_text')

    peak_plot_link = pydeg_df.at[row, 'peak_plot_link']
    if isinstance(peak_plot_link, str):
        peak_plot = html.Img(src=peak_plot_link,
                             style={'width': '100%',
                                    'height': 'auto'},
                             **{'data-table': transcript})
    else:
        peak_plot = html.P(f'No plot was produced for transcript {transcript}',
                           className='alt_text')

    return header, gene_plot, peak_plot


@app.callback(
    Output('miRNA_div', 'children'),
    Output(miRNA_alignment_results, 'data'),
    Output(miRNA_alignment_results, 'page_size'),
    Output('miRNA_table', 'style'),
    Input('pydeg_settings_item', 'value'),
    Input('tx_name', 'children')
)
def miRNA_result(pydeg_settings, itx):
    transcript = re.findall(tx_pattern, itx)[0]
    df = miRNA_alignment_df.copy()
    dff = df.loc[(df['pydeg_settings'].eq(pydeg_settings)) &
                 (df['Transcript'] == transcript)]
    page_size = len(dff)
    if page_size > 0:
        # Candidates exist
        miRNA_div = ' '
        miRNA_alignment_results = dff.to_dict('records')
        display_status = {'display': 'block'}
    else:
        # No candidates found
        miRNA_div = html.P(
            f'No alignment was obtained for the peak sequence {transcript}',
            className='alt_text')
        miRNA_alignment_results = dff.to_dict('records')
        display_status = {'display': 'none'}

    return miRNA_div, miRNA_alignment_results, page_size, display_status


@app.callback(
    Output('miRNA_alignment_plot', 'children'),
    Input(miRNA_alignment_results, 'active_cell'),
    Input(miRNA_alignment_results, 'page_size'),
    Input('tx_name', 'children')
)
def cell_clicked_miRNA(active_cell, page_size, itx):
    category = re.findall(r'[0-4]-[A-C]', itx)[0]
    if page_size == 0 or category != '1-A':
        miRNA_alignment_plot = ' '
    else:
        row = active_cell['row_id']
        transcript = re.findall(tx_pattern, itx)[0]
        transcript_mirna = miRNA_alignment_df.at[row, 'Transcript']

        # Check if transcript from pydeg is same as the first from miRNA
        # This is needed when a miRNA table is shown for the first time
        if transcript != transcript_mirna:
            ibool = miRNA_alignment_df['Transcript'].str.contains(transcript)
            row = miRNA_alignment_df[ibool].index[0]

        mirmap_link = miRNA_alignment_df.at[row, 'mirmap_link']
        globalAln_link = miRNA_alignment_df.at[row, 'global_link']

        if (isinstance(mirmap_link, str) & isinstance(globalAln_link, str)):
            miRNA_alignment_plot = html.Div(
                [html.P('Peak alignment (mirmap)'),
                 html.Img(src=mirmap_link,
                          style={'width': '100%',
                                 'height': 'auto'}),
                 html.Br(),
                 html.P('Peak region alignment (global)'),
                 html.Img(src=globalAln_link,
                          style={'width': '100%',
                                 'height': 'auto'})],
                className='miRNA_alignment')
        elif isinstance(mirmap_link, str):
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
        else:
            print('There is no plot')
            miRNA_alignment_plot = ''

    return miRNA_alignment_plot


@app.callback(
    Output('factor_description', 'children'),
    Input('x_axis_value', 'value'),
    Input('group_value', 'value')
)
def display_group_info(group1, group2):
    group1_description = pydeg_group_description[
        pydeg_group_description['group'].eq(group1)]
    group2_description = pydeg_group_description[
        pydeg_group_description['group'].eq(group2)]

    if group1_description.empty and group2_description.empty:
        idiv = ' '
    elif group1_description.empty:
        info = group2_description['description']
        idiv = dash_dangerously_set_inner_html.DangerouslySetInnerHTML(
            info)
    elif group2_description.empty:
        info = group1_description['description']
        idiv = dash_dangerously_set_inner_html.DangerouslySetInnerHTML(
            info)
    else:
        info1 = group1_description['description']
        info2 = group2_description['description']
        idiv = html.Div([
            dash_dangerously_set_inner_html.DangerouslySetInnerHTML(
                info1),
            html.Br(),
            dash_dangerously_set_inner_html.DangerouslySetInnerHTML(
                info2)])
    return idiv


if __name__ == '__main__':
    app.run_server()
