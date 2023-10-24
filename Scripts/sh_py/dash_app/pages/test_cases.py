import re
import dash
from dash import dcc, html, Input, Output, callback, \
    dash_table, no_update, State, clientside_callback
import dash_bootstrap_components as dbc
from dash_bootstrap_templates import ThemeChangerAIO, \
    template_from_url
import dash_dangerously_set_inner_html
import dash_breakpoints
import plotly.express as px
import pandas as pd
# import numpy as np
# import json
from . import bibsearch as bib

dash.register_page(__name__, name='Test cases')

# Define active cell in peak data frame.
initial_active_cell = {'row': 0, 'column': 3, 'row_id': 0}
initial_active_cell_bib = {'row': 0, 'column': 1}

# Class definition for headers
icon_show = "fas fa-angle-down me-4 align-self-center icon"
icon_hide = "fas fa-angle-up me-4 align-self-center icon"
hdr_hide = 'd-flex bg-opacity-25 bg-primary'
hdr_show = 'd-flex bg-opacity-10 bg-primary'
hdr_div = 'd-flex bg-opacity-25 bg-primary'

notes_btn = html.Div([
    html.I(className="fas fa-pencil-square text-primary fs-2 opacity-50",
           id="btn_notes",
           n_clicks=0
           ),
    dbc.Tooltip(
        "Click to take notes",
        target="btn_notes",
        placement='left',
        className="bg-opacity-50 bg-primary text-secondary"
    ),
    dbc.Offcanvas(
        id="offcanvas_notes",
        scrollable=True,
        backdrop=False,
        placement='bottom',
        title="Notes",
        is_open=False,
    )
], style={'position': 'fixed',
          'bottom': 0,
          'right': '3vw',
          'marginBottom': '30vh',
          })

dataSet_cards = dbc.CardGroup([
    dbc.Card(
        dbc.CardBody([
            html.H5(id="dataSet_title",
                    className="mb-3"),
            html.P(id="dataSet_authors",
                   className="description_h4 mb-2"),
            html.P(id="dataSet_date",
                   className="d-flex description_h4 justify-content-end"),
            dbc.Button("Selected",
                       id="zhang_btn",
                       n_clicks=0,
                       color="info",
                       outline=True),
        ]),
        id="zhang_card",
        color="primary",
        outline=True,
        className="bg-primary bg-opacity-75")
])

layout = html.Main([
    html.Div([
        html.Div(id='viewport-container', style={'display': 'none'}),
        dash_breakpoints.WindowBreakpoints(
            id="breakpoints",
            widthBreakpointThresholdsPx=[576, 768, 992, 1200, 1400],
            widthBreakpointNames=["xs", "sm", "md", "lg", "xl", "xxl"],
        )
    ]),
    html.Section([
        dbc.Alert([
            html.H2('Datasets',
                    className='border-bottom border-3 border-secondary pb-2 mb-4'),
            dataSet_cards
        ], color='primary'),
        dcc.Store(
            data='Zhang-2021',
            id='dataSet_name'
        ),
        dcc.Store(
            id="ivars"
        ),
        dcc.Store(
            id="pydeg_data"
        ),
        dcc.Store(
            id="miRNA_data"
        )
    ]),
    html.Section([
        dbc.Alert([
            html.Div([
                html.H2('Test case description',
                        className='flex-fill'),
                html.I(className=icon_hide,
                       id="description_btn",
                       n_clicks=0)
            ], id="description_header",
                     className=hdr_div),
            dbc.Collapse(
                is_open=True,
                id="description_html",
                className="mt-4",
            ),
        ], color='secondary')
    ]),

    html.Section([
        dbc.Alert([
            html.Div([
                html.H2('PyDegradome settings',
                        className='flex-fill'),
                html.I(className=icon_hide,
                       id="settings_btn",
                       n_clicks=0)
            ], id="settings_header",
                     className=hdr_div),
            dbc.Collapse(
                is_open=True,
                id="settings_description_html",
                className="mt-4",
                children=[
                    html.Div([
                        dbc.RadioItems(
                            id='pydeg_settings_item',
                            options=[
                                {'label': "MF: 2 - CI: 0.95", 'value': 0},
                                {'label': "MF: 4 - CI: 0.95", 'value': 1},
                                {'label': "MF: 3 - CI: 0.99", 'value': 2},
                            ],
                            value=0,
                            className='pysettings_items'
                        ),
                        html.Div(
                            bib.get_description('pydeg_settings'),
                            className='description_h3 pysettings_description',
                            id='pydeg_settings'
                        )], className="d-flex")
                ])
        ], color='secondary')
    ]),

    #  [S] Summary of peak classification
    html.Section([
        dbc.Alert([
            html.Div([
                html.H2('Summary of peak classification',
                        className='flex-fill'),
                html.I(className=icon_hide,
                       id="summary_btn",
                       n_clicks=0)
            ], id="summary_header",
                     className=hdr_div),
            dbc.Collapse(
                is_open=True,
                id="summary_description_html",
                className="mt-4",
                children=[
                    html.Div(
                        bib.get_description('peak_count'),
                        className='description_h3 mt-3 mb-4'
                    ),
                    dbc.Row([
                        dbc.Col([
                            dbc.Label('X-axis variable',
                                      className='description_h4'),
                            dcc.Dropdown(
                                id='x_axis_value',
                                value='comparison',
                                placeholder='Value for x axis',
                                options=[
                                    {'label': bib.new_columns[x_axis],
                                     'value': x_axis} for x_axis in
                                    bib.selected_cols[1:4]
                                ], className="dropdownFont"
                            )
                        ], className='four.columns'),
                        dbc.Col([
                            dbc.Label('Group data by',
                                      className='description_h4'),
                            dcc.Dropdown(
                                id='group_value',
                                value='category_1',
                                placeholder='Select a value divide bars',
                                options=[
                                    {'label': bib.new_columns[group],
                                     'value': group} for group in
                                    bib.selected_cols[4:10]
                                ], className="dropdownFont"
                            ),
                        ], className='four.columns'),
                        dbc.Col([
                            dbc.Label('Exclude from classification 1',
                                      className='description_h4'),
                            dcc.Dropdown(
                                id='class1_value',
                                placeholder='Remove one or more',
                                multi=True,
                                persistence=False,
                                className="dropdownFont"
                            ),
                        ], className='four.columns'),
                        dbc.Col([
                            dbc.Label('Remove comparison',
                                      className='description_h4'),
                            dcc.Dropdown(
                                id='comparison_value',
                                className="dropdownFont"
                            ),
                        ], className='four.columns')
                    ], className="mb-3"),
                    html.Div([
                        dcc.Graph(id='peak_count_barplot',
                                  figure=bib.make_empty_fig(),
                                  ),
                    ]),
                    html.Div(
                        id='factor_description',
                        className='mt-4 description_h3'
                    )
                ])
        ], color='secondary')
    ]),

    # # [S] List of candidates
    html.Section([
        dbc.Alert([
            html.Div([
                html.H3('List of candidates',
                        className='flex-fill'),
                html.I(className=icon_hide,
                       id="candidates_btn",
                       n_clicks=0)
            ], id="candidates_header",
                     className=hdr_div),
            dbc.Collapse(
                is_open=True,
                id="candidates_description_html",
                className="mt-4",
                children=[
                    html.Div(
                        bib.get_description('peak_list'),
                        className='description_h3',
                        id='peak_list'
                    ),
                    html.Div([
                        html.Div([
                            html.Div(id="dropdown_pytable",
                                     children=bib.pytable_dropdown_default),
                            html.Div([
                                dash_table.DataTable(
                                    id='py_table',
                                    columns=[
                                        {'name': ['', 'Rank'],
                                         'id': 'img_rank'},
                                        {'name': ['Classification', '1'],
                                         'id': 'category_1'},
                                        {'name': ['Classification', '2'],
                                         'id': 'category_2'},
                                        {'name': ['', 'Transcript'],
                                         'id': 'tx_name'},
                                        {'name': ['', 'Repr.?'],
                                         'id': 'rep_gene'},
                                        {'name': ['', 'Feature'],
                                         'id': 'feature_type'},
                                        {'name': ['', 'Gene name'],
                                         'id': 'gene_name'},
                                        {'name': ['Peak(T):', 'Tx Max(C)'],
                                         'id': 'ratioPTx'},
                                        {'name': ['Peaks', '2>'],
                                         'id': 'MorePeaks'},
                                        {'name': 'Comparison',
                                         'id': 'comparison'},
                                        {'name': 'Plot link',
                                         'id': 'plot_link'},
                                        {'name': 'Alignment link',
                                         'id': 'miRNA_link'},
                                    ],
                                    style_table={'overflowX': 'auto'},
                                    merge_duplicate_headers=True,
                                    page_size=15,
                                    page_current=0,
                                    row_selectable='multi',
                                    selected_rows=[],
                                    hidden_columns=['comparison',
                                                    'plot_link', 'miRNA_link'],
                                    active_cell=initial_active_cell,
                                    style_cell={
                                        'font-size': '0.575rem',
                                        'font-family': 'sans-serif',
                                        'padding': '0.5em 0.5em'
                                    },
                                    style_data_conditional=([
                                        {'if': {
                                            'filter_query':
                                            '{plot_link} = Yes',
                                            'column_id': 'tx_name'
                                        },
                                         'backgroundColor': 'var(--bs-light)',
                                         'fontWeight': 'bold'},
                                        {'if': {
                                            'filter_query':
                                            '{miRNA_link} = Yes',
                                            'column_id': 'tx_name'
                                        },
                                         'color': 'var(--bs-primary)',
                                         'fontWeight': 'bold'}
                                    ]),
                                    css=[{"selector": ".show-hide",
                                          "rule": "display: none"},
                                         {'selector':
                                          '.previous-next-container',
                                          'rule': 'font-size: 0.625rem;'}]
                                )
                            ]),  # Table
                        ], className='py_table'),
                        html.Div([
                            dcc.Tabs(
                                id='tab_plots',
                                value='gene_plot',
                                children=[
                                    dcc.Tab(
                                        label='Decay Plot [Gene level]',
                                        value='gene_plot',
                                        className="dropdownFont text-wrap"),
                                    dcc.Tab(
                                        label='Decay Plot [Peak level]',
                                        value='peak_plot',
                                        className="dropdownFont text-wrap"),
                                    dcc.Tab(
                                        label='miRNA alignment',
                                        value='miRNA_tab',
                                        className="dropdownFont text-wrap",
                                        children=[
                                            dash_table.DataTable(
                                                id='miRNA_datatable',
                                                selected_cells=[])
                                        ]
                                    )
                                ]
                            ),
                            html.Div(id='tabs-content-decay-plots')
                        ], className='plots'),
                        dcc.Store(
                            id='miRNAplot_df',
                            data=None
                        ),
                    ], className='py_results')
                ]),
        ], color='secondary')
    ]),  # Section
    # # [S] Related literature
    html.Section([
        dbc.Alert([
            html.Div([
                html.H2('Related literature for the selected transcripts',
                        className='flex-fill'),
                html.I(className=icon_hide,
                       id="search_refs_btn",
                       n_clicks=0)
            ], id="search_refs_header",
                     className=hdr_div),
            dbc.Collapse(
                is_open=True,
                id="search_refs_description_html",
                className="mt-4",
                children=[
                    notes_btn,
                    dbc.Row([
                        dbc.Col([
                            dbc.Label("NCBI valid e-mail",
                                      html_for="ncbi_email",
                                      className='description_h4'),
                            dbc.Input(
                                type="email",
                                id="ncbi_email",
                                placeholder="Enter email",
                                debounce=True,
                                persistence=True,
                                persistence_type='session',
                                required=True,
                                style={"fontSize": "0.750rem"}
                            ),
                        ],
                                width=4
                                ),
                        dbc.Col([
                            dbc.Label("Search term",
                                      html_for="op_term",
                                      className='description_h4'),
                            dbc.Input(
                                type="text",
                                id="op_term",
                                placeholder="Search term",
                                debounce=True,
                                persistence=True,
                                persistence_type='session',
                                required=False,
                                style={"fontSize": "0.750rem"}
                            ),
                        ],
                                width=3
                                ),
                        dbc.Col([
                            dbc.Label("Results",
                                      html_for="n_results",
                                      className='description_h4'),
                            dbc.Input(
                                type="number",
                                id="n_results",
                                min=1,
                                max=99,
                                value=5,
                                debounce=True,
                                style={"fontSize": "0.750rem"}
                            )
                        ],
                                width=2
                                ),
                        dbc.Col([
                            dbc.Button(
                                id='search_biblio',
                                outline=True,
                                color="primary",
                                n_clicks=0,
                                children='Search',
                                disabled=False,
                                className="description_h4 position-absolute bottom-0 start-0",
                            ),
                        ], className="align-self-end position-relative",
                                width=2
                                ),
                        dbc.Col(html.Div(id="animate_search"),
                                className="d-flex align-self-end")
                    ], className="d-flex g-3"),

                    # Biblio output
                    html.Div(
                        id='biblio_search_output',
                        children=[
                            html.Table(
                                id='biblio_table',
                                children=[
                                    dash_table.DataTable(
                                        id='biblio_datatable',
                                        columns=[
                                            {'name': 'Author(s)',
                                             'id': 'Author(s)'},
                                            {'name': 'Title',
                                             'id': 'Title'},
                                            {'name': 'Year',
                                             'id': 'Year'},
                                            {'name': 'Journal',
                                             'id': 'Journal'},
                                            {'name': 'Transcript',
                                             'id': 'Transcript'},
                                            {'name': 'Doi',
                                             'id': 'Doi',
                                             'presentation': 'markdown'},
                                            {'name': 'Abstract',
                                             'id': 'Abstract'}],
                                        hidden_columns=['Abstract'],
                                        markdown_options={"html": True},
                                        active_cell=initial_active_cell_bib,
                                        selected_cells=[],
                                        page_current=0,
                                        style_data={
                                            'whiteSpace': 'normal',
                                            'height': 'auto',
                                        },
                                        style_cell={
                                            'font-size': '0.625rem',
                                            'font-family': 'sans-serif',
                                            'padding': '0.5em 1em'
                                        },
                                        css=[{"selector": ".show-hide",
                                              "rule": "display: none"},
                                             {'selector':
                                              '.previous-next-container',
                                              'rule': 'font-size: 0.625rem;'}],
                                        page_size=5,
                                        row_selectable='multi',
                                        selected_rows=[],
                                                )  # Datatable
                                ], className="mt-3"),  # Table
                            dbc.Button(
                                "Search log",
                                id="btn_biblio_log",
                                outline=True,
                                className="w-auto description_h4 mb-3",
                                color="primary",
                                size="sm",
                                n_clicks=0,
                            ),
                            dbc.Collapse(
                                dbc.Card(
                                    dbc.CardBody(id='biblio_log'),
                                    className='mb-3'
                                ),
                                id="biblio_card",
                                is_open=False,
                            ),
                            dbc.Row([
                                dbc.Col([
                                    dbc.Button(
                                        id='btn_save_biblio',
                                        outline=True,
                                        color="secondary",
                                        children="Download bibliography",
                                        className="description_h4",
                                        size="sm"
                                    )
                                ]),
                                dbc.Col([
                                    dcc.Download(
                                        id="download_biblio"
                                    )
                                ]),
                            ])
                        ], style={'display': 'none'}),
                    dcc.Store(
                        data=[],
                        id='bib_records'
                    ),
                    dcc.Store(
                        data=True,
                        id='active_searching'
                    ),
                    dcc.Store(
                        data=[],
                        id='notes_input'
                    ),
                    html.Div([
                        html.Div(id='entry_title',
                                 className='bibOutputTitle'),
                        html.Div(id='entry_author',
                                 className='description_h4'),
                        html.Div(id='entry_date',
                                 className='bibOutputDate'),
                        html.Div(id='entry_abstract',
                                 className='bibOutputAbstract')
                    ])]),
        ], color='secondary')
    ]),  # Section
    # References
    html.Section([
        dbc.Alert([
            html.Div([
                html.H2('References',
                        className='flex-fill'),
                html.I(className=icon_hide,
                       id="references_btn",
                       n_clicks=0)
            ], id="references_header",
                     className=hdr_div),
            dbc.Collapse(
                is_open=True,
                id="references_description_html",
                className="mt-4",
                children=[
                    html.Div(
                        bib.get_description('references',
                                            'test_cases'),
                        className='description_h3',
                        id='study_description'
                    )
                ]
            )
        ], color='secondary')
    ])
])

# ------------------------------
# [F] Call back functions
clientside_callback(
    """(w) => {
    return `${w}`
    }""",
    Output("viewport-container", "children"),
    Input("breakpoints", "width"),
)


# START Callbacks for header display
@callback(
    [Output("description_html", "is_open"),
     Output("description_btn", "className"),
     Output("description_header", "className")],
    [Input("description_btn", "n_clicks")],
    [State("description_html", "is_open")],
    prevent_initial_call=True
)
def toggle_collapse_description(n, is_open):
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


@callback(
    [Output("settings_description_html", "is_open"),
     Output("settings_btn", "className"),
     Output("settings_header", "className")],
    [Input("settings_btn", "n_clicks")],
    [State("settings_description_html", "is_open")],
    prevent_initial_call=True
)
def toggle_collapse_settings(n, is_open):
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


@callback(
    [Output("summary_description_html", "is_open"),
     Output("summary_btn", "className"),
     Output("summary_header", "className")],
    [Input("summary_btn", "n_clicks")],
    [State("summary_description_html", "is_open")],
    prevent_initial_call=True
)
def toggle_collapse_summary(n, is_open):
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


@callback(
    [Output("candidates_description_html", "is_open"),
     Output("candidates_btn", "className"),
     Output("candidates_header", "className")],
    [Input("candidates_btn", "n_clicks")],
    [State("candidates_description_html", "is_open")],
    prevent_initial_call=True
)
def toggle_collapse_candidates(n, is_open):
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


@callback(
    [Output("search_refs_description_html", "is_open"),
     Output("search_refs_btn", "className"),
     Output("search_refs_header", "className")],
    [Input("search_refs_btn", "n_clicks")],
    [State("search_refs_description_html", "is_open")],
    prevent_initial_call=True
)
def toggle_collapse_references(n, is_open):
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
# END callback for header display


@callback(
    [Output("zhang_btn", "children"),
     Output("zhang_btn", "color"),
     Output("zhang_card", "className"),
     Output("zhang_card", "color"),
     Output("oliver_btn", "children"),
     Output("oliver_btn", "color"),
     Output("oliver_card", "className"),
     Output("oliver_card", "color"),
     Output("dataSet_name", "data")],
    [Input("oliver_btn", "n_clicks"),
     Input("zhang_btn", "n_clicks")],
    [State("oliver_btn", "children"),
     State("zhang_btn", "children")],
    prevent_initial_call=True
)
def toggle_dataset(o_clicks, z_clicks, o_state, z_state):
    if int(o_clicks) > 0 or int(z_clicks) > 0:
        if o_state == "Select":
            btn_state_on = 'Select'
            btn_color_on = 'primary'
            card_class_on = "bg-secondary bg-opacity-10 text-secondary"
            card_color_on = "secondary"
            btn_state_off = 'Selected'
            btn_color_off = 'info'
            card_class_off = "bg-primary bg-opacity-75"
            card_color_off = "primary"
            data_set = 'Oliver-2022'
        elif z_state == "Select":
            btn_state_off = 'Select'
            btn_color_off = 'primary'
            card_class_off = "bg-secondary bg-opacity-10 text-secondary"
            card_color_off = "secondary"
            btn_state_on = 'Selected'
            btn_color_on = 'secondary'
            card_class_on = "bg-primary bg-opacity-75"
            card_color_on = "primary"
            data_set = 'Zhang-2021'

        return btn_state_on, btn_color_on, card_class_on, card_color_on, \
            btn_state_off, btn_color_off, card_class_off, card_color_off, \
            data_set


@callback(
    Output("ivars", "data"),
    Input("dataSet_name", "data")
)
def set_variables(name):
    ivars = bib.import_vars(name)
    return ivars


@callback(
    [Output("dataSet_title", "children"),
     Output("dataSet_authors", "children"),
     Output("dataSet_date", "children")],
    Input("ivars", "data")
)
def dataSet_info(ivars):
    dataSet_title = bib.get_description('title',
                                        ivars['ibase'], False)
    dataSet_authors = bib.get_description('authors',
                                        ivars['ibase'], False)
    dataSet_date = bib.get_description('date',
                                        ivars['ibase'], False)
    return dataSet_title, dataSet_authors, dataSet_date


@callback(
    [Output("pydeg_data", "data"),
     Output("miRNA_data", "data")],
    [Input("dataSet_name", "data"),
     Input("ivars", "data"),
     Input("pydeg_settings_item", "value")]
)
def import_data(name, ivars, pysettings):
    data_dict = bib.import_data(name, ivars, pysettings)
    pydeg_data = data_dict["pydeg_df"]
    miRNA_data = data_dict["miRNA_df"]
    return pydeg_data.to_dict('records'), miRNA_data.to_dict('records')


@callback(
    [Output("class1_value", "options"),
     Output("comparison_value", "options")],
    [Input('x_axis_value', 'value'),
     Input('group_value', 'value'),
     Input("pydeg_data", "data"),
     Input('ivars', 'data')]
)
def dropdown_barplot(x1, x2, pydeg_data, ivars):
    pydeg_df = pd.DataFrame.from_records(pydeg_data)
    comparisons = [x for x in
                   sorted(pydeg_df.comparison.unique())]
    class1_values = [x for x in
                     sorted(pydeg_df.category_1.unique())]
    if x1 == 'comparison':
        options_2 = [{'label': icomp,
                      'value': icomp,
                      'disabled': True} for icomp in comparisons]
    else:
        options_2 = [{'label': icomp,
                      'value': icomp} for icomp in comparisons]

    if x2 == 'category_1':
        options_1 = [{'label': icat1,
                      'value': icat1,
                      'disabled': True} for icat1 in class1_values]
    else:
        options_1 = [{'label': icat1,
                      'value': icat1} for icat1 in class1_values]

    return options_1, options_2


@callback(
    Output('peak_count_barplot', 'figure'),
    [Input('pydeg_data', 'data'),
     Input('ivars', 'data'),
     Input('x_axis_value', 'value'),
     Input('group_value', 'value'),
     Input("class1_value", "value"),
     Input("comparison_value", "value"),
     Input(ThemeChangerAIO.ids.radio("theme"), "value"),
     Input('viewport-container', 'children')
     ]
)
def peak_count_barplot(pydeg_data, ivars,
                       x1, x2, x3, x4, theme, viewport):
    wd = int(viewport)
    pydeg_df = pd.DataFrame.from_records(pydeg_data)
    if x3:
        pydeg_df_filtered = pydeg_df[(~pydeg_df['category_1'].isin(x3)) &
                                     (pydeg_df['comparison'] != x4)]
    else:
        pydeg_df_filtered = pydeg_df[pydeg_df['comparison'] != x4]

    if x2:
        counts = pydeg_df_filtered.groupby([x1, x2]).size()
        df_counts = counts.reset_index().rename(
            columns={'level_0': x1, 'level_1': x2, 0: 'Peak number'}
        )
        if df_counts.columns[1] == 'category_1':
            df_counts.loc[:, 'category_1'] = df_counts['category_1'].\
                astype(str).replace(ivars['cat1_dict'])
            df_counts = df_counts.sort_values(by='category_1',
                                              ascending=True)
        elif df_counts.columns[1] == 'category_2':
            df_counts.loc[:, 'category_2'] = df_counts['category_2'].\
                    replace(ivars['cat2_dict'])

        df_counts = df_counts.rename(columns=bib.new_columns)
        x1_plot = bib.new_columns[x1]
        x2_plot = bib.new_columns[x2]
        fig = px.bar(df_counts, x=x1_plot, y='Peak number',
                     color=x2_plot, barmode='group',
                     template=template_from_url(theme),
                     )
    else:
        counts = pydeg_df_filtered.groupby([x1]).size()
        df_counts = counts.reset_index().rename(
            columns={'level_0': x1, 0: 'Peak number'}
        )
        df_counts = df_counts.rename(columns=bib.new_columns)
        x1_plot = bib.new_columns[x1]
        fig = px.bar(df_counts, x=x1_plot, y='Peak number',
                     template=template_from_url(theme),
                     )

    axis_title = bib.calcFontSize(wd)
    x_tick = axis_title-axis_title*0.14
    legend_font = axis_title-axis_title*0.28
    y_tick = axis_title-axis_title*0.21
    fig.update_layout(
        legend=dict(
            orientation="h",
            entrywidth=100,
            yanchor="bottom",
            y=-0.6,
            xanchor="right",
            x=1,
            font=dict(size=legend_font)
        )
    )

    fig.update_xaxes(
        title_font=dict(size=axis_title),
        tickfont=dict(size=x_tick)
    )

    fig.update_yaxes(
        title_font=dict(size=axis_title),
        tickfont=dict(size=y_tick)
    )

    return fig


@callback(
    Output('factor_description', 'children'),
    Input('x_axis_value', 'value'),
    Input('group_value', 'value')
)
def display_group_info(group1, group2):
    """
    Based on the selected option from two drop-down menus,
    search on a tsv file 'pydeg_group_description' for
    the associated information.
    """
    group1_description = bib.pydeg_group_description[
        bib.pydeg_group_description['group'].eq(group1)]
    group2_description = bib.pydeg_group_description[
        bib.pydeg_group_description['group'].eq(group2)]

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


@callback(
    Output("dropdown_pytable", "children"),
    [Input('pydeg_data', 'data'),
     Input("pydeg_settings_item", "value")]
)
def dropdown_pytable(pydeg_data, pysettings):
    pydeg_df = pd.DataFrame.from_records(pydeg_data)
    dropdown_row = dbc.Row([
        dbc.Col([
            dcc.Dropdown(
                id="class1_drop",
                placeholder='Classification 1',
                options=[x for x in
                         sorted(pydeg_df.category_1.unique())],
                className="dropdownFont"
            )], className='five.columns'),
        dbc.Col([
            dcc.Dropdown(
                id="class2_drop",
                placeholder='Classification 2',
                options=[x for x in
                         sorted(pydeg_df.category_2.unique())],
                className="dropdownFont"
            )], className='five.columns'),
        dbc.Col([
            dcc.Dropdown(
                id="feat_drop",
                placeholder='Feature',
                options=[x for x in
                         sorted(pydeg_df.feature_type.unique())],
                className="dropdownFont"
            )], className='five.columns'),
        dbc.Col([
            dcc.Dropdown(
                id="plot_drop",
                placeholder='Has peak plot',
                options=[x for x in
                         sorted(pydeg_df.plot_link.unique())],
                className="dropdownFont"
            )], className='five.columns'),
        dbc.Col([
            dcc.Dropdown(
                id="mirna_drop",
                placeholder='Has miRNA alignment',
                options=[x for x in
                         sorted(pydeg_df.miRNA_link.unique())],
                className="dropdownFont"
            )], className='five.columns')
    ])

    return dropdown_row


@callback(
    Output('py_table', 'data'),
    [Input('pydeg_data', 'data'),
     Input("class1_drop", 'value'),
     Input("class2_drop", 'value'),
     Input("feat_drop", 'value'),
     Input("plot_drop", 'value'),
     Input("mirna_drop", 'value')],
    prevent_initial_call=True
)
def update_dropdown_options(pydeg_data, class_1,
                            class_2, feat, plot, mirna):
    pydeg_df = pd.DataFrame.from_records(pydeg_data)

    if class_1:
        pydeg_df = pydeg_df[pydeg_df.category_1 == class_1]
    if class_2:
        pydeg_df = pydeg_df[pydeg_df.category_2 == class_2]
    if feat:
        pydeg_df = pydeg_df[pydeg_df.feature_type == feat]
    if plot:
        pydeg_df = pydeg_df[pydeg_df.plot_link == plot]
    if mirna:
        pydeg_df = pydeg_df[pydeg_df.miRNA_link == mirna]

    return pydeg_df.to_dict('records')


@callback(
    Output("description_html", "children"),
    Input("ivars", "data")
)
def dataSet_description(ivars):
    dataSet_description = html.Div(
        bib.get_description('test_case_description',
                            ivars['ibase']),
        className='description_h3',
        id='study_description'
    )
    return dataSet_description


@callback(
     Output('py_table', 'active_cell'),
     [Input('py_table', 'page_current'),
      Input('py_table', 'page_size')],
     prevent_initial_call=True
 )
def reset_active_cell(page, size):
    row_id = page*size
    active_cell = {'row': 0, 'column': 3, 'row_id': row_id}
    return active_cell


@callback(
    [Output('tabs-content-decay-plots', 'children'),
     Output('miRNAplot_df', 'data')],
    [Input('py_table', 'active_cell'),
     Input('pydeg_data', 'data'),
     Input('miRNA_data', 'data'),
     Input('tab_plots', 'value')],
)
def render_tab_content(active_cell, pydeg_data, miRNA_data, tab):
    if active_cell is None:
        return no_update

    pydeg_df = pd.DataFrame.from_records(pydeg_data)
    row = active_cell['row_id']
    transcript = pydeg_df.at[row, 'tx_name']
    comparison = pydeg_df.at[row, 'comparison']
    cat1 = str(pydeg_df.at[row, 'category_1'])
    cat2 = str(pydeg_df.at[row, 'category_2'])
    category = cat1 + "-" + cat2
    header = f'Decay plot for {transcript} ({category}) \
from comparison {comparison}'

    df = pd.DataFrame.from_records(miRNA_data)
    dff = df.loc[df['Transcript'] == transcript]
    if tab == 'gene_plot':
        gene_plot_link = pydeg_df.at[row, 'gene_plot_link']
        if isinstance(gene_plot_link, str):
            decay_plot = html.Div([
                html.Div(id='miRNA_plot'),
                html.Div(children=header,
                         className="header_dPlot"),
                html.Img(src=gene_plot_link,
                         className='centered_image',
                         **{'data-table': transcript})
            ])
        else:
            decay_plot = html.P(
                f'No plot was produced for transcript {transcript}',
                className='no_output_w')
    elif tab == 'peak_plot':
        peak_plot_link = pydeg_df.at[row, 'peak_plot_link']
        if isinstance(peak_plot_link, str):
            decay_plot = html.Div([
                html.Div(id='miRNA_plot'),
                html.Div(children=header,
                         className="header_dPlot"),
                html.Img(src=peak_plot_link,
                         className='centered_image',
                         **{'data-table': transcript})
            ])
        else:
            decay_plot = html.P(
                f'No plot was produced for transcript {transcript}',
                className='no_output_w')
    elif tab == 'miRNA_tab':
        page_size = len(dff)
        if page_size > 0:
            row_id = dff.loc[dff.index[0], 'id']
            decay_plot = html.Div([
                dash_table.DataTable(
                    id="miRNA_datatable",
                    columns=[
                        {'name': ['', 'Transcript'],
                         'id': 'Transcript'},
                        {'name': ['', 'miRNA'],
                         'id': 'miRNA'},
                        {'name': ['Alignment', 'peak region'],
                         'id': 'Score_x'},
                        {'name': ['Alignment', 'Peak'],
                         'id': 'Score_y'},
                        {'name': ['', 'Comparison'],
                         'id': 'Comparison'},
                    ],
                    data=dff.to_dict("records"),
                    merge_duplicate_headers=True,
                    style_cell_conditional=[
                        {"if": {"column_id": "Transcript"}, "width": "20%"},
                        {"if": {"column_id": "miRNA"}, "width": "20%"},
                        {"if": {"column_id": "Score_x"}, "width": "9%"},
                        {"if": {"column_id": "Score_y"}, "width": "15%"},
                        {"if": {"column_id": "Comparison"}, "width": "30%"},
                    ],
                    style_cell={"fontSize": "0.625rem",
                                'font-family': 'sans-serif',
                                'padding': '0.5em 0.7em'},
                    active_cell={'row': 0, 'column': 0, 'row_id': row_id},
                    # selected_cells=[],
                    page_size=3
                ),
                html.Div(id='miRNA_plot')
            ], className="mt-3")
        else:
            # No candidates found
            decay_plot = html.Div(
                html.P(f'No alignment was obtained for\
the peak of transcript {transcript}',
                       className='no_output_w')
            )

    return decay_plot, dff.to_dict('records')


@callback(
    Output('miRNA_plot', 'children'),
    Input('miRNAplot_df', 'data'),
    Input('tab_plots', 'value'),
    prevent_initial_call=True
)
def render_miRNA_plot(miRNA_data, tab):
    if tab == 'miRNA_tab':
        miRNA_df = pd.DataFrame.from_records(miRNA_data)
        mirmap_link = miRNA_df.loc[0, 'mirmap_link']
        globalAln_link = miRNA_df.loc[0, 'global_link']
        miRNA_alignment_plot = bib.draw_miRNAplot(mirmap_link, globalAln_link)
        return miRNA_alignment_plot
    else:
        return None


@callback(
    Output('miRNA_plot', 'children', allow_duplicate=True),
    Input("miRNA_datatable", "active_cell"),
    State('miRNAplot_df', 'data'),
    prevent_initial_call=True
)
def update_miRNA_plot(active_cell, miRNA_data):
    miRNA_df = pd.DataFrame.from_records(miRNA_data)
    mirmap_link = miRNA_df.loc[active_cell['row'], 'mirmap_link']
    globalAln_link = miRNA_df.loc[active_cell['row'], 'global_link']
    miRNA_alignment_plot = bib.draw_miRNAplot(mirmap_link, globalAln_link)

    return miRNA_alignment_plot


@callback(
    Output('animate_search', 'children', allow_duplicate=True),
    [Input('search_biblio', 'n_clicks')
     ],
    prevent_initial_call=True
)
def animate_search(clicks):
    if clicks:
        search_status = dbc.Spinner(color="primary", type="grow")
    return search_status


@callback(
    [Output('biblio_datatable', 'data'),
     Output('bib_records', 'data'),
     Output('biblio_search_output', 'style'),
     Output('biblio_log', 'children'),
     Output('animate_search', 'children')],
    [Input('search_biblio', 'n_clicks')],
    [State('py_table', 'selected_rows'),
     State('py_table', 'data'),
     State('op_term', 'value'),
     State('ncbi_email', 'value'),
     State('n_results', 'value')],
    prevent_initial_call=True
)
def result_biblio(search_biblio, selected_tx, pydeg_data,
                  op_term, email, n_results):
    if search_biblio > 0:
        loaded_df = pd.DataFrame.from_records(pydeg_data)
        itx_list = loaded_df.loc[selected_tx, 'tx_name']
        bib.setEmail(email)

        entries_list = []
        biblio_list = []
        entries_found = 0
        not_found = []
        for transcript in itx_list:
            print('------------------------')
            print(f'Working on {transcript}')
            transcript_2 = re.sub(r'\.[0-9]$', '', transcript)

            if op_term:
                search_term = f"{transcript} OR {transcript_2} AND {op_term}"
            else:
                search_term = f"{transcript} OR {transcript_2}"

            print("Searching")
            ref_ids = bib.esearch('pmc', search_term, n_items=n_results)
            if ref_ids:
                ref_ids_pubmed = bib.getPubmedId(ref_ids)

                print("Parsing records")
                entries = bib.getRefRecords('pubmed', ref_ids_pubmed)
                entries_list.append(entries)
                n_entries = len(entries)
                entries_found += n_entries
                biblio_df = bib.getBibDF(entries, transcript)
                biblio_list.append(biblio_df)
                print("Finished")
            else:
                not_found.append(transcript)

        entries_combined = [item for sublist in
                            entries_list for item in sublist]
        biblio_df_combined = pd.concat(biblio_list, axis=0)
        display_status = {'display': 'block'}
        # Log results
        if not_found:
            txs = ', '.join(not_found)
            biblio_log = html.Div([
                html.P(f"Entries found: {entries_found}",
                       className="text-success"),
                html.P(f"No results for the following transcripts: \
    {txs}", className="text-warning")
            ], className='description_h4')
        else:
            biblio_log = html.Div([
                html.P(f"Entries found: {entries_found}",
                       className="text-success")
            ], className='description_h4')
        return biblio_df_combined.to_dict('records'), \
            entries_combined, display_status, biblio_log, None


@callback(
    [Output('entry_title', 'children'),
     Output('entry_author', 'children'),
     Output('entry_date', 'children'),
     Output('entry_abstract', 'children'),],
    [Input('biblio_datatable', 'active_cell'),
     Input('bib_records', 'data'),
     Input('biblio_datatable', 'page_current'),
     Input('biblio_datatable', 'page_size')
     ],
    prevent_initial_call=True
)
def cell_clicked_bib(active_cell, entries, page_current, page_size):
    if active_cell:
        biblio_df = entries
        row = active_cell['row'] + (page_current*page_size)
        key_list = ["title", "authors", "journal", "publication_date",
                    "doi", "abstract", "keywords"]
        ititle = bib.printBibSection(biblio_df[row], key_list[0:1], "\n")
        iauthor = bib.printBibSection(biblio_df[row], key_list[1:2], "\n")
        idate = bib.printBibSection(biblio_df[row], key_list[2:5], "; ")
        iabstract = bib.printBibSection(biblio_df[row], key_list[5:], "\n")
    else:
        ititle = iauthor = idate = iabstract = None
    return ititle, iauthor, idate, iabstract


@callback(
    Output("download_notes", "data"),
    [Input("btn_save_notes", "n_clicks")],
    [State('biblio_notes', 'value')],
    prevent_initial_call=True,
)
def save_notes(n_clicks, value):
    return dict(content=value, filename="notes.txt")


@callback(
    Output("download_biblio", "data"),
    [Input("btn_save_biblio", "n_clicks")],
    [State('biblio_datatable', 'selected_rows'),
     State('bib_records', 'data')],
    prevent_initial_call=True
)
def save_refs(n_clicks, selected_refs, bib_records):
    pubmed_ids = [bib_records[i]['pubmed_id'] for i in selected_refs]
    biblio_data = bib.fetchRefs('pubmed', pubmed_ids)
    return dict(content=biblio_data, filename="References.medline")


@callback(
    Output("biblio_card", "is_open"),
    Input("btn_biblio_log", "n_clicks"),
    [State("biblio_card", "is_open")],
)
def toggle_biblio_log(n, is_open):
    if n:
        return not is_open
    return is_open


@callback(
    [Output("offcanvas_notes", "is_open"),
     Output('offcanvas_notes', 'children')],
    Input("btn_notes", "n_clicks"),
    [State('notes_input', 'children'),
     State("offcanvas_notes", "is_open")],
)
def toggle_notes(n1, text, is_open):
    if n1:
        notes_div = html.Div([
            dcc.Textarea(
                id='biblio_notes',
                placeholder='Write notes here...',
                style={'width': '100%', 'fontSize': '0.875rem'},
                value=text
            ),
            dbc.Button(
                id="btn_save_notes",
                children="Save notes",
                outline=True,
                color="secondary",
                size="sm"
            ),
            dcc.Download(
                id="download_notes"
            )
        ])
        return not is_open, notes_div
    return is_open, None


@callback(
    Output('notes_input', 'children'),
    Input('biblio_notes', 'value'),
    prevent_initial_call=True
)
def update_notes(value):
    return value
