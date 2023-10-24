import dash
from dash import html, Input, Output, callback, State
# import plotly.express as px
import dash_bootstrap_components as dbc
# import json
from dash_dangerously_set_inner_html import DangerouslySetInnerHTML
from . import bibsearch as bib

dash.register_page(__name__, path='/', name='Description')

# Class definition for headers
icon_hide = "fas fa-angle-up me-4 align-self-center icon"
icon_show = "fas fa-angle-down me-4 align-self-center icon"
hdr_hide = 'd-flex bg-opacity-25 bg-primary'
hdr_show = 'd-flex bg-opacity-10 bg-primary'
hdr_div = 'd-flex bg-opacity-25 bg-primary'
layout = html.Main([
    # Summary
    html.Section([
        dbc.Alert([html.H2('Summary',
                           className='border-bottom border-3 border-secondary pb-2 mb-4'),

                   DangerouslySetInnerHTML('''
<p>
This Dash app summarizes the analysis of mRNA degradation fragments using a third party python script (<code>PyDegradome</code> <a href="#citeproc_bib_item_2">Gaglia, Rycroft, and Glaunsinger 2015</a>) followed by a custom classification. The result is a table of transcripts (containing signals of degradation fragments, <i>peaks</i>) sorted according to two main classification criteria. Plots for the accumulation of degradation fragments along the transcript are linked for the user to inspect. Additionally, alignments of sequences around peaks and known miRNA are shown for alignment scores above certain threshold. Finally, an option to search for related literature based on the transcript name and other search criteria is presented. As a proof of concept the analysis of two data sets is presented.
</p>
        ''')], color='primary')
    ]),
    # Description Pydegradome
    html.Section([
        dbc.Alert([
            html.Div([
                html.H2('PyDegradome description',
                        className='flex-fill'),
                html.I(className=icon_hide,
                       id="pydeg_btn",
                       n_clicks=0)
            ], id="pydeg_header",
                     className=hdr_div),
            dbc.Collapse(
                is_open=True,
                id="pydeg_description_html",
                className="mt-4",
                children=DangerouslySetInnerHTML('''
<p>
<code>PyDegradome</code> is a python script developed by <a href="#citeproc_bib_item_2">Gaglia, Rycroft, and Glaunsinger 2015</a> for the study of mRNA targets of a viral endonuclease upon infection of human cells. Broadly, the script compares two samples, one having the factor thought to be involved in the production of mRNA degradation fragments (<i>test</i>) and another devoided of such factor (<i>control</i>); from this comparison a table is presented with genomic coordinates for regions where degradation intermediates are accumulated differently in the test sample. Differential accumulation of degradation intermediates in the control sample is not reported unless the test is repeated and the order of samples is reversed.
</p>

<p>
Statistical significance for the identified regions (referred as <i>peaks</i>)  is determined upon comparison with a threshold value calculated using the ratio of reads along the exon (<i>test:control</i>), multiplied by a user-defined <b><i>multiplication factor</i></b> for the production of degradation fragments. <b><i>Confidence levels</i></b> for the identification of significant differences are also defined by the user. 
</p>

<p>
Figure <a href="#org6aa7b74">1</a> illustrates how peaks are identified, briefly, reads for the test and control samples are mapped to a single exon (panel 1), collapsed to the 5' end and counted (panel 2), with the total count and the multiplicative factor, the threshold at the exon level is calculated (panel 3), finally, the significance of the different regions is be considered based on the selected confidence level (panel 4). Bars on panels 3 and 4 represent reads grouped over a 4-nt window (user-defined); neighboring bars that clear the significance threshold are grouped into a single <i>peak</i>.
</p>

<figure id="org6aa7b74">
<img class="img_full" src="./assets/Images/pydeg-method.png" alt="Illustration of PyDegradome analysis" class="image" />
<figcaption class="caption_full"><span>Figure 1: </span>Illustration of how peaks are identified by the <code>pyDegradome</code> script. Based on figure 1 <a href="#citeproc_bib_item_2">(Gaglia, Rycroft, and Glaunsinger 2015)</a></figcaption>
</figure>

<p>
Significant peaks are reported in a text file (see below); genomic coordinates for the peak boundaries (start and end) and counts for the reads present in it are also indicated however, note that no information about the associated transcript is present. 
</p>

<div class="description_h4">
<label class="org-src-name"><span class="listing-number">Listing 1: </span>Sample of <code>PyDegradome</code> output showing the information for 2 regions (<i>peaks</i>) with significant differences. For each peak 3 lines are reported. The tab-separated line maps to the information provided in the header (line starting with '#'). The next 2 are the read counts found 10nt upstream and downstream of the highest point of the peak.</label><pre class="src src-nil" id="org6b27960">
#chr	strand	peak_start	peak_stop	max_peak_test	ratio	factor	max_count_test	max_count_ctl	tot_count_test	tot_count_ctl
5	+	10015540	10015543	10015540	4.79591836735	14.387755102	3	0	9	0
[1, 0, 0, 1, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 1, 0]
[2, 1, 1, 3, 0, 2, 0, 2, 0, 0, 3, 3, 3, 0, 0, 0, 0, 0, 0, 1, 0, 5, 0]
5	+	10136622	10136623	10136622	0.546218487395	1.63865546218	3	0	3	0
[0, 1, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0]
[0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 3, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0]
</pre>
</div>        
                     '''))
        ], color='secondary')
    ]),
    # Post Pydegradome
    html.Section([
        dbc.Alert([
            html.Div([
                html.H2('Classification post PyDegradome',
                        className='flex-fill'),
                html.I(className=icon_hide,
                       id="classification_btn",
                       n_clicks=0)
            ], id="classification_header",
                     className=hdr_div),
            dbc.Collapse(
                is_open=True,
                id="classification_description_html",
                className="mt-4",
                children=DangerouslySetInnerHTML("""
<p>
As shown above the results of the <code>PyDegradome</code> script require some level of post processing by the user; at least, the annotation of the genomic coordinates to the corresponding gene and (following the original work,  <a href="#citeproc_bib_item_2">Gaglia, Rycroft, and Glaunsinger 2015</a> Figure 1D) the identification of shared peaks between replicates. Additionally, there are not recommended settings for the analysis rather, the authors advice to test different values for the <i>multiplicative factor (MF)</i> and the <i>confidence interval (CI)</i>) in order to identify those that would result in peaks with the highest confidence  (See <a href="#citeproc_bib_item_2">Gaglia, Rycroft, and Glaunsinger 2015</a> p. 5 and Figure S2B).
</p>

<p>
For the present classification, the confidence of the identified peaks was assessed by inspection of the accumulation of degradation intermediates along the whole transcript. Figure <a href="#org0b7f53b">2</a> shows a plot of all the reads mapped to a given transcript. Such plot is referred as a decay plot but also as <i>t-plots</i> (<a href="#citeproc_bib_item_3">German et al. 2008</a>) or <i>D-plots</i> (<a href="#citeproc_bib_item_6">Nagarajan et al. 2019</a>). Reads found on the test samples are plotted on the top panel (with replicates identified by different shades of red) whereas, on the bottom panel are those reads corresponding to the control sample (with replicates identified by different shades of blue). A vertical yellow line highlights the identified peak and the corresponding gene model is illustrated at the bottom.  
</p>


<figure id="org0b7f53b">
<img src="./assets/Images/05_Gene_1-C_AT5G17070_1_0_95_4_2.jpg" alt="05_Gene_1-C_AT5G17070_1_0_95_4_2.jpg" class="image" />
<figcaption><span class="figure-number">Figure 2: </span>Example of a degradome plot (D-plot). Reads corresponding to test samples are plotted on the top whereas, those for control samples are on the bottom. A red arrowhead marks the region identified by <code>PyDegradome</code> as significantly different between samples (<i>i.e.</i> <i>peak</i>). Note that eventhough there is a clear difference between samples, this is not the highest signal along the transcript. Values were scaled using <code>DESeq</code> scale factor.</figcaption>
</figure>

<p>
While conducting these observations it was noted that in a number of cases the identified peak do not correspond to the highest read along the gene region and it is even possible to find reads corresponding to the control sample with larger values than those of the test sample. Such cases arise due to the 'exon-focused' way by which <code>PyDegradome</code> identifies peaks, and it is not evident unless all gene reads are plotted. Figure <a href="#orgd79daf9">3</a> illustrates the above point, panel A aims to depict how the script identifies a region on a given exon as significantly different without considering the information on neighboring exons whereas, panel B shows a complete picture where it is evident that reads corresponding to the control sample (on exon 2) have the largest value along the gene. The horizontal black line represents the significance threshold that is calculated per-exon. 
</p>


<figure id="orgd79daf9">
<img src="./assets/Images/pydeg-limitations.png" alt="pydeg-limitations.png" class="image" />
<figcaption><span class="figure-number">Figure 3: </span>Limitations of <code>PyDegradome</code> script in identifing the region with the highest read along the gene-region. This is due to the exon-focused nature on which peaks are identied (A) and how this may ignore the information of neighboring exons where reads with larger values could be found (B)</figcaption>
</figure>


<p>
To account for the above observations a filtering and classification step that looks at the whole transcript was included. Peaks identified in the <code>PyDegradome</code> output were classified into 4(+1) categories based on the ratio between the highest read value in the peak region and the highest read found outside the peak region (<i>peak:background</i> Figure <a href="#org84d25ac">4</a>). For the first 3 categories an additional requirement was that peaks should be present in both replicates (<i>i.e.</i> their coordinates should overlap by at least 1nt) whereas, the last category considered peaks found in only one replicate. Peak coordinates were used to calculate the number of reads in the remaining replicate.
</p>

<ul class="org-ul">
<li><b>Category 1:</b> Peaks with a peak:background ratio higher than 1</li>
<li><b>Category 2:</b> Peaks with a ratio between 0.8 and 1</li>
<li><b>Category 3:</b> Peaks with a ratio between 0.7-0.8.</li>
<li><b>Category 4:</b> Peaks with a ratio higher than 1 in one replicate and higher than 0.8 in the other were selected.</li>
<li><b>Category 0:</b> Peaks that could not be classified into one of the above categories</li>
</ul>


<figure id="org84d25ac">
<img src="./assets/Images/pydeg-cat1-4.png" alt="pydeg-cat1-4.png" class="image" />
<figcaption><span class="figure-number">Figure 4: </span>Criteria for the classification/filtering of peaks obtained from the <code>PyDegradome</code> analysis. A ratio between the peak's reads (highlighted yellow) and those outside (background) was used to define 4 categories (Circled numbers). Categories 1 to 3 consider peaks found in both replicates while,  category \textcircled{4} analyzed peaks found only in one replicate. In the latter case, coordinates of the identified peak were applied on the other replicate as it were a significant 'peak'</figcaption>
</figure>

<p>
The above classification focused on reads found only in the test sample and thus, ignored the possibility of reads in the control sample with larger values (although the identification of the peak region partly considered this). To account for this possibility, a second classification step was included where the mean value of the peaks identified in test samples was divided by the mean of the highest read along the gene region in the control samples (Figure <a href="#org66d24be">5</a>). Based on this ratio, 3 categories were defined.
</p>
<ul class="org-ul">
<li><b>Category A:</b> Transcripts with a ratio higher than 1.</li>
<li><b>Category B:</b> Transcripts with two or more peaks and whose, ratio was higher than 50% of the control reads (median).</li>
<li><b>Category C:</b>, Transcripts that could not be classified in the above categories</li>
</ul>


<figure id="org66d24be">
<img src="./assets/Images/Pydeg-catA-C.png" alt="Pydeg-catA-C.png" class="image" />
<figcaption><span class="figure-number">Figure 5: </span>Criteria for the second classification of peaks obtained from the <code>PyDegradome</code> analysis. A ratio between the highest read in the peak region of test samples(mean of scaled values, red bars) and the highest read along the transcript from control samples  (light blue bar) were used to classify the peaks into 3 categories. Category A for peaks with a ratio greater than 1. Category B was designed for those transcripts with more than 1 peak, where the ratio was higher than the median value of control samples. Finally category C included all the candidates that could not be classified in the two categories above. No candidates were removed from further analysis.</figcaption>
</figure>

<p>
The two classification steps along with the annotation, summary and other steps were written in a series of <code>R</code> scripts (described <a href="https://ssl-blog.netlify.app/posts/degradome-analysis/degradome-code/">elsewhere</a>). The output is a number of tables summarizing the classification filtering steps. For the purpose of the Dash app, these were further formatted in python. The code is available on Github:
</p>
<ul class="org-ul">
<li><a href="https://github.com/ssl-bio/Degradome-analysis">Degradome analysis</a></li>
<li><a href="https://github.com/ssl-bio/R_postpydeg">R package</a></li>
<li><a href="https://github.com/ssl-bio/R_postpydeg">Dash app</a></li>
</ul>
                     """)),
        ], color='secondary')
    ]),
    # Alignment of miRNA species
    html.Section([
        dbc.Alert([
            html.Div([
                html.H2('Alignment of sequences around peaks with known miRNA',
                        className='flex-fill'),
                html.I(className=icon_hide,
                       id="alignment_btn",
                       n_clicks=0)
            ], id="alignment_header",
                     className=hdr_div),
            dbc.Collapse(
                is_open=True,
                id="alignment_description_html",
                className="mt-4",
                children=DangerouslySetInnerHTML("""
<p>
As stated above, the <code>PyDegradome</code> script (<a href="#citeproc_bib_item_2">Gaglia, Rycroft, and Glaunsinger 2015</a>) was developed to study the targets of an endonuclease; this, driven by the fact that most of the software available for the study of mRNA degradation fragments is aimed to the identification of potential miRNA targets. For plants that is particularly the case as the endonuclease involved in the ribosome-associated mRNA control remains to be discovered (<a href="#citeproc_bib_item_4">Karamyshev and Karamysheva 2018</a>). In this context, it was tested whether sequences found around the identified peaks could align to known miRNA species. To this end, two different approaches were tested. The first one, performed a <i>global pairwise-alignment</i> (implemented in <code>biopython</code> V. 1.75 (<a href="#citeproc_bib_item_1">Cock et al. 2009</a>) using a 22nt region centered around the peak whereas, the second one used a program that predicts the strength of the interaction between a miRNA and a putative target (<code>Mirmap</code>, <a href="#citeproc_bib_item_7">Vejnar and Zdobnov 2012</a>); in this case a longer sequence around peaks (40nt) was tested against to known miRNAs.
</p>

<p>
To set a threshold for the alignment score, which would account for the mismatches that are known to occur between miRNAs and their target sequences, an alignment file for miRNA species and their known (or predicted) targets was downloaded from <code>TarDB</code> (<a href="#citeproc_bib_item_5">Liu et al. 2021</a>) and used as input for both alignment programs. Alignment scores were obtained and their lowest values used as threshold for the selection of candidates. For the Arabidopsis data, the lowest <i>global pairwise-alignment</i> score found was 33 whereas, for the <code>Mirmap</code> alignment the lowest thermodynamic value was -13.8.
</p>

<p>
Note that due to storage limitation the alignment was conducted on sequences from peaks classified in one of the following categories, 1-A, 1-B, 2-A or 2-B
</p>      
                      """))
        ], color='secondary')
    ]),
    # Search for related literature
    html.Section([
        dbc.Alert([
            html.Div([
                html.H2('Search for related literature',
                        className='flex-fill'),
                html.I(className=icon_hide,
                       id="search_btn",
                       n_clicks=0)
            ], id="search_header",
                     className=hdr_div),
            dbc.Collapse(
                is_open=True,
                id="search_description_html",
                className="mt-4",
                children=DangerouslySetInnerHTML("""
<p>
A function to select literature mentioning transcripts of interest is included in the app. Each classified peak has a checkbox for selecting those transcripts with potential interest (Figure <a href="#orgfaf92ea">6</a>, upper panel). Then, on the section <i>Related literature for the selected transcripts</i> three fields for setting the query are presented (Figure <a href="#orgfaf92ea">6</a>, lower panel). The first is a registered e-mail address for carrying a query on NCBI (<a href="https://account.ncbi.nlm.nih.gov/">login link</a>); next is a field to include additional search terms (using an <code>AND</code> keyword in the background) and finally there is a field for specifying the number of bibliographic entries to be retrieved per selected transcript with a default of 5 and a maximum of 100 (an arbitrary value to avoid overloading the server).
</p>


<figure id="orgfaf92ea">
<img src="./assets/Images/search_references.png" alt="search_references.png" class="image" />
<figcaption><span class="figure-number">Figure 6: </span>Searching for related literature. Marked candidates (upper panel) are used to search for literature mentioning them (bottom panel). Search is carried on NCBI's Pub Med Central (PMC) for which a registered e-mail address is needed. Additional search terms can be included and the maximum number of retrieved references can be modified from a default of 5 up to 100.</figcaption>
</figure>

<p>
The search is carried using <code>biopython</code>'s <a href="https://biopython.org/docs/1.75/api/Bio.Entrez.html#"><code>Bio.entrez</code></a> querying NCBI's PubMed Central (PMC) database. Search results are displayed in a table with basic information of each retrieved reference and a section below with the reference abstract which can be accessed by clicking on a row. Besides the basic bibliographic information, the table also includes a column with the transcript's name used for the query, a link to the publication's web page and a checkbox for selecting and downloading references in <code>medline</code> (a format that is similar to <code>RIS</code> and can be parsed by most reference management software). 
</p>    
                     """))
        ], color='secondary')
    ]),
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
                                            'intro'),
                        className='description_h3',
                        id='study_description'
                    )
                ]
            )
        ], color='secondary')
    ])
])


# [F] Callback functions
@callback(
    [Output("pydeg_description_html", "is_open"),
     Output("pydeg_btn", "className"),
     Output("pydeg_header", "className")],
    [Input("pydeg_btn", "n_clicks")],
    [State("pydeg_description_html", "is_open")],
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
    [Output("alignment_description_html", "is_open"),
     Output("alignment_btn", "className"),
     Output("alignment_header", "className")],
    [Input("alignment_btn", "n_clicks")],
    [State("alignment_description_html", "is_open")],
    prevent_initial_call=True
)
def toggle_collapse_alignment(n, is_open):
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
    [Output("classification_description_html", "is_open"),
     Output("classification_btn", "className"),
     Output("classification_header", "className")],
    [Input("classification_btn", "n_clicks")],
    [State("classification_description_html", "is_open")],
    prevent_initial_call=True
)
def toggle_collapse_classification(n, is_open):
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
    [Output("search_description_html", "is_open"),
     Output("search_btn", "className"),
     Output("search_header", "className")],
    [Input("search_btn", "n_clicks")],
    [State("search_description_html", "is_open")],
    prevent_initial_call=True
)
def toggle_collapse_search(n, is_open):
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
    [Output("references_description_html", "is_open"),
     Output("references_btn", "className"),
     Output("references_header", "className")],
    [Input("references_btn", "n_clicks")],
    [State("references_description_html", "is_open")],
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
