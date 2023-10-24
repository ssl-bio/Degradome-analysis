import dash
from dash import html
import dash_bootstrap_components as dbc
from dash_bootstrap_templates import ThemeChangerAIO

doc_title = 'Dashboard for the analysis of mRNA degradation fragments'

e_stylesheets = [dbc.themes.JOURNAL, dbc.icons.FONT_AWESOME]
theme_change = ThemeChangerAIO(aio_id="theme",
                               radio_props={"className": "description_h4",
                                            'value': dbc.themes.JOURNAL},
                               button_props={"className": "rounded-pill me-4"})

app = dash.Dash(__name__, use_pages=True,
                assets_ignore='.#custom.css',
                meta_tags=[
                    {"name": "viewport",
                     "content": "width=device-width, initial-scale=1,\
                     maximum-scale=1"}
                ],
                external_stylesheets=e_stylesheets)
app.config.suppress_callback_exceptions = True

inavbar = dbc.Nav(
    [
        dbc.NavLink(
            [
                html.Div(page["name"]),
            ],
            href=page["path"],
            active="exact",
            className="border rounded-pill description_h3"
        )
        for page in dash.page_registry.values()
    ], pills=True,
    className="ms-4 gap-2"
)

app.layout = dbc.Container([
    dbc.Row([
        dbc.Col(html.H1(doc_title, className="text-light mt-4 mb-4 text-wrap"))
    ], className="bg-primary"),

    # Navbar
    dbc.Row([
        dbc.Col([inavbar], className="d-flex justify-content-start"),
        dbc.Col([theme_change], className="d-flex justify-content-end")
    ], className="inavbar bg-opacity-25 bg-secondary mt-4 mb-4"),

    dbc.Row([
        dbc.Col([
            dash.page_container
        ])
    ], className="bg-light")

], fluid=True)


if __name__ == "__main__":
    app.run()
