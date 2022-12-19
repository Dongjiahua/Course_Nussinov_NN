from algo import Nussinov
import numpy as np
import dash
from dash.dependencies import Input, Output, State
import dash_bio as dashbio
from dash import dcc, html,dash_table
from dash.exceptions import PreventUpdate
from table import example
import pandas as pd
import os
os.environ["MKL_THREADING_LAYER"]="GNU"
nu = Nussinov()
import subprocess
import tempfile

def build_seq(rna, stru):
    seq = dict()
    for st in stru:
        seq[st]={
            "sequence":rna,
            "structure":st
        }
    return seq 
ini_seq = build_seq(nu.rna,nu.stru) 

app = dash.Dash(__name__)

app.layout = html.Div([
html.Div([
    html.Div(
        children=[
            html.Label('Input your RNA sequence.'),
            html.Br(),
            dcc.Input(id="rna_input", value='GGUCCAC', type='text',style={"width":"90%"}),
            html.Br(),
            html.Button('Submit', id='submit-val', n_clicks=0),
            html.Br(),
            html.Br(),
            html.Label('Select the structure to display below.'),
            dcc.Dropdown(
                id='my-default-forna-sequence-display',
                options=[
                    {'label': name, 'value': name} for name in ini_seq.keys()
                ],
                multi=True,
                value=["((..).)"]
            ),
            html.Br(),
            html.Br(),
            html.Label('Select the structure to backtrace'),
            dcc.Dropdown(
                id='table_choice',
                options=[
                    {'label': name, 'value': name} for name in ini_seq.keys()
                ],
                multi=False,
                value=["((..).)"]
            ),
        ],style={'padding': 10, 'flex': 1}
    ),
    html.Div(
        children=[
            dash_table.DataTable(editable=False,id="dp_table")
        ],style={'padding': 20, 'flex': 1}
    ),

    ], style={'display': 'flex', 'flex-direction': 'row'}),

    
    html.Div(
        children=[
            html.Div(children=[
                html.Label('Selected Structure visualization'),
                html.Div(children=[
                    dashbio.FornaContainer(id='my-default-forna')],style={"border":"2px black solid"})
                        ],style={'padding': 10, 'flex': 1}),
            html.Div(children=[
                html.Label('MXFOLD2 Predicted'),
                html.Div(children=[
                dashbio.FornaContainer(id='mxfold2-forna')],style={"border":"2px black solid"})
                    ],style={'padding': 10, 'flex': 1}),
        ], style={'display': 'flex', 'flex-direction': 'row'}
    )
]
)
# app.layout = html.Div([
#     html.Label('Input your RNA sequence below.'),
#     dcc.Input(id="rna_input", value='GGUCCAC', type='text'),
#     html.Button('Submit', id='submit-val', n_clicks=0),
#     html.Label('Select the sequences to display below.'),
#     dcc.Dropdown(
#         id='my-default-forna-sequence-display',
#         options=[
#             {'label': name, 'value': name} for name in ini_seq.keys()
#         ],
#         multi=True,
#         value=["((..).)"]
#     ),
#     dashbio.FornaContainer(id='my-default-forna'),
#     dcc.Dropdown(
#         id='table_choice',
#         options=[
#             {'label': name, 'value': name} for name in ini_seq.keys()
#         ],
#         multi=False,
#         value=["((..).)"]
#     ),
#     dash_table.DataTable(sort_action='native',editable=False,id="dp_table")

# ])

@app.callback(
    Output('my-default-forna', 'sequences'),
    Input('my-default-forna-sequence-display', 'value'),
)
def show_selected_sequences(value):
    ini_seq = build_seq(nu.rna,nu.stru) 
    if value is None:
        raise PreventUpdate
    return [
        ini_seq[selected_sequence]
        for selected_sequence in value
    ]
    
@app.callback(
    Output('my-default-forna-sequence-display', 'options'),
    Output('table_choice', 'options'),
    Output('mxfold2-forna', 'sequences'),
    Output('table_choice',"value"),
    Output('my-default-forna-sequence-display',"value"),
    Input('submit-val', 'n_clicks'),
    State('rna_input', 'value')
)
def show_selected_sequences(n_clicks,value):
    rna=value.upper()
    for s in rna:
        if s not in ["U","A","C","G"]:
            raise dash.exceptions.PreventUpdate
    nu.rna = rna
    sc, nu.matrix = nu.solve(rna)
    nu.stru, nu.trace = nu.backward(rna, nu.matrix)
    options=[{'label':opt, 'value':opt} for opt in nu.stru]

    with tempfile.NamedTemporaryFile(suffix='.fa') as tempf:
        tempf.write(f">S\n{rna}".encode())
        tempf.seek(0)
        # os.system("mxfold2 predict test.fa")
        output = subprocess.check_output(f'mxfold2 predict {tempf.name}', shell=True).decode()
        # print(rna)
        # print(output)
        structure = output.split("\n")[2].split(" ")[0]
    seq = [{
            "sequence":rna,
            "structure":structure
    }]

    return options, options,seq,[nu.stru[0]],[nu.stru[0]]

@app.callback(
    Output("dp_table","columns"),
    Output("dp_table","data"),
    Output("dp_table","style_data_conditional"),
    Input("table_choice","value")
)
def show_dp_table(value):
    print(value)
    matrix = nu.matrix.tolist()
    
    rna = nu.rna 
    # print(matrix)
    for i in range(len(matrix)):
        matrix[i].insert(0,rna[i]) 
    data=matrix
    df = pd.DataFrame(data)

    # df['id'] = df.index
    
    columns = [{"name":"RNA", "id":"dp_s"}]
    columns += [{"name":rna[i], "id":f"dp_{i}"} for i in range(len(rna))]
    df_columns = [c["id"] for c in columns]
    # print(df)
    df.columns = df_columns
    ret_data = df.to_dict("records")
    if type(value)==list and len(value)>0:
        value = value[0]
    cells = nu.trace[value]
    condition=[]
    for cell in cells:
        condition.append(
            {
                'if': {
                    'row_index': cell[0],
                    'column_id': f'dp_{cell[1]}'
                },
                'backgroundColor': 'tomato',
                'color': 'white'
            }
        )
    # print(nu.trace)
    return columns, ret_data, condition

if __name__ == '__main__':
    app.run_server(debug=True)
    print(1)