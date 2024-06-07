import seaborn as sns
import pandas as pd
import networkx as nx
import matplotlib.pyplot as plt
from seqanalysis.util.logging import config_logging
from IPython import embed
import plotly.express as px
from dash import Dash, dcc, html, callback, Output, Input, ctx, State
from dash.exceptions import PreventUpdate
import dash_cytoscape as cyto

log = config_logging()


plt.rcParams["svg.fonttype"] = (
    "none"  # this is so that svg figures save text as text and not the single letters
)


def plot_transition_matrix(matrix, labelx, labely, save_path, title):
    """
    Plot a heatmap of a transition matrix.

    Parameters:
    - matrix (array-like): The transition matrix to be visualized.
    - labels (list): Labels for the x and y axes.
    - save_path (str): File path to save the generated plot.
    - title (str): Title of the plot.
    """
    fig, ax = plt.subplots()
    hm = sns.heatmap(
        matrix,
        ax=ax,
        annot=True,
        vmin=0,
        vmax=100,
        fmt="d",
        cmap="Greys",
        xticklabels=labelx,
        yticklabels=labely,
    )
    ax.set_yticklabels(hm.get_yticklabels(), rotation=0)
    ax.set_xticklabels(hm.get_xticklabels(), rotation=45)
    ax.tick_params(left=False, bottom=False)
    sns.despine(top=False, right=False, left=False, bottom=False)
    ax.set_title(title)
    fig.tight_layout()
    log.info(f"Saving plot to {save_path}")
    fig.savefig(save_path, dpi=300)


def plot_transition_diagram(matrix, labels, node_size, matrix_labels):
    """
    Plot a transition diagram based on the given matrix and labels.

    Parameters:
    - matrix (array-like): The transition matrix to be visualized.
    - labels (list): Labels for the nodes in the diagram.
    - node_size (float): Size of the nodes in the diagram.
    - edge_width (float): Width scaling factor for the edges in the diagram.
    - save_path (str): File path to save the generated plot.
    - title (str): Title of the plot.
    """

    # Create a directed graph from the given matrix
    embed()
    fig = px.imshow(
        matrix,
        color_continuous_scale="Magma",
    )
    fig.update_traces(text=matrix_labels, texttemplate="%{text}")
    fig.update_xaxes(side="top")

    Graph = nx.from_numpy_array(matrix, create_using=nx.DiGraph)

    # Set the positions of nodes in a circular layout
    positions = nx.circular_layout(Graph)
    # mulitply by 10 to make the plot bigger
    positions = {node: [pos[0] * 10, pos[1] * 10] for node, pos in positions.items()}
    # get the x and y coordinates of the nodes

    nodes = [
        {
            "data": {
                "id": str(node),
                "label": label,
                "weight": w * 100,
                "height": w * 100,
            },
            "position": {"x": pos[0], "y": pos[1]},
            "classes": "bignode",
        }
        if w > 0.7
        else {
            "data": {"id": str(node), "label": label},
            "position": {"x": pos[0], "y": pos[1]},
            "classes": str(node),
        }
        for node, label, pos, w in zip(
            Graph.nodes, labels, positions.values(), node_size
        )
    ]

    edges = [
        {
            "data": {
                "source": str(source),
                "target": str(target),
                "weight": weight / 10,
                "label": str(weight),
            }
        }
        for source, target, weight in Graph.edges(data="weight")
    ]
    elements = nodes + edges

    # node_sizes = [
    #     {
    #         "selector": "." + str(node),
    #         "style": {
    #             "width": w * 100,
    #             "height": h * 100,
    #         },
    #     }
    #     for node, w, h in zip(*nodes, node_size, node_size)
    # ]
    app = Dash()
    app.layout = (
        html.Div(
            [
                html.Div(
                    className="TransitionDiagram",
                    children=[
                        cyto.Cytoscape(
                            id="cytoscape",
                            layout={"name": "circle"},
                            style={
                                "width": "49%",
                                "height": "90vh",
                                "display": "inline-block",
                            },
                            elements=elements,
                            stylesheet=[
                                # *node_sizes,
                                {
                                    "selector": "nodes",
                                    "style": {
                                        "label": "data(label)",
                                        "shape": "circle",
                                        "z-index": "10",
                                        "text-halign": "right",
                                        "width": "data(weight)",
                                        "height": "data(weight)",
                                    },
                                },
                                {
                                    "selector": "*",
                                    "style": {
                                        "target-distance-from-node": "10px",
                                        "source-distance-from-node": "10px",
                                        "z-index-compare": "manual",
                                        "text-margin-y": "-10px",
                                    },
                                },
                                {
                                    "selector": "edges",
                                    "style": {
                                        "curve-style": "bezier",
                                        "label": "data(label)",
                                        "target-arrow-shape": "triangle",
                                        "z-index": "1",
                                        "width": "data(weight)",
                                        "source-text-offset": 10,
                                    },
                                },
                                {
                                    "selector": ":loop",
                                    "style": {
                                        "control-point-step-size": 70,
                                    },
                                },
                            ],
                        ),
                        dcc.Graph(
                            figure=fig,
                            style={
                                "width": "50%",
                                "height": "90vh",
                                "display": "inline-block",
                            },
                        ),
                    ],
                ),
                html.Div(
                    [
                        html.P(id="cytoscape-tapNodeData-output"),
                        html.Div("Download graph:"),
                        html.Button("Save as PNG", id="btn-get-svg"),
                        html.Button(
                            "Remove Node", id="btn-remove-node", n_clicks_timestamp=0
                        ),
                    ],
                ),
            ],
        ),
    )

    @callback(
        Output("cytoscape", "generateImage"),
        Input("btn-get-svg", "n_clicks"),
    )
    def save_svg(btn):
        if btn is None:
            raise PreventUpdate
        else:
            return {"type": "png", "action": "download"}

    @callback(
        Output("cytoscape", "elements"),
        Input("cytoscape", "tapNodeData"),
        Input("btn-remove-node", "n_clicks"),
        State("cytoscape", "elements"),
    )
    def remove_node(data, n, elements):
        if n is None:
            raise PreventUpdate
        else:
            elements = [elem for elem in elements if elem["data"]["id"] != data["id"]]

        return elements

    # @callback(
    #     Output("cytoscape-tapNodeData-output", "children"),
    #     Input("cytoscape", "tapNodeData"),
    # )
    # def displayTapNodeData(data):
    #     if data:
    #         return f"You clicked on node {data['label']}"

    app.run_server(debug=True, use_reloader=False)

    exit()
