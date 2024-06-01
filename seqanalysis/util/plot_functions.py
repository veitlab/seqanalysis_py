import seaborn as sns
import pandas as pd
import networkx as nx
import matplotlib.pyplot as plt
from seqanalysis.util.logging import config_logging
from IPython import embed
from dash import Dash, dcc, html, callback, Output, Input, ctx
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


def plot_transition_diagram(matrix, labels, node_size, edge_width, save_path, title):
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
    embed()
    # Create a directed graph from the given matrix
    Graph = nx.from_numpy_array(matrix, create_using=nx.DiGraph)

    # Map node labels to corresponding nodes
    node_labels = dict(zip(Graph, labels))

    # Get edge labels from the graph
    edge_labels = nx.get_edge_attributes(Graph, "weight")

    # Set the positions of nodes in a circular layout
    positions = nx.circular_layout(Graph)
    # mulitply by 10 to make the plot bigger
    positions = {node: [pos[0] * 10, pos[1] * 10] for node, pos in positions.items()}
    # get the x and y coordinates of the nodes

    nodes = [
        {
            "data": {"id": str(node), "label": label},
            "position": {"x": pos[0], "y": pos[1]},
        }
        for node, label, pos in zip(Graph.nodes, labels, positions.values())
    ]
    edges = [
        {"data": {"source": str(source), "target": str(target), "weight": weight}}
        for source, target, weight in Graph.edges(data="weight")
    ]
    elements = nodes + edges
    styles = {
        "output": {
            "overflow-y": "scroll",
            "overflow-wrap": "break-word",
            "height": "calc(100% - 25px)",
            "border": "thin lightgrey solid",
        },
        "tab": {"height": "calc(98vh - 115px)"},
    }
    app = Dash()
    app.layout = (
        html.Div(
            [
                html.Div(
                    className="eight columns",
                    children=[
                        cyto.Cytoscape(
                            id="cytoscape",
                            layout={"name": "circle"},
                            style={"width": "calc(100%-500px)", "height": "95vh"},
                            elements=elements,
                            stylesheet=[
                                {
                                    "selector": "edges",
                                    "style": {
                                        "curve-style": "bezier",
                                        "label": "data(weight)",
                                        "target-arrow-shape": "triangle",
                                    },
                                },
                                {
                                    "selector": "nodes",
                                    "style": {
                                        "shape": "circle",
                                        "label": "data(label)",
                                    },
                                },
                            ],
                        )
                    ],
                ),
                html.Div("Download graph:"),
                html.Button("Save as SVG", id="btn-get-svg"),
            ],
        ),
    )

    @callback(
        Output("cytoscape", "generateImage"),
        Input("btn-get-svg", "n_clicks"),
    )
    def save_svg(n_clicks):
        if ctx.triggered:
            return {"type": "svg", "action": "download"}

    app.run_server(debug=True, use_reloader=False)
    exit()
