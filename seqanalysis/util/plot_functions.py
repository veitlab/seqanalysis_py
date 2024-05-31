import seaborn as sns
import pandas as pd
import networkx as nx
import matplotlib.pyplot as plt
from seqanalysis.util.logging import config_logging
import plotly.graph_objects as go
from IPython import embed
from dash import Dash, dcc, html
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

    x_pos, y_pos = zip(*positions.values())
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
    app = Dash()
    app.layout = html.Div(
        [
            cyto.Cytoscape(
                id="cytoscape",
                layout={"name": "circle"},
                style={"width": "100%", "height": "800px"},
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
        ]
    )
    app.run_server(debug=True, use_reloader=False)  #
    embed()
    exit()

    #
    # node_trace = go.Scatter(
    #     x=x_pos,
    #     y=y_pos,
    #     mode="markers",
    #     marker=dict(size=node_size, color="orange"),
    #     text=labels,
    #     hoverinfo="text",
    # )
    # edge_trace = go.Scatter(
    #     x=x_pos,
    #     y=y_pos,
    #     line=dict(width=edge_width, color="black"),
    #     hoverinfo="none",
    #     mode="lines",
    # )
    # fig = go.Figure(
    #     data=[edge_trace, node_trace],
    #     layout=go.Layout(
    #         title="<br>Network graph made with Python",
    #         titlefont_size=16,
    #         showlegend=False,
    #         hovermode="closest",
    #         margin=dict(b=20, l=5, r=5, t=40),
    #         annotations=[
    #             dict(
    #                 text="Transition Diagram",
    #                 showarrow=False,
    #                 xref="paper",
    #                 yref="paper",
    #                 x=0.005,
    #                 y=-0.002,
    #             )
    #         ],
    #         # xaxis=dict(showgrid=False, zeroline=False, showticklabels=False),
    #         # yaxis=dict(showgrid=False, zeroline=False, showticklabels=False),
    #     ),
    # )

    app.layout = html.Div([dcc.Graph(figure=fig)])

    app.run_server(debug=True, use_reloader=False)  #
    embed()
    exit()

    # Create a subplot with a specified size and margins
    fig, ax = plt.subplots(figsize=(21 / 2.54, 19 / 2.54))
    fig.subplots_adjust(left=0.05, bottom=0.05, right=0.95, top=0.95)

    # Draw nodes with specified size, color, and transparency
    nx.draw_networkx_nodes(
        Graph,
        pos=positions,
        node_size=node_size,
        node_color="tab:orange",
        ax=ax,
        alpha=0.9,
    )

    # Draw node labels
    nx.draw_networkx_labels(Graph, pos=positions, labels=node_labels)

    # Draw edges with specified width, arrows, and style
    edge_width = [x / edge_width for x in [*edge_labels.values()]]
    nx.draw_networkx_edges(
        Graph,
        pos=positions,
        node_size=node_size,
        width=edge_width,
        arrows=True,
        arrowsize=20,
        min_target_margin=25,
        min_source_margin=10,
        connectionstyle="arc3,rad=0.2",
        ax=ax,
    )

    # Draw edge labels at the midpoint of the edges
    nx.draw_networkx_edge_labels(
        Graph, positions, label_pos=0.5, edge_labels=edge_labels, ax=ax, rotate=False
    )

    # Remove spines for a cleaner appearance
    ax.spines["top"].set_visible(False)
    ax.spines["bottom"].set_visible(False)
    ax.spines["right"].set_visible(False)
    ax.spines["left"].set_visible(False)

    # Set the title of the plot
    plt.title(title)

    # Save the plot to the specified file path with a specified DPI
    log.info(f"Saving plot to {save_path}")
    fig.savefig(save_path, dpi=300)
