import os
import glob
import numpy as np
import matplotlib.pyplot as plt
from IPython import embed
import pandas as pd
import scipy.io as sio
import networkx as nx


def plot_transition_matrix(matrix, labels):
    fig, ax = plt.subplots()
    im = ax.imshow(matrix)
    ax.set_xticks(list(range(len(matrix))))
    ax.set_yticks(list(range(len(matrix))))
    ax.set_xticklabels(labels)
    ax.set_yticklabels(labels)
    fig.colorbar(im)
    fig.savefig('matrix.jpg', dpi=300)



def plot_transition_diagram(matrix, labels, node_size):
    fig, ax = plt.subplots(figsize=(21/2.54,19/2.54))
    fig.subplots_adjust(left=0.05, bottom=0.05, right=0.95, top=0.95)
    Graph = nx.from_numpy_matrix(matrix, create_using=nx.DiGraph)

    node_labels = dict(zip(Graph, labels))
    # Graph = nx.relabel_nodes(Graph, node_labels)

    edge_labels = nx.get_edge_attributes(Graph, 'weight')
    positions = nx.circular_layout(Graph)

    # ToDo: nodes
    nx.draw_networkx_nodes(Graph, pos=positions, node_size=node_size, ax=ax)
    nx.draw_networkx_labels(Graph, pos=positions, labels=node_labels)
    # google networkx drawing to get better graphs with networkx
    # nx.draw(Graph, pos=positions, node_size=node_size, label=labels, with_labels=True, ax=ax)
    # ToDo: edges
    edge_width = [x / 20 for x in [*edge_labels.values()]]
    nx.draw_networkx_edges(Graph, pos=positions, width=edge_width,
                           arrows=True, arrowsize=20,
                           min_target_margin=25, min_source_margin=10, connectionstyle="arc3,rad=0.1",
                           ax=ax)
    nx.draw_networkx_edge_labels(Graph, positions, label_pos=0.5, edge_labels=edge_labels, ax=ax, rotate=False)

    ax.spines["top"].set_visible(False)
    ax.spines["bottom"].set_visible(False)
    ax.spines["right"].set_visible(False)
    ax.spines["left"].set_visible(False)

    fig.savefig('graph.pdf', dpi=300)
