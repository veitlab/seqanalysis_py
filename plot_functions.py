import seaborn as sns
import pandas as pd
import networkx as nx
import matplotlib.pyplot as plt

from IPython import embed


# im = ax.imshow(matrix)
# ax.set_xticks(list(range(len(matrix))))
# ax.set_yticks(list(range(len(matrix))))
# ax.set_xticklabels(labels)
# ax.set_yticklabels(labels)
# fig.colorbar(im)

def plot_transition_matrix(matrix, labels, save_path, title):
    fig, ax = plt.subplots()
    hm = sns.heatmap(matrix, ax = ax, annot=True, vmin=0, vmax=100, fmt="d", cmap='Greys',
               xticklabels=labels, yticklabels=labels)
    ax.set_yticklabels(hm.get_yticklabels(), rotation=0)
    ax.tick_params(left=False, bottom=False)
    sns.despine(top=False, right=False, left=False, bottom=False)
    ax.set_title(title)
    fig.savefig(save_path, dpi=300)


def plot_transition_diagram(matrix, labels, node_size, edge_width, save_path, title):
    fig, ax = plt.subplots(figsize=(21 / 2.54, 19 / 2.54))
    fig.subplots_adjust(left=0.05, bottom=0.05, right=0.95, top=0.95)
    Graph = nx.from_numpy_array(matrix, create_using=nx.DiGraph)

    node_labels = dict(zip(Graph, labels))
    # Graph = nx.relabel_nodes(Graph, node_labels)

    edge_labels = nx.get_edge_attributes(Graph, 'weight')
    positions = nx.circular_layout(Graph)

    # ToDo: nodes
    nx.draw_networkx_nodes(Graph, pos=positions, node_size=node_size, node_color="tab:orange", ax=ax, alpha=0.9)
    nx.draw_networkx_labels(Graph, pos=positions, labels=node_labels)
    # google networkx drawing to get better graphs with networkx
    # nx.draw(Graph, pos=positions, node_size=node_size, label=labels, with_labels=True, ax=ax)
    # ToDo: edges
    edge_width = [x / edge_width for x in [*edge_labels.values()]]
    embed()
    nx.draw_networkx_edges(Graph, pos=positions, node_size=node_size, width=edge_width,
                           arrows=True, arrowsize=20,
                           min_target_margin=25, min_source_margin=10, connectionstyle="arc3,rad=0.2",
                           ax=ax)
    nx.draw_networkx_edge_labels(Graph, positions, label_pos=0.5, edge_labels=edge_labels, ax=ax, rotate=False)

    ax.spines["top"].set_visible(False)
    ax.spines["bottom"].set_visible(False)
    ax.spines["right"].set_visible(False)
    ax.spines["left"].set_visible(False)

    plt.title(title)

    fig.savefig(save_path, dpi=300)
    #
    # embed()
    # quit()
