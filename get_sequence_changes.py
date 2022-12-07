import re
import os
import glob
import numpy as np
import pandas as pd
import networkx as nx
import plot_functions as pf
import helper_functions as hf
import matplotlib.pyplot as plt

from IPython import embed

# folder_path = 'D:/Birds/wh113gr7'
folder_path = 'D:/Birds/wh45wh40'
chunks = ['dmm', 'afc', 'hhcg', 'jjnn', 'i+']
intro_notes = ['S', 'E', 'l', 'k']
edge_threshold = 5  # edges < 5 % aren't shown
# move_path = 'C:/Users/User/Documents/Doktorarbeit/data/wh09pk39'
folder = ['screen', 'target_a2']

sequence_length = 10
seq_start = ['gj', 'gb']


os.chdir(folder_path)


for folderidx in range(len(folder)):
    os.chdir(folder[folderidx])
    print(os.getcwd())
    days = [d for d in os.listdir(os.getcwd()) if os.path.isdir(d)]

    all_bouts = []

    for dayidx in range(len(days)):
        os.chdir(days[dayidx])
        print(os.getcwd())
        seqs = hf.get_labels(intro_notes)
        bouts, _ = hf.get_bouts(seqs, chunks[0])
        all_bouts.append(bouts)
        os.chdir('..')

    # ---------------------------------------------------------------------------------------------------------------

    all_bouts = ''.join(set(all_bouts))
    reg_exp = r'' + seq_start[folderidx] + '{1}(' + '.' * sequence_length + ')'

    follow_seq = re.findall(reg_exp, all_bouts)

    for follseq_idx in range(len(follow_seq)):
        follow_seq[follseq_idx] = re.split('S', follow_seq[follseq_idx])[0]
        follow_seq[follseq_idx] = re.split('b', follow_seq[follseq_idx])[0]

    # tabelle =  [[source], [target], [source_syl], [target_syl], [1]]
    tabelle = [[], [], [], [], []]
    for i in range(len(follow_seq)):
        for j in range(len(follow_seq[i]) - 1):
            tabelle[0].append(j)
            tabelle[1].append(j + 1)
            tabelle[2].append(follow_seq[i][j])
            tabelle[3].append(follow_seq[i][j + 1])
            tabelle[4].append(1)

    # make pandas Data Frame to get value of how many connections are between specific source and target syllable
    tpd = pd.DataFrame(np.array(tabelle).T, columns=['A', 'B', 'C', 'D', 'E'])
    tpiv = pd.pivot_table(tpd, index=['A', 'B', 'C', 'D'], values='E', aggfunc=len)
    DataTab = tpiv.reset_index()

    # to normalize on the sequence position
    for i in range(len(np.unique(DataTab.A))):
        DataTab.E[DataTab.A == str(i)] = DataTab.E[DataTab.A == str(i)] / np.sum(DataTab.E[DataTab.A == str(i)]) * 100

    d = {'source': np.array(DataTab.A + DataTab.C),
         'target': np.array(DataTab.B + DataTab.D),
         'value': np.array(DataTab.E)}

    df = pd.DataFrame(data=d)
    unique_source_target = list(pd.unique(df[['source', 'target']].values.ravel('K')))
    mapping_dict = {k: v for v, k in enumerate(unique_source_target)}

    df['source'] = df['source'].map(mapping_dict)
    df['target'] = df['target'].map(mapping_dict)

    # make final things or graph
    # ToDo: node
    node_pos = hf.get_node_positions(unique_source_target)
    node_label = [string[1] for string in unique_source_target]
    node_idx = [*range(len(node_label))]
    node_label_dict = {n: nodeL for n, nodeL in zip(node_idx, node_label)}

    node_size = []
    for i in range(len(node_pos)):
        if np.any(df.source == i):
            node_size.append(np.ceil(df.value[df.source == i].sum()) * 50)
        elif np.any(df.target == i):
            node_size.append(np.ceil(df.value[df.target == i].sum()) * 50)

    # ToDo: egdes
    edge_label_dict = dict()
    for n, edgeL in zip(np.column_stack((df.source, df.target)), np.round(df.value, 2)):
        if edgeL >= 1:
            edge_label_dict[tuple(n)] = str(edgeL)

    # ToDo: Graph
    G = nx.Graph()
    G.add_nodes_from(node_idx)
    G.add_edges_from(np.column_stack((df.source, df.target)))

    fig, ax = plt.subplots()
    nx.draw_networkx_nodes(G, pos=node_pos, node_size=node_size, ax=ax)
    nx.draw_networkx_edges(G, pos=node_pos, ax=ax)
    nx.draw_networkx_labels(G, pos=node_pos, labels=node_label_dict)
    nx.draw_networkx_edge_labels(G, pos=node_pos, label_pos=0.5, edge_labels=edge_label_dict, ax=ax, rotate=False)
    plt.show()

    # embed()
    # quit()

    os.chdir('..')

