import os
import glob
import numpy as np
import matplotlib.pyplot as plt
from IPython import embed
import pandas as pd
import scipy.io as sio
import networkx as nx
from markovchain import MarkovChain
import graphviz as gv
import helper_functions as hf
import plot_functions as pf
import re

# -------------------------------------------------------------------------------------------------------------------
# paths
folder_path = 'D:/Birds/wh45wh40'
save_path = 'C:/Users/User/Documents/Doktorarbeit/data/wh45wh40/py_data/'

# -------------------------------------------------------------------------------------------------------------------
# folders to analyse
folder = ['recovery1']
# folder = ['screen', 'target_a2', 'recovery1']

# -------------------------------------------------------------------------------------------------------------------
#
bout_chunk = 'dmm' # chunk has to be in a bout to be analysed as a bout
intro_notes = ['S', 'E', 'l', 'k']  # all intro nodes are replaced by 'i',
                                    # 'S' and 'E' have to be in the list so that the code works

# -------------------------------------------------------------------------------------------------------------------
# constants
edge_threshold = 5  # edges < 5 % aren't shown

# -------------------------------------------------------------------------------------------------------------------
# unique labels; list of labels for each folder
unique_labels = [['S', 'I', 'D', 'd', 'm', 'o', 'A', 'a', 'f', 'c', 'L', 'K', 'H', 'h', 'r', 'e', 'g', 'j', 'S', 's', 'n', 'p', 'E'],
                 ['S', 'I', 'D', 'd', 'm', 'o', 'A', 'a', 'f', 'c', 'L', 'K', 'H', 'h', 'r', 'e', 'g', 'b', 'j', 'S', 's', 'n', 'p', 'E'],
                 ['S', 'I', 'D', 'd', 'm', 'o', 'A', 'a', 'f', 'c', 'L', 'K', 'H', 'h', 'r', 'e', 'g', 'j', 'S', 's', 'n', 'p', 'E']]

# -------------------------------------------------------------------------------------------------------------------
# double syllables, replaces_double_syllables
double_syl = ['mm', 'nn', 'hhc', 'jj', 'bj']
rep_double_syl = ['mo', 'np', 'hre', 'js', 'bs']
chunks = ['dmo', 'afc', 'hreg', 'snp', 'i+', 'l+', 'k+', 'ES']

# -------------------------------------------------------------------------------------------------------------------


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
        bouts, _ = hf.get_bouts(seqs, bout_chunk)
        all_bouts.append(bouts)
        os.chdir('..')

    # ---------------------------------------------------------------------------------------------------------------

    all_bouts = ''.join(set(all_bouts))
    # unique_labels = sorted(list(set(all_bouts)))

    for i in range(len(double_syl)):
        all_bouts = re.sub(double_syl[i], rep_double_syl[i], all_bouts)

    np.save(save_path+folder[folderidx]+'_bouts', all_bouts)
    # trans_matrix, trans_matrix_prob = hf.get_transition_matrix(unique_labels[folderidx], all_bouts)
    # node_size = np.round(np.sum(trans_matrix, axis=1)/np.min(np.sum(trans_matrix, axis=1)), 2) * 100
    #
    # trans_matrix_prob = np.around(trans_matrix_prob, 2) * 100
    # trans_matrix_prob = trans_matrix_prob.astype(int)
    # trans_matrix_prob[trans_matrix_prob < edge_threshold] = 0
    #
    # 'Plot Transition Matrix and Transition Diagram'
    # node_labels = unique_labels[folderidx]
    # # pf.plot_transition_matrix(trans_matrix_prob, node_labels)
    # pf.plot_transition_diagram(trans_matrix_prob, node_labels, node_size)

    # ---------------------------------------------------------------------------------------------------------------

    bouts = hf.replace_chunks(all_bouts, chunks)

    # ToDo: put them in a good order, based on the best succession
    # unique_labels = sorted(list(set(bouts)))
    # unique_labels = ['S', 'I', 'D', 'd', 'm', 'o', 'A', 'a', 'f', 'c', 'L', 'K', 'H', 'h', 'r', 'e', 'g', 'b', 'j', 'S', 's', 'n', 'p', 'E']

    tm, tmp = hf.get_transition_matrix(unique_labels[folderidx], bouts)
    tm = tm.astype(int)

    # node_size = np.round(np.sum(tm, axis=1)/np.min(np.sum(tm, axis=1)), 2)
    #
    # tmp = np.around(tmp, 2) * 100
    # tmp = tmp.astype(int)
    # tmp[tmp < edge_threshold] = 0

    'Plot Transition Matrix and Transition Diagram'
    # pf.plot_transition_matrix(tmp, unique_labels)
    # pf.plot_transition_diagram(tmp, unique_labels, node_size)

    # embed()
    # quit()
    # ---------------------------------------------------------------------------------------------------------------

    'Plot Transition Matrix and Transition Diagram'
    # pf.plot_transition_matrix(tm, unique_labels)

    k = np.where(sum(tm) / sum(sum(tm)) * 100 <= 1)
    tmd = np.delete(tm, k, axis=1)
    tmd = np.delete(tmd, k, axis=0)

    node_size = np.round(np.sum(tmd, axis=1)/np.min(np.sum(tmd, axis=1)), 2)*100
    tmpd = (tmd.T / np.sum(tmd, axis=1)).T

    tmpd = np.around(tmpd, 2) * 100
    tmpd = tmpd.astype(int)
    tmpd[tmpd < edge_threshold] = 0

    ul = np.delete(unique_labels[folderidx], k)
    # ul = ['i+', 'dmm', 'afc', 'l+', 'k+', 'hhcg', 'b', 'jnn', 'Start']

    'Plot Transition Matrix and Transition Diagram'
    # pf.plot_transition_matrix(tmpd, ul)
    pf.plot_transition_diagram(tmpd, ul, node_size*5)
    plt.title(folder[folderidx])

    os.chdir('..')


plt.show()

embed()
quit()



