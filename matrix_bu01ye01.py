import os
import glob
import numpy as np
import matplotlib.pyplot as plt
from IPython import embed
import helper_functions as hf
import plot_functions as pf
import re

# -------------------------------------------------------------------------------------------------------------------
# paths
folder_path = 'D:/Birds/bu01ye01'
save_path = 'C:/Users/User/Documents/Doktorarbeit/data/bu01ye01/py_data/'

# -------------------------------------------------------------------------------------------------------------------
# folders to analyse
folder = ['screen', 'post_surgery']

# -------------------------------------------------------------------------------------------------------------------
#
bout_chunk = 'pjbj' # chunk has to be in a bout to be analysed as a bout
intro_notes = ['S', 'E', 'k', 'c']  # all intro nodes are replaced by 'i',
                                    # 'S' and 'E' have to be in the list so that the code works

# -------------------------------------------------------------------------------------------------------------------
# constants
edge_threshold = 5  # edges < 5 % aren't shown

# -------------------------------------------------------------------------------------------------------------------
# unique labels; list of labels for each folder
unique_labels = [['E', 'S', 'I', 'P', 'p', 'j', 'B', 'b', 's', 'H', 'h', 'e', 'G', 'g', 'd', 'F', 'f', 'l', 'c'],
                 ['E', 'S', 'I', 'P', 'p', 'j', 'B', 'b', 's', 'H', 'h', 'e', 'G', 'g', 'd', 'F', 'f', 'l', 'c', 'k', 'o']]

# -------------------------------------------------------------------------------------------------------------------
# double syllables, replaces_double_syllables
double_syl = [ 'bj']
rep_double_syl = ['bs']
chunks = ['pj', 'bs', 'he', 'gd', 'he', 'fl', 'i+', 'c+', 'ES']

# -------------------------------------------------------------------------------------------------------------------

# embed()
# quit()
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
    for i in range(len(double_syl)):
        all_bouts = re.sub(double_syl[i], rep_double_syl[i], all_bouts)

    unique_labels = sorted(list(set(all_bouts)))
    # np.save(save_path+folder[folderidx]+'_bouts', all_bouts)
    trans_matrix, trans_matrix_prob = hf.get_transition_matrix(unique_labels, all_bouts)
    node_size = np.round(np.sum(trans_matrix, axis=1)/np.min(np.sum(trans_matrix, axis=1)), 2) * 100

    trans_matrix_prob = hf.get_node_matrix(trans_matrix_prob, edge_threshold)

    'Plot Transition Matrix and Transition Diagram'
    node_labels = unique_labels[folderidx]
    # pf.plot_transition_matrix(trans_matrix_prob, node_labels)
    pf.plot_transition_diagram(trans_matrix_prob, node_labels, node_size)

    # ---------------------------------------------------------------------------------------------------------------

    bouts = hf.replace_chunks(all_bouts, chunks)

    # ToDo: put them in a good order, based on the best succession

    tm, tmp = hf.get_transition_matrix(unique_labels[folderidx], bouts)
    # tm = tm.astype(int)
    # embed()
    # quit()
    # node_size = np.round(np.sum(tm, axis=1)/np.min(np.sum(tm, axis=1)), 2)
    #
    # tmp = get_node_matrix(tmp, edge_threshold)

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

    tmpd = hf.get_node_matrix(tmpd, edge_threshold)

    ul = np.delete(unique_labels[folderidx], k)
    # ul = ['Start', 'i+', 'pj', 'p', 'bj', 'he', 'gd', 'fl']

    'Plot Transition Matrix and Transition Diagram'
    pf.plot_transition_matrix(tmpd, ul)
    pf.plot_transition_diagram(tmpd, ul, node_size*5)
    plt.title(folder[folderidx])

    os.chdir('..')


plt.show()

embed()
quit()
