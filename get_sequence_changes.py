import re
import os
import glob
import yaml
import numpy as np
import pandas as pd
import networkx as nx
import plot_functions as pf
import helper_functions as hf
import matplotlib.pyplot as plt

from IPython import embed

folder_path = 'D:/Birds/wh45wh40'
chunks = ['dmm', 'afc', 'hhcg', 'jjnn', 'i+']
intro_notes = ['S', 'E', 'l', 'k']
edge_threshold = 5  # edges < 5 % aren't shown
# move_path = 'C:/Users/User/Documents/Doktorarbeit/data/wh09pk39'

# ---------------------------------------------------------------------------------------------------------------
bout_221027 = 'YI3b18gd9flW0b18gd9flW0b18gd9flW0b18gd9flW0b18gd9flW0b18gdpk5W0pkYIW0b18g73b18gd9flW0b18gU5W0b1pYIW0b18gd9flW0b18gd9flW0b18gd9flW0b18gd9flW0pb18gd9flW0b18gd9flW0b18gdp5W0b18gdZYIW0b18gd9flW0b18gd9flW0b18gd9flW0b18gd9flW0b18gd9flW0b18gdpb18flW06ZYIW0b18gd9flW0b18gd9flW0b18gd9flW0b18gd9flW0b18gdpYIW0b18gd9flW0b18gd9flW0b18gd9flW0b18gd9flW0b18gdp'
unique_labels_pre = sorted(list(set(bout_221027)))

bout_221027 = re.sub('9', 'h', bout_221027)
bout_221027 = re.sub('8', 'H', bout_221027)

bout_221027 = re.sub('7', 'C', bout_221027)
bout_221027 = re.sub('6', 'C', bout_221027)
bout_221027 = re.sub('5', 'C', bout_221027)

bout_221027 = re.sub('3', 'B', bout_221027)
bout_221027 = re.sub('1', 'J', bout_221027)
bout_221027 = re.sub('0', 'j', bout_221027)

bout_221027 = re.sub('Z', 'K', bout_221027)

bout_221027 = re.sub('lW', 'lA', bout_221027)
bout_221027 = re.sub('W', 'p', bout_221027)
bout_221027 = re.sub('U', 'P', bout_221027)

labels_post_221027 = ['Y', 'I', 'p', 'j', 'b', 'J', 'H', 'g', 'd', 'h', 'f', 'l', 'A', 'P', 'C', 'K', 'B']
tm_221027, _ = hf.get_transition_matrix(bout_221027, labels_post_221027)
embed()
# ---------------------------------------------------------------------------------------------------------------
k = np.where(sum(tm_221027) / sum(sum(tm_221027)) * 100 <= 0.1)
tmd = np.delete(tm_221027, k, axis=1)
tmd = np.delete(tmd, k, axis=0)
print(np.delete(labels_post_221027, k))

node_labels = np.delete(labels_post_221027, k)
node_labels = ['Y', 'i+', 'p1', 'j1', 'b', 'j2', 'he1', 'g', 'd', 'he2', 'f', 'l', 'p2', 'p3', 'c+', 'k+', 'j3']
tmpd = (tmd.T / np.sum(tmd, axis=1)).T
tmpd = hf.get_node_matrix(tmpd, 3)
node_size = np.round(np.sum(tmd, axis=1) / np.min(np.sum(tmd, axis=1)), 2) * 100
'Plot Transition Matrix and Transition Diagram'
# pf.plot_transition_matrix(tmpd,
#                           node_labels,
#                           'C:/Users/User/Documents/Doktorarbeit/data/bu01ye01/py_data/' + '221027' + '_matrix.pdf',
#                           '221027')
pf.plot_transition_diagram(tmpd,
                           node_labels,
                           node_size,
                           10,
                           'C:/Users/User/Documents/Doktorarbeit/data/bu01ye01/py_data/' + '221027' + '_graph.pdf',
                           '221027')

# ---------------------------------------------------------------------------------------------------------------

bout_screen = 'YI45b68gd9fl25b68gd9fl25b68gd9fl25b68gd9fl25b68gd9fl25b68gjb68gdCYI45b68gd9fl25b68gd9fl25b68gd9fl25b68gd9fl25b68gd9fl25b68gd9fl25b68gd9fl25b68gd3CYI45b68gd9fl25b68gd9fl25b68gd9fl25b68gd3b68gdYI45b68gd9fl25b68gd9fl25b68gd9fl25b68gd9fl25b68gd9fl25b68gd3b68gpYI45b6b68gd9fl25b68gd9fl25b68gd9fl25b68gd3b68gd3b68gdjYI45p5b68gd9fl25b68gd9fl25b68gpb68gd9jYI45b68gd9fl25b68gd9fl25b68gd9fl25b68gd9fl25b68gd9fl25b68gd3b68gd9jpb68gdjYI45b68gd9fl25b68gd9fl25b68gd9fl25b68gd3b68gd9fl25b68gd3b68gdYI45p5b68gd9fl25b68gd9fl25b68gd9fl25b68gd9fl25b68gd9fl25b68gd9fl25b68gd3b68gd3YI45p5b68gd9fl25b68gd9fl25b68gd9fl25b68gd9fl25b68gd9fl25b68gd9fl2YI45p5b68gd9fl25b68gd9fl25b68gd3b68gdCp5b68gd3b68gdCp5b68gd3b68gd9fl25b68gd3YI45b68gd9fl25b68gd9fl25b68gjb68gd9fl25b68gd9fl25b68gjKYI45b68gd9fl25b68gd9fl25b68gd9fl25b68gd9fl25b68gd9fl25b68gd3b68gd9fl25b68gjb68gdCp5b68gdjYI45p5b68gd9fl25b68gd9fl25b68gd9fl25b68gd9fl25b68gd9fl25b68gdYI45p5b68gd9fl25b68gd9fl25b68gpb68gd9fl25b68gd3b68gdCp5p5b68gdjYI45p5p5b68gd9fl25b68gd9fl25b68gd9fl25b68gd9fl25b68gd9fl25b68gd3b68gd3b68gdjYI45b68gjb68gd9fl25b68gd9fl25b68gd9fl25b68gd9fl25b68gd3b68gjb68gjYI45b68gd9fl25b68gd9fl25b68gd9fl25b68gd9fl25b68gd9flYI45b68gd9fl25b68gd9fl25b68gd9fl25b68gd9fl25b68flCp5pb68gdCYI45p5p5b68gd9fl25b68gd9fl25b68gd9fl25b68gd3b68gd3b68gjYI45p5b68gd9fl25b68gd9fl25b68gd9fl25b68gd9fl25b68gd3b68gCK7p5b68gdYI45p5b68gd9fl25b68gd9fl25b68gd9fl25b68gd9fl25b68gd3b68gd9fl25b68gjKYI45b68gd9fl25b68gd9fl25b68gd9fl25b68gd9fl25b68gd9fl25b68gd35b68gd9fl25b68gjK7YI45b68gd9fl25b68gd9fl25b68gd9fl25b68gdYI45p5b68gd9fl25b68gd9fl25b68fl25b68gd9fl25b68gd9fl25b68gd3b68gd9fl25b68gd3CYI45p5b68gd9fl25b68gd9fl25b68gd9fl25b68gd9fl25b68gd9fl25b68gd3b68gd3YIgjp5b68gd9fl25b68gd9fl25b68gd9fl25b68gd9fl25b68gd9fl25b68gd9fl25b68gd9fl25b68gd3b68gCYI45b68gd9fl25b68gd9fl25b68gd9fl25b68gd9fl25b68gd9fl25b68gd9fl25b68gd9fl25b68gd3b68gdCKYI45p5b68gd9fl25b68gd9fl25b68gd9fl25b68gd9fl25b68gd9fl25b68gd3b68gdCKYI45p5b68gd9fl25b68gd9fl25b68fl25b68gd3b68gdCb68gd9fl25b68gd9fl25b68gpK7YI45p5b68gd9fl25b68gd9fl25b68gd9fl25b68gd9fl25b68gd9fl25b68gd9fl25b68gd3b68gCK7p5b68gdjYI45b68gd9fl25b68gd9fl25b68gd9fl25b68gd9fl25b68gd3b68gdCKYI45b68gd9fl25b68gd9fl25b68gd9fl25b68gd358gd9fl25b68gdCK7p5b68gdjYI45b68gd9fl25b68gd9fl25b68fl25b68gd9fl25b68gd9fl25b68gd9fl25b68gd3CYI45p5b68gd9fl25b68gd9fl25b68gd9fCYI45b68gd9fl25b68gd9fl25b68gd3Cp5b68gd9fl25pb68gd9fl25b68gd9fl25b68gd9fl25b68gdYI45p5b68gd9fl25b68gd9fl25b68gd9flCp5b68gd3b68gd9fl25b68gd9fl25b68gdCp5b68gpYI45p5b68gd9fl25b68gd3b68gd9fl25b68gd9fl25b68gdCp5b68gd3b68gd9fl25b68gd3b68gdbYI45p5b68gd9fl25b68gd3b68gdCp5b68gd9fl25b68gd9fl25b68gd3b68gd3YI45b68gd9fl25b68gd9fl25b68gd9fl25b68gd9fl25b68gd9fl25b68gd3b68gjb68fl25b68gdjCYI45b68gj8gd9fl25b68gd9fl25b68gd9fl25b68gd9fl25b68gd3b68gd9fl25b68gd3b68gdjYI45p5b68gd9fl25b68gd9fl25b68gd9fl25b68gd3b68gd9fl25b68gjKYI45b68gd9fl25b68gd9fl25b68gd9fl25b68gd9fl25b68gdCp5b68gd9fl25b68gdjYI45p5b68gd9fl25b68gd9fl25b68gd9fl25b68gjYI45p5b68gd9fl25b68gd9fl25b68gd9fl25b68gd3b68gd9fl25b68gd9fl25b68gd9fl25b68gd9fl25b68gdCKYI45b68gd9fl25b68gd9fl25b68gd9jd9fl25b68gd9fl25b68gd9fl25b68gd9fl25b68gdCKYI45p5b68gd9fl25b68gd9fl25b68gd9fl25b68gd9fl25b68gd3b68gd9fl25b68gd9fl25b68gd9fl25b68gjKYI45b68gdjp5b68gd9fl25b68gd9fl25b68gd9fl25b68gd9fl25p5b68gd3b68gd3YI45b68gd9fl25b68gd9fl25b68gd9fl25b68gd9fl25b68gd9fl25b68gd9fl25b68gd3b68gp5YI45p5b68gd9fl25b68fl25b68gd9fl25b68gd9fl25b68gd9fl25b68gd9fl25b6b68gd9fl25b68gdCYI45p5b68gd9fl25b68gd9fl25b68gd9fl25b68gd9fl25b68gd9fl25b68gd9fl25b68gdCKYI45b68gd9fl25b68gd9fl25b68gd9fl25b68gd9fl25b68gd9fl25b68gdCp5b68gdCb68gd9fl25b68gdCKYI45b68gd9fl25b68gd9fl25b68gd9fl25b68gd9fl25b68gd9fl25b68gd3b68gd9fl25b68gdjYI45b68gd9fl25b68fl25b68fl25b68gd9fl25b68gd9fl25b68gd9fl25b68gdCp5b68gd9fl25b68gdCKYI45p5b68gd9fl25b68gd9fl25b68gd9fl25b68gd9fl25b68gd9fl25b68gd9fl25b68gdCp5b68gdCKYI45b68gd9fl25b68gjYI45b68gd9fl25b68gd9fl25b68gd9fl25b68gd9fl25b68gd9fl25b68gd9fl25b68gCKYI45p5p5b68gd9fl25b68gpb68gCp5b68gd9fl25b68gd9fl25b68gdCYI45p5b68gd9fl25b68fl25b68gd9fl25b68gd9fl25b68gd9fl25b68gd3b68gd9fl25b68gjYI45p5p5b68gd9fl25b68gd9fl25b68gd9fl25b68gd9fl25b68fl25b68gd9fl25b68gd3b68gdCKYI45p5b68gd9fCp5b6b68gd9fl25b68gd9fl25b68gd9fl25b68gd9fl25b68gd9fl25b68gdYI45p5b68gd9fl25b68gd9fl25b68gd9fl25b68fl25b68gd9fl25b68gd9fl25b68gd9fl25b68gdCKYI45p5b68gd9fl25b68gd9fl25b68gd9fl25b68gd9fl25b68gd9fl25b68gd9fl25b68gd9fl25b68gd9flYI45p5b68gd9fl25b68gd9fl25b68gd9gd9fl25b68gd9fl25b68gd9fl25b68gd9fl25b68gd9fl25b68gd9fl25bYI45p5b68gd9fl25b68gd9fl25b68gd9fl25b68gd9fl25b68gd9fl25b68gd9fl25b68gd9fl25b68gd9fl25b68CYI45b68gd9fl25b68gd9fl25b68gdCp5b68gd9fl25b68gd9fl25b68gd9fl25pYI45p5b68gd9fl25b68gd9fl25b68gd9fl25b68gd9fl25b68gdCKYI45p5b68gd9fl25b68gd9fl25b68gd9fl25b68gd9fl25b68gd9fl25b68gd9fl25b68gd3b68gd9fl25b68gjYI45p5b68gd9fl25b68gd9fl25b68gd9fl25b68gd9fl25b68gdCp5b68fYI45b68gd9fl25b68gd9fl25b68gd9fl25b68gd9fl25b68gd3b68gdYI45b68gd9fl25b68gd9fl25b68gd9fl25b68gd9gd9fl25b68gd3b68gdCp5b68gd3b68gd9fl25b68gdCp5b68gjYI45p5b68gd9fl25b68gd9fl25b68gd9fl25b68gd9fl25b68gd9fl25b68gd3b68gdCKYI45p5b68gd9fl25b68gd9fl25b68gd9fl25b68gd9fl25b68gd9fl25b68gd9gdCp5b68gd3b68gdClYI45b68gd9fl25b68gd9fl25b68gd9fl25b68gd9fl25b68gd9fl25b68gd9fl25b68gd9fl25b68gdCKYI45b68gd9fl25b68gd9fl25b68gd9fl25b68gd9fl25b68gd9fl25b68gd9fl25b68gdCK7p5b68gd9fl25b68gd9fljYI45b68gd9fl25b68gd9fl25b68gd9fl25b68gd9fl25b68gd9fl25b68gd9fl25b68gd9fl25b68gd9fl25b68gd9fpKYI458gd9fl25b68gd9fl25b68gd9fl25b68gd9fl25b68gd9fl25b68gd9fl25b68gd9fl25b68gd9fl25b68gdCYI45p5b68gd9fl25b68gd9fl25b68gd9fl25b68gd9fl25b68gd9gd3K7p5b68gd9fl25b68gd3b68gd3b68gdCKYI45b68gd3b68gd9fl25b68gd9fl25b68gd9fl25b68gd9fl25b68gd9fl25b68gd9fl25b68gYI45p5b68gd9fl25b68gd9fl25b68gd9fl25b68gd9fl25b68gd3b68gd9fjYI45p5b68gd9fl25b68fl25b68gd9fl25b68gd9fl25b68gd9fl25b68gd9fl25b68gd9fl25b68gd3b68gpYI45b68gd9fl25b68gd9fl25b68gd9fl25b68gd9fl25b68gd9fl25b68gd3b68gd9fl25b68gjK7p5b68gd9fl25b68gdCp5b68gdjYI45b68gd9fl25b68gd9fl25b68gd9fl25b68gd9fl25b68gjb68gd9fl25b68gdCp5b68gdjYI45b68gd9fl25b68gd9fl25b68gd9fl25pCp5b68gdjYI45b68gd9fl25b68gd9fl25b68gd9fl25b68gd9fl25b68gd9fl25b68gd9fl25b68gdCpb68gd9fl25b6pb68gd9fl25CYI45b68gd9fl25b68gd9fl25b68gd9fl25b68gd9fl25b68gd9fl25b68gdCK7p5b68gd9fl25b68gdCYI45p5b68gd9fl25b68gd9fl25b68gd9fl25b68gd9fl25b68gd9gCYIlCp5p5b68gd9fl25b68gd9fl25b68gd9fl25b68gd3b68gd9fl25b68gdCKYI45p5b68gd9fl25b68gd9fl25b68gd9fl25b68gd9fl25b68gd3b68gdYI45b68gd9fl25b68gd9fl25b68gd9fl25b68gd9fl25b68gd9fl25b68gd9fl25b68gd9fl25b68gdCp5b68gd9fl25b68gdCKYI45b68gd9fl25b68gd9fl25b68gd9fl25b68gd9fl25bKYI45b68gd9fl25b68gd9fl25b68gd9fl25b68gd9fl25b68gd9fl25b68gd9fl25b68gd9fl25b68gd9fl25b68gd9fl25b68gjKYI45p5b68gd9fl25b68gd9fl25b68gd9fl25b68gd9fl25b68gd9fl25b68gdCK7p5b68gdCKYI45p5b68gd9fl25b68gd9fl25b68gd9fl25b68gd9fl25b68gd9fl25b68gd9fl25b68gdCKYI45b68gd9fl25b68gd9fl25b68gd9fl25b68gd9fl25b68gd9fl25b68gd3KYI45p5b68gd9fl25b68gd9fl25b68gd9fl25b68gd3b68gd9fl25b68gd3KYI45p5b68gd9fl25b68gd9fl25b68gd9fl25b68gd9fl25b68gd9fl25b68gd9fl25b68gd9fl25b68gd9fl25b68gd9fl25b68gdCKYI45b68gd9fl25b68gd9fl25b68gd9fl25b68gd3b68gd9fl25b68gd9fl25b68gd3b68gdCYI45p5b68gd9fl25b68gd9fl25b68gd9fl25b68gd9fl25b68gd9fl25b68gd3b68gd9fl25b68gd9fl25pCKYI45p5b68gd9fl25b68gd9fl25b68gd3b68gd9fl25b68gd9fl25b68gd9fl25b68gd9fl25b68gd9fl25b68gdCKYI45p5b68gd9fl25b68gd9fl25b68gd9fl25b68gd9fl25b68gd9fl25b68gd9fl25b68gd9fl25b68gdCKYI45b68gd9fl25b68gd9fl25b68gd9fl25b68gd9fl25b68gd9fl25b68gd9fl25b68gd9fl25b68gd9fl25b68gdCKYI45p5b68gd9fl25b68gd9fl25b68gd9fl25b68gd9fl25b68gd9gjK7p5b68gd9fl25b68gd9fl25b68gd9fl25b68gpKYI45b68gd9fl25b68gd9fl25b68gd9fl25b68gd9fl25b68gd9fl25pYI45p5b68gd9fl25b68gd9fl25b68gd9fl25b68gd9fl25b68gd9fl25b68gd9fl25b68gd9fl25b68gjYI45p5p5b68gd9fl25b68gd9fl25b68fl25b68gd9fl25b68gd9fl25b68gd9fl25b68gdCKYI45p5b68gd9fl25b68gd9fl25b68gd9gd9fl25b68gd9fl25b68gd9fl25b68gd9fljlYI45p58gd9fl25b68gd9fl25b68gd9fl25b68gd9fl25b68gjK7p5b68gd9fl25b68gd9fl25b68gd3b68gd9fl25b68gdCYI45p5p5b68gd9fl25b68gd9fl25b68gd9gd9fl25b68gd9fl25b68gd9fl25b68gdCp5b68gd3b68gdCp5b68gd9fl25b68gCYI45p5b68gd9fl25b68gd9fl25b68gd9fl25b68gd9fl25b68gd9fl25b68gd9fl25b68gd9fl25b68gd9fl25b6pKYI45p5b68gd9fl25b68gd9fl25b68gd9fl25b68fl25b68fl25b68gd9fl25b68gd9fl25b68gpK7p5b68gd9fl25bYI45b68gd9fl25b68gd9fl258gpb68gd9fl25b68gd9fl25b68gd9fl25b68gd9flCYI45p58gd9fl25b68gd9fl25b68gd9fl25b68gd9fl25b68gdYI45p5b68gd9fl25b68gd9fl25b68gd9fl25b68gd9fl25b68gd9fl25b68gdCp5b68gd9fl25b68gdjYI45p5b68gd9fl25b68gd9fl25b68gd9fl25b68gd9fl25b68gd9fl25b68gd9fl25b68gd9fl25b68gdCK7p5b68gd9fl25b68gdCKYI45b68gd9fl25b68gd9fl25b68gd9fl25b68gd9fl25b68gd9fl25b68gd9fl25b68gd9fl25b68gjKYI45b68gd9fl25b68gd9fl25b68gd9fl25b68gd9fl25b68gd9fl25b68gd9fl25b68gd9fl25b68gd3b68gdCKYI45b68gd9fl25b68gd3b68gd9fl25b68gd9fl25b68gdCK7p5b68gd9fl25b68gd9fl25b68gd9fl25b68gd3CKYI45p5b68gd9fl25b68gd9fl25b68gd9fl25b68gd9fl25b68gd9gd9fl25b68gd9fl25b68gKYI45b68gd9fl25b68gd9fl25b68gd9fl25b68gd9fl25b68gd9fl25p5b68gjYI45p5b68gd9fl25b68gd9fl25b68gd9fl25b68gd9fl25b68gd9fl25b68gd3b68gdCK7p5b68gd9jKYI45p5b68gd9fl25b68gd9fl25b68gd9fl25b68gd9fl25b68gdjYI45p5b68gd9fl25b68gd9fl25b68gd9fl25b68gd9fl25b68gd9fl25b68gd9fl25b68gd9fl25b68gd9CKYI45p5p5b68gd9fl25b68fl25b68gd9fl25b68gd9fl25b68gd9fl25b68gd9fl25b68gd9fl25b68gd9fjK7p5b68gd9fl25b68gd3b68gd9gYI45p5b68fl25b68fl25b68gd9fl25b68fl25b68fl25b68gd9fl25b68gd3b68gd9fl25b68gd9fl25b68gpKYI45b68gd9fl25b68gd9fl25b68gd9fl25b68gd9fl25b68gd9fl25b68gd9fl25b68gd9fl25b68gd3K7p5b68gd9fl25b68gdCK7YI45b68gd9fl25b68gd9fl25b68gd9fl25b68gd9fl25b68gd9fl25b68gd9fl25b68gd9fl25b68gd9fl25b68gd9fl25b68gjKYI45p5b68gd9fl25b68gd9fl25b68fl25b68gd9fl25b68gd9fl25b68gd9fl25b68gd9fl25b68gdCK7p5b68gd9flCYI45p5b68gd9fl25b68gd9fl25b68gd9fl25b68gd9fl25b68gd9fl25b68gd9fl25b68gd3YI45b68gd9fl25b68gd9fl25b68gd9fl25b68gd9fl25b68fl25b68gd9fl25b68gd9fl25b68gd9fl25b68gd9fl25b6pKYI45p5b68gd9fl25b68gd9fl25b68gd9fl25b68gd9fl25b68gd3b68gdCp5b68gd9fl25b68gd3b68gdjYI45p5b68gd9fl25b68gd9fl25b68gd9fl25b68gd9fl25b68gd9fl25b68gd9fl25b68gd9fl25b68gd9fl25b68gdCK7p5b68gdCKYI45p5b68gd9fl25b68gd9fl25b68gd9fl25b68gd9fl25b68gd9fl25b68gdCYI45b68gd9fl25b68gd9fl25b68gd9fl25b68gd9fl25b68gd9fl25b68gd9fl25b68gd3b68gd9fl25b68gd9gpYI45p5b68gd9fl25b68gd9fl25b68gd9fl25b68gd9fl25b68gd9flCp5b68gd9fl25b68gd9fl25b68gd9fl25b68gdCYI45p5b68gd9fl25b68fl25b68gd9fl25b68gd9fl25b68gd9fl25b68gd3b68gd9flCYI45p58gd9fl25b68gd9fl25b68gd9fl25b68gd9fl25b68gd9fl25b68gd9fl25b68gCp5b68gd9fl25b68gdCYI45p5b68gd9fl25b68gd9fl25b68gd9fl25b68gd9fl25b68gd9fl25b68gd9fl25b68gd9gdCK7p5b68gd9fl25b68gdCKYI45p58gd9fl25b68gd9fl25b68gd9fl25b68gd9fl25b68gd9fl25CYI45p5b68gd9fl25b68gd9fl25b68gd9fl25b68gd9fl25b68gd9fl25b68gd9fl25b68gdCK7p5b68gdCKYI458gd9fl25b68gd9fl25b68gd9fl25b68gd9fl25b68gd9fl25b68gd9fl25b68gjb68gd9fl25b68gdCKYI45b68gd9fl25b68gd9fl25b68gd9fl25b68gd9fl25b68gd3b68gd9fl25b68gd3b68gd9fpYI45b68gd9fl25b68gd9fl25b68gd9fl25b68gd9fl25b68gd9fl25b68gd9fl25b68gdCYI45b68gd9fl25b68gd9fl25b68gd9fl25b68gd9fl25b68gd9fl25b68gd9fl25b68gd3b68gd9fl25b68gdCp5b68gd9fljYI45b68gd9fl25b68gd9fl25b68gd9fl25b68gd9fl25b68gd9fl25b68gd3b68gd9fl25b68gd9fl25b68gdCKYI45b68fl25b68gd9fl25b68gd9fl25b68gd9fl25b68gd9fl25b68gd9fl2YI45b68gd9fl25b68gd9fl25b68gd9fl25b68gd9fl25b68gd9fl25b68gd9fl25b68gd9Cp5b68gdCKYI45b68gd9fl25b68gd9fl25b68gd9fl25b68gd9fl25b68gd9fl25b68gd9fl25b68gd3b68fl25b68gd9fl25b68gdjYI45b68fl25b68gd9fl25b68gd9fl25bYI45b68gd9fl25b68gd9fl25b68gd3bYI45b68gd9fl25b68gd9fl25b68gd9fl25b68gd9fl25b68gd9fl25b68gd9fl25b68gd9fl25b68gdCp5b68gd9fl25b68gdjYI45b68gd9fl258gd9fl25b68gd9fl25b68gd9fl258gd9fl25b68gd9fl25b68gd3b68gdCp5b68gdCp5b68gd3b68gd3b68gd9fl25b68gjYI45b68gd9fl25b68gd9fl25b68gd9fl25b68gd9fl25b68gd9fl25b68gdCp5b68gYI45p5b68gd9fl25b68fl25b68gd9fl25b68gd9fl25b68gd9fl25b68gd9fl25b68gd3K7p5b68gd9fl25b68gd3KYI45gjb68gd9fl25b68gd9fl25b68gd9fl25b68gd9fl25b68gd9fl25b68gdjYI45p5b68gd9fl25b68gd9fl25b68gd9gd9fl25b68gd9fl25b68gd9fl25b68gd9fl25b68gd9fl25b68gd9flCKYI45p5b68gd9fl25b68gd9fl25b68gd9fl25b68gd9fl25b68gd9fl25b68gd9fl25b68gd9gd9fl25b68gdCKYI45p5b68gd9fl25b68gd9fl25b68gd9fl25b68gd9fl25b68gd9fl25b68gd9fl25b68gd3b68gd9fl25b68gd9fl25b68gpKYI45p5b68gd9fl25b68gd9fl25b68gd9fl25b68gd9fl25b68gd9fl25b68gd9fl25b68gd9fl25b68gd9fl25b68gdCK7YI45b68gd9fl25b68gd9fl25b68gd9fl25b68gdCKYI45p5b68gd9fl25b68gd9fl25b68gd9fl25b68gd9fl25b68gd9fl25b68gdCKYI45p5p5b68gd9fl25b68gd9fl25b68gd9fl25b68gd9fl25b68gd9fl25b68gd9fl25b68gd9fl25b68gdjYI45p5b68gd9fl25b68gd9fl25b68gd9fl25b68gd3b68gd9fl25b68gd9fl25b68gd9fl25b68gd9fl25b68gd3YI45b68gd9fl25b68gd9fl25b68gd9fl25b68gd9fl25b68gd9fl25b68gd9fl25b68gd9fl25b68gdCKYI45p5b68gd9fl25b68gd9fl25b68gd9fl25b68gd9fl25b68gd9fl25b68gdCK7YI45p5b68gd9fl25b68gd9fl25b68gd9fl25b68gd9fl25b68gd9fl25b68gd9fl25pYI45p5b68gd9fl25b68gd9fl25b68gd9fl25b68gd9fl25b68gd9fl25b68gd9fl25b68gd9fl25b68gd9fl25b6Cp5b68gdCKYI45b68gd9fl25b68gd9fl25b68gd9fl25b68gd9fl25b68gd9gd9fl25b68gdCp5b68gd9fl25b68gjKYI45p5b68gd9fl25b68gd9fl25b68gd9fl25b68gd9fl25b68gd9fl25b68gd9fl25b68gd9fl25b68gdCK7YI45b68fl25b68gd9fl25b68gd9fl25b68gd9fl25b68gd9fl25b68gd9fl25b68gd3b68gdCKYI45b68gd9fl25b68gd9fl25b68gd9fl25b68gd9fl25b68gdCKYI45p5b68gd9fl25b68gd9fl25b68gd9fl25b68gd9fl25b68gd9fl25b68gd9fl25b68gdCp5b68gd9fl25b68gd9fl25b68gd3b68gpKYI45p5b68gd9fl25b68gd9fl25b68gd9fl25b68gd9fl25b68gd9fl25b68gd9fl25b68gd3b68gd9fl25bYI45p5p5b68gd9fl25b68gd9fl25b68gd9fl25b68gd9fl25b68gdCKpb68gd9fl25b68gd3b68gdCKYI45b68gd9fl25b68gd9fl25b68gd9gd9fl25b68gd9fl25b68gd9fl25b68gd9fl25b68gd9fl25b68gd9fl25b68gdCK7p5b68fl25b6KYI45b68gd9fl25b68gd9fl25b68gd9fl25b68gd9fl25b68gd9fl25b68gd9fl25b68gd9fjb68gjYI45p5b68gd9fl25b68gd9fl25b68gd9fl25b68gd9fl25b68gd9fl25b68gjKYI45p5b68gd9fl25b68gd9fl25b68gd9fl25b68gd9fl25b68gd9fl25b68gd9fl25b68gdCp5b68gd9fl25b68gd3YI45p5b68gd9fl25b68gd9fl25b68gd9fl25b68gd9fl25b68gd9fl25b68gd9fl25b68gd9fl25b68gd3b68gd9fl25b68gdCKYI45b68gd9fl25b68gd9fl25b68gd9fl25b68gd9gd9fl25b68gd9fl25b68gd9fl25b68gd9fl25b68gd9gdCK7p5b68gd9fl25b68gd9fl25b68gd9fl25b68gdjYI45b68gd9fl25b68gd9fl25b68gd9fl25b68gd9fl25b68gd9fl25b68gd9fl25b68gd3b68gYI45b68gd9fl25b68gd9fl25b68gd9fl25b68gd9gd3b68gd9fl25b68gd9fl25b68gd9fl25b68gd9fl25b68gjKYI45p5p5b68gd9fl25b68gd9fl25b68gd9fl25b68gd9fl25b68gd9fl25b68gd9gd9fl25b68gjKYI45p5b68gd9fl25b68gd9fl25b68gd9fl25b68gd9fl25b68gd9fl25b68gd9fl25b68gd9fl25b68gd9fl25b68gjYI45p5b68gd9fl25b68gd9fl25b68gd9fl25b68gd9fl25b68gd9fl25b68gd9gd3b68gd3KYI45b68gd9fl25b68gd9fl25b68gd9fl25b68gd9fl25b68gd9fl25b68gd9fl25b68gdCK7p5b68gdCp5Cp5b68gjKYI45b68gd9fl25b68gd9fl25b68gd9fl25b68gd9fl25b68gd9fl25b68gdCp5b68gd9fl25b68gd9fl25b68gd3b68gjYI45b68gd9fl25b68gd9fl25b68gd9fl25b68gd9fl25b68gdYI45gjb68gd9fl25b68gd9fl25b68gd9fl25b68gd9fl25b68gd9fl25b68gdCKYI45p5b68gd9fl25b68gd9fl25b68gd9fl25b68gd9fl25b68gjb68gd9fl25b68gd3b68gd9fl25b68gdCK7p5b68gd9fl25b68gjb68gd3YI45b68gd9fl25b68gd9fl25b68gd9fl25b68gd9fl25b68gd9fl25bYI45p5b68gd9fl25b68gd9fl25b68gd9fl25b68gd9fl25b68gd9fl25b68gd3b68gd9fl25b68gd3YI45b68gd9fl25b68gd9fl25b68gd9fl25b68gd9fl25b68gd3b68gd3YI45b68gd9fl25b68gd9fl25b68gd9fl25b68gd9gd9fl25b68gd9fl25b68gd9fl25b68gd3b68gjb68gd9fl25b68gd9fl25b68gdCp5b68gdj'
unique_labels_pre = sorted(list(set(bout_screen)))

bout_screen = re.sub('9', 'h', bout_screen)
bout_screen = re.sub('8', 'H', bout_screen)

bout_screen = re.sub('7', 'C', bout_screen)

bout_screen = re.sub('6', 'J', bout_screen)
bout_screen = re.sub('5', 'j', bout_screen)

bout_screen = re.sub('4', 'p', bout_screen)
bout_screen = re.sub('3', 'P', bout_screen)
bout_screen = re.sub('2', 'A', bout_screen)


# labels_post_screen = sorted(list(set(bout_screen)))
labels_post_screen = ['Y', 'I', 'p', 'j', 'b', 'J', 'H', 'g', 'd', 'h', 'f', 'l', 'A', 'P', 'C', 'K']
tm_screen, _ = hf.get_transition_matrix(bout_screen, labels_post_screen)
# ---------------------------------------------------------------------------------------------------------------
k = np.where(sum(tm_screen) / sum(sum(tm_screen)) * 100 <= 0.1)
tmd = np.delete(tm_screen, k, axis=1)
tmd = np.delete(tmd, k, axis=0)
print(np.delete(labels_post_screen, k))

node_labels = np.delete(labels_post_screen, k)
node_labels = ['Y', 'i+', 'p1', 'j1', 'b', 'j2', 'he1', 'g', 'd', 'he2', 'f', 'l', 'p2', 'p3', 'c+', 'k+']
tmpd = (tmd.T / np.sum(tmd, axis=1)).T
tmpd = hf.get_node_matrix(tmpd, 2)
node_size = np.round(np.sum(tmd, axis=1) / np.min(np.sum(tmd, axis=1)), 2) * 200
'Plot Transition Matrix and Transition Diagram'
# pf.plot_transition_matrix(tmpd,
#                           node_labels,
#                           'C:/Users/User/Documents/Doktorarbeit/data/bu01ye01/py_data/' + 'screen' + '_matrix.pdf',
#                           'screen')
pf.plot_transition_diagram(tmpd,
                           node_labels,
                           node_size,
                           10,
                           'C:/Users/User/Documents/Doktorarbeit/data/bu01ye01/py_data/' + 'screen' + '_graph.pdf',
                           'screen')
plt.show()
