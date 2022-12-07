import os
import glob
import numpy as np
import matplotlib.pyplot as plt
from IPython import embed
import pandas as pd
import scipy.io as sio
import networkx as nx
import helper_functions as hf
import plot_functions as pf




folder_path = 'D:/Birds/wh45wh40'
edge_threshold = 5 # edges < 5 % aren't shown
# move_path = 'C:/Users/User/Documents/Doktorarbeit/data/wh09pk39'

# folder = os.listdir(folder_path)
folder = ['screen']
os.chdir(folder_path)

all_bouts = []
label_trans = []
all_syldurs = []
all_gapdurs = []
for folderidx in range(len(folder)):
    os.chdir(folder[folderidx])
    days = [d for d in os.listdir(os.getcwd()) if os.path.isdir(d)]

    for dayidx in range(len(days)):
        os.chdir(days[dayidx])

        mat_list = glob.glob('*cbin.not.mat')
        durs = []
        for matidx in mat_list:
            mat = sio.loadmat(matidx)


            dur = mat['onsets'][1:] - mat['offsets'][:-1]
            all_gapdurs.append(dur)

            labels = mat['labels'][0]
            for i in range(len(labels) - 1):
                label_trans.append([labels[i] + labels[i + 1]])

        os.chdir('..')


    durs = np.vstack(np.array(all_gapdurs))
    lt = np.array(label_trans)

    unique_lt = np.unique(lt)
    for i in range(len(unique_lt)):
        x = durs[lt == unique_lt[i]]
        hist, bins = np.histogram(x, np.linspace(0, 200, 81))
        plt.plot(bins[:-1]+2.5, hist)

    plt.show()

    embed()
    quit()

