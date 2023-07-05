import re
import sys
import glob
import yaml
import numpy as np
import plot_functions as pf
import helper_functions as hf
import matplotlib.pyplot as plt
from IPython import embed

def get_data(path, intro_notes, bout_chunk):
    file_list = glob.glob(path)

    seqs = hf.get_labels(file_list, intro_notes)
    bouts, _ = hf.get_bouts(seqs, bout_chunk)

    return bouts


def main(folder_path, folder, syl, labels):

    folder_lengths = []
    for folder_idx in range(len(folder)):
        bouts = get_data(folder_path+folder[folder_idx]+'/*cbin.not.mat', ['E', 'S'], 'f')
        all_b = re.findall(syl, bouts)
        lengths = [len(item)-1 for item in all_b]
        folder_lengths.append(lengths)

    max_steps = np.max(np.hstack(folder_lengths))
    steps = np.linspace(0,max_steps+1,max_steps+2)

    fig1, ax1 = plt.subplots(1, 1, figsize=(21 / 2.54, 12 / 2.54), sharex=True)
    fig1.subplots_adjust(left=0.13, bottom=0.17, top=0.99, right=0.99, wspace=0.5, hspace=0.7)

    for plot_idx in range(len(folder_lengths)):
        counts, bins = np.histogram(folder_lengths[plot_idx], steps)

        ax1.plot(bins[:-1], counts/np.sum(counts), label=labels[plot_idx], linewidth=4)

    ax1.set_xlabel('Repeat number', fontsize=25)
    ax1.set_ylabel('Rel. frequency', fontsize=25)
    plt.rcParams['font.size'] = 25
    plt.rcParams['axes.linewidth'] = 2
    ax1.spines["top"].set_visible(False)
    ax1.spines["right"].set_visible(False)
    ax1.tick_params(width=2)
    ax1.tick_params(axis='both', which='major', labelsize=20)
    plt.legend()
    plt.show()

    embed()
    quit()


if __name__ == '__main__':
    # this script plots transition matrix and diagrams
    #
    # INPUT:
    # sys.argv[1] = yaml file of the bird, example: example_yaml.yaml
    # sys.argv[2] = analysis catch or all files: input: catch, all
    #
    # OUTPUT:
    # figures
    syllable = 'db+'
    path = 'D:/GP2022/data/or05pk04/'
    labels = ['Before training', 'After training']

    folders = ['label_or/221114_bp', 'label_or/221121_lastbase']

    main(path, folders, syllable, labels)