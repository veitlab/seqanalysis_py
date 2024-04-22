import re
import sys
import glob
import numpy as np
import seaborn as sns
import jacquiutil.helper_functions as hf
import matplotlib.pyplot as plt
from IPython import embed


def get_data(path, intro_notes, bout_chunk):
    file_list = glob.glob(path)

    seqs = hf.get_labels(file_list, intro_notes)
    bouts, _ = hf.get_bouts(seqs, bout_chunk)

    return bouts


def get_catch_data(path, intro_notes, bout_chunk):
    # assembles the files in the path
    file_list = []
    list = glob.glob(path)

    for i in range(len(list)):
        with open(list[i] + 'batch.catch', 'r') as file:
            line_list = file.readlines()
            file_list.extend([list[i] + item.rstrip() + '.not.mat' for item in line_list])
    # embed()
    # quit()
    seqs = hf.get_labels(file_list, intro_notes)
    bouts, _ = hf.get_bouts(seqs, bout_chunk)

    return bouts



if __name__ == '__main__':

    # ------------------------------------------------------------------------------------------------------------
    # or05pk04
    plt.rcParams['font.size'] = 15
    plt.rcParams['axes.linewidth'] = 2

    bouts = get_data('D:/Birds/or05pk04/Gulnur/exp*/*/*cbin.not.mat', ['E', 'S', 'k', 'f'], 'b')

    matrix = np.zeros([5, 11])
    matrix[:] = np.nan
    for repeat_counter in range(5):
        search_regex = '(?<=[^bRT])b{' + str(repeat_counter + 1) + '}[RT]b*(?<![^bRT])'
        b_PB_b = re.findall(search_regex, bouts)
        only_b = [sub[repeat_counter + 2:] for sub in b_PB_b]
        letters, numbers = np.unique(only_b, return_counts=True)
        print(np.unique(only_b, return_counts=True))
        lengths = [len(item) for item in letters]
        print(lengths)
        matrix[repeat_counter, lengths] = numbers

    matrix_perc = matrix.T / np.nansum(matrix, axis=1).T
    matrix_perc = matrix_perc.T
    yticklbl = ['b{1}RT', 'b{2}RT', 'b{3}RT', 'b{4}RT', 'b{5}RT']
    fig, ax = plt.subplots(1, 1, figsize=(25 / 2.54, 14 / 2.54))
    fig.subplots_adjust(left=0.15, bottom=0.17, top=0.9, right=0.99, wspace=0.5, hspace=0.7)
    hm = sns.heatmap(matrix_perc, ax=ax, annot=True, vmin=0, vmax=1, fmt='.2f', cmap='Greys',
                     xticklabels=range(11), yticklabels=yticklbl)
    ax.set_yticklabels(hm.get_yticklabels(), rotation=0)
    ax.tick_params(left=False, bottom=False)
    sns.despine(top=False, right=False, left=False, bottom=False)
    ax.set_title('playback matrix or05pk04')
    ax.set_xlabel('number of b')

    plt.show()
    embed()
    quit()
    # fig.savefig('C:/Users/veitlab/Documents/or05pk04/figures/playback_matrix.pdf', dpi=300)

    # ------------------------------------------------------------------------------------------------------------
