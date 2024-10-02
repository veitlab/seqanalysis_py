import re
import glob
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from seqanalysis.util.get_data_transition_diagram import get_labels, get_bouts


def get_catch_data(path, intro_notes, bout_chunk):
    list = glob.glob(path)

    file_list = []
    for i in range(len(list)):
        with open(list[i], 'r') as file:
            line_list = file.readlines()
            file_list.extend([list[i].rstrip('batch.notcatch') + item.rstrip() + '.not.mat' for item in line_list])

    seqs = get_labels(file_list, intro_notes)
    bouts, _ = get_bouts(seqs, bout_chunk)

    return bouts


def plot_sensitivity_bar_plot(syllable, always_syl, root_dir, folders, days, condition):

    # Apply a fancy colormap to the figure with the numbers of colors equal to the possible syllables
    cmap_jet = plt.get_cmap('jet')
    colors = cmap_jet(np.linspace(0, 1, len(condition)))

    sensitivity_matrix = np.full([len(syllable), len(folders)], np.nan)

    for folder_idx, folder in enumerate(folders):
        bouts = get_catch_data(root_dir + folder, ['_'], always_syl)
        for syl_idx, syl in enumerate(syllable):
            all_syl = re.findall(syl, bouts)
            sensitivity_matrix[syl_idx, folder_idx] = len(all_syl)

    sensitivity_matrix_percent = np.round(sensitivity_matrix / np.sum(sensitivity_matrix, axis=0) * 100, 2)
    df = pd.DataFrame(sensitivity_matrix_percent.T, columns=condition)
    df.insert(0, 'days', days, True)
    fig, ax = plt.subplots(1, 1, figsize=(21 / 2.54, 12 / 2.54))
    fig.subplots_adjust(left=0.13, bottom=0.18, top=0.99, right=0.94, wspace=0.5, hspace=0.7)

    df.plot(x='days',
            kind='bar',
            stacked=True,
            color=colors,
            ax=ax)

    ax.set_xlabel('training day')
    ax.set_ylabel('probability [%]')
    plt.legend()

    return fig


if __name__ == "__main__":

    file_path_glob = 'C:/Users/User/Documents/Doktorarbeit/data/bu02bk02/'
    file_path_save = 'C:/Users/User/Documents/Doktorarbeit/project2_repeats/bu02bk02/figures/'

    folders = ['/exp1_b4+_WN/240507/batch.notcatch',
               '/exp1_b4+_WN/240508/batch.notcatch',
               '/exp1_b4+_WN/240509/batch.notcatch',
               '/exp1_b4+_WN/240510/batch.notcatch']

    x_axis_labels = ['T1', 'T2', 'T3', 'T4']

    # this syllable is in every song file you want to analyse (is from the do_transition_matrix: bout_chunk)
    always_syl = 'd'
    # use regular expressions to find the branch point
    syllable = ['dccb{,3}w+', 'dccb{4,}(?!w)']

    line_labels = ['WN', '*']

    fig2 = plot_sensitivity_bar_plot(syllable, always_syl, file_path_glob, folders, x_axis_labels, line_labels)

    fig2.savefig(f'{file_path_save}sensitivity_bar_plot.svg')
    plt.show()