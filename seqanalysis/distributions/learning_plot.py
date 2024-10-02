import re
import glob
import pathlib
import numpy as np
import matplotlib.pyplot as plt
from seqanalysis.util.get_data_transition_diagram import get_labels, get_bouts, get_analyse_files_data
from IPython import embed


def get_catch_data(folder_path, analyse_files, intro_notes, bout_chunk):

    file_list = get_analyse_files_data(folder_path, analyse_files)
    seqs = get_labels(file_list, intro_notes, 'i')
    bouts, _ = get_bouts(seqs, bout_chunk)

    return bouts


def main(branch_syll, target_syll, always_syll, root_dir, folders, analyse_file, x_labels):
    # Apply a fancy colormap to the figure with the numbers of colors equal to the possible syllables
    cmap_jet = plt.get_cmap('jet')
    colors = cmap_jet(np.linspace(0, 1, len(branch_syll)))

    branching_from_syl_matrix = np.zeros([len(folders), len(branch_syll)])
    for folder_idx in range(len(folders)):
        folder_path = pathlib.Path(root_dir + folders[folder_idx])
        bouts = get_catch_data(folder_path, analyse_file, ['_'], always_syll)

        for trans_syl_idx in range(len(branch_syll)):
            exp_occur = re.findall(target_syll + branch_syll[trans_syl_idx], bouts)
            branching_from_syl_matrix[folder_idx, trans_syl_idx] = len(exp_occur)

    matrix_prob = (branching_from_syl_matrix.T / np.sum(branching_from_syl_matrix, axis=1)).T * 100

    fig1, ax1 = plt.subplots(1, 1, figsize=(21 / 2.54, 12 / 2.54))
    fig1.subplots_adjust(left=0.13, bottom=0.18, top=0.99, right=0.94, wspace=0.5, hspace=0.7)

    for plot_idx in range(len(branch_syll)):
        ax1.plot(range(len(folders)), matrix_prob[:, plot_idx],
                 label=target_syll + branch_syll[plot_idx], marker='o', linewidth=4, color=colors[plot_idx])

    ax1.set_xticks(range(len(folders)), x_labels)
    ax1.set_ylabel('probability [%]')

    plt.legend()

    return fig1


if __name__ == "__main__":
    file_path_glob = 'C:/Users/User/Documents/Doktorarbeit/data/'
    file_path_save = 'C:/Users/User/Documents/Doktorarbeit/data/Franzi/figures/'

    data_folders = ['Franzi/']

    x_axis_labels = ['pre']

    # this syllable is in every song file you want to analyse (is from the do_transition_matrix: bout_chunk)
    always_syl = '4'
    # use regular expressions to find the branch point
    target_syl = '23'
    branch_syl = ['2', '(?!2).']

    fig = main(branch_syl, target_syl, always_syl, file_path_glob, data_folders, 'keep', x_axis_labels)

    fig.savefig(f'{file_path_save}trans_prob_from_[ef]b.svg')
    fig.savefig(f'{file_path_save}trans_prob_from_[ef]b.png')

    plt.show()
