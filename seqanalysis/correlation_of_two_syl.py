import re
import sys
import glob
import yaml
import numpy as np
import jacquiutil.plot_functions as pf
import jacquiutil.helper_functions as hf
import matplotlib.pyplot as plt
from IPython import embed
from scipy.stats import pearsonr, linregress


def get_data(path, intro_notes, bout_chunk):
    """Retrieve bout data from files in the specified path."""
    file_list = glob.glob(path)
    seqs = hf.get_labels(file_list, intro_notes)
    bouts, _ = hf.get_bouts(seqs, bout_chunk)
    return bouts


def main(folder_path, folder, syl1, syl2, labels):


    for folder_idx, folder_name in enumerate(folder):
        bouts = get_data(folder_path + folder_name + '*cbin.not.mat', ['E', 'S'], 'b')

        embed()
        quit()


        all_syl1 = re.findall(syl1, bouts)
        all_syl2 = re.findall(syl2, bouts)
        all_syl1_lengths = [len(item)-1 for item in all_syl1]
        all_syl2_lengths = [len(item)-1 for item in all_syl2]

        fig1, ax1 = plt.subplots(1, 1, figsize=(21 / 2.54, 12 / 2.54), sharex=True)
        fig1.subplots_adjust(left=0.13, bottom=0.17, top=0.99, right=0.99, wspace=0.5, hspace=0.7)

        ax1.scatter(all_syl1_lengths, all_syl2_lengths)

        # Perform linear regression
        slope, intercept, r_value, p_value, std_err = linregress(all_syl1_lengths, all_syl2_lengths)
        # Calculate the regression line
        regression_line = slope * range(7) + intercept
        ax1.plot(range(7), regression_line)

        ax1.set_ylabel('Rel. frequency', fontsize=20)
        ax1.spines["top"].set_visible(False)
        ax1.spines["right"].set_visible(False)
        ax1.tick_params(width=2)
        ax1.tick_params(axis='both', which='major', labelsize=18)
        plt.show()


if __name__ == '__main__':
    # Example usage
    plt.rcParams['svg.fonttype'] = 'none'  # this is so that svg figures save text as text and not the single letters
    plt.rcParams['font.size'] = 20
    plt.rcParams['axes.linewidth'] = 2

    syl2 = 'db+'
    syl1 = 'd+b'

    path = 'D:/Birds/or05pk04/*/*/'
    folders = ['22110*/', '221114/', '221115/', '221116/', '221117/', '221118/', '221121/']
    labels = ['BL', 'T1', 'T2', 'T3', 'T4', 'T5', 'BL']

    main(path, folders, syl1, syl2, labels)
