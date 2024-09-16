import re
import sys
import glob
import yaml
import numpy as np
import seqanalysis.util.plot_functions as pf
import seqanalysis.util.helper_functions as hf
import matplotlib.pyplot as plt
from IPython import embed


def get_data(path, intro_notes, bout_chunk):
    """Retrieve bout data from files in the specified path."""
    file_list = glob.glob(path)
    seqs = hf.get_labels(file_list, intro_notes)
    bouts, _ = hf.get_bouts(seqs, bout_chunk)
    return bouts


def calc_transition_matrix(folder_path, folder, target_syl, trans_syl, transition_func):
    """Calculate the transition matrix based on the specified transition function."""
    transition_matrix = np.zeros([len(folder), len(trans_syl)])

    for folder_idx, folder_name in enumerate(folder):
        bouts = get_data(
            folder_path + folder_name + "*cbin.not.mat", ["E", "S"], target_syl
        )

        for trans_syl_idx, trans_syl_label in enumerate(trans_syl):
            exp_occur = transition_func(bouts, target_syl, trans_syl_label)
            transition_matrix[folder_idx, trans_syl_idx] = len(exp_occur)

    return transition_matrix


def calc_convergence_to_syl(bouts, target_syl, trans_syl):
    """Calculate occurrences of convergence to a syllable in bouts."""
    return re.findall(trans_syl + target_syl, bouts)


def calc_branching_from_syl(bouts, target_syl, trans_syl):
    """Calculate occurrences of branching from a syllable in bouts."""
    return re.findall(target_syl + trans_syl, bouts)


def main(folder_path, folder, target_syl, trans_syl, labels, transition):
    """Main function to analyze and plot transition matrix."""
    if transition == "converg":
        transition_func = calc_convergence_to_syl
    elif transition == "branch":
        transition_func = calc_branching_from_syl
    else:
        raise ValueError("Invalid transition type. Use 'converg' or 'branch'.")

    matrix = calc_transition_matrix(
        folder_path, folder, target_syl, trans_syl, transition_func
    )
    matrix_prob = (matrix.T / np.sum(matrix, axis=1)).T

    embed()
    # Plotting the transition matrix
    fig1, ax1 = plt.subplots(1, 1, figsize=(21 / 2.54, 12 / 2.54), sharex=True)
    fig1.subplots_adjust(
        left=0.13, bottom=0.17, top=0.99, right=0.99, wspace=0.5, hspace=0.7
    )

    for trans_syl_idx, trans_syl_label in enumerate(trans_syl):
        label = (
            f"{trans_syl_label}{target_syl}"
            if transition == "converg"
            else f"{target_syl}{trans_syl_label}"
        )
        ax1.plot(
            range(len(folder)), matrix_prob[:, trans_syl_idx], label=label, linewidth=4
        )

    ax1.set_xticks(range(len(folder)), labels)
    ax1.set_ylabel("Rel. frequency", fontsize=20)

    ax1.spines["top"].set_visible(False)
    ax1.spines["right"].set_visible(False)
    ax1.tick_params(width=2)
    ax1.tick_params(axis="both", which="major", labelsize=18)
    plt.legend()
    fig1.savefig(f"..\\..\\..\\data\\or05pk04\\figures\\{transition}.svg")

    plt.show()


if __name__ == "__main__":
    # Example usage
    plt.rcParams["svg.fonttype"] = (
        "none"  # this is so that svg figures save text as text and not the single letters
    )
    plt.rcParams["font.size"] = 20
    plt.rcParams["axes.linewidth"] = 2

    # transition_type = 'branch'  # use branch or converg
    transition_type = "converg"  # use branch or converg
    target_syllables = [["f", "d", "a", "c"], ["b"]]
    # target_syllables = [['b', 'e', 'f', 'E'], ['d']]

    path = "D:/Birds/or05pk04/*/*/"
    folders = [
        "22110*/",
        "221114/",
        "221115/",
        "221116/",
        "221117/",
        "221118/",
        "221121/",
    ]
    labels = ["BL", "T1", "T2", "T3", "T4", "T5", "BL"]

    main(
        path,
        folders,
        target_syllables[1][0],
        target_syllables[0],
        labels,
        transition_type,
    )
