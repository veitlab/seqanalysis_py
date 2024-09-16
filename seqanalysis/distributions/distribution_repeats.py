import re
import sys
import glob
import yaml
import pathlib
import numpy as np
import seqanalysis.util.plot_functions as pf
import seqanalysis.util.helper_functions as hf
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
        with open(list[i] + "batch.catch", "r") as file:
            line_list = file.readlines()
            file_list.extend(
                [list[i] + item.rstrip() + ".not.mat" for item in line_list]
            )
    seqs = hf.get_labels(file_list, intro_notes)
    bouts, _ = hf.get_bouts(seqs, bout_chunk)

    return bouts


def main(folder_path, syl, labels, bout_chunk):
    folder_path = pathlib.Path(folder_path)
    files = list(folder_path.glob("**/[!syll*][!check]*.not.mat"))

    seqs = hf.get_labels(
        files,
        [
            "_",
            "j",
            "a",
            "d",
            "c",
            "e",
        ],
    )
    bouts, _ = hf.get_bouts(seqs, bout_chunk)

    # bouts = get_catch_data(folder_path+folder[folder_idx], ['E', 'S', 'k', 'f'], bout_chunk)
    all_b = re.findall(syl, bouts)
    folder_lengths = [len(item) - 0 for item in all_b]

    max_steps = np.max(folder_lengths)
    steps = np.linspace(0, max_steps + 1, max_steps + 2)

    fig1, ax1 = plt.subplots(1, 1, figsize=(21 / 2.54, 12 / 2.54), sharex=True)
    fig1.subplots_adjust(
        left=0.13, bottom=0.17, top=0.99, right=0.99, wspace=0.5, hspace=0.7
    )

    counts, bins = np.histogram(folder_lengths, steps)

    ax1.plot(bins[:-1], counts / np.sum(counts), label=labels, linewidth=4)

    ax1.set_xlabel("Repeat number", fontsize=25)
    ax1.set_ylabel("Rel. frequency", fontsize=25)
    plt.rcParams["font.size"] = 25
    plt.rcParams["axes.linewidth"] = 2
    ax1.spines["top"].set_visible(False)
    ax1.spines["right"].set_visible(False)
    ax1.tick_params(width=2)
    ax1.tick_params(axis="both", which="major", labelsize=20)
    # ax1.plot(bins[1:-1], counts[1:]/np.sum(counts[1:]), label=labels[plot_idx], linewidth=4)

    plt.legend()
    fig1.savefig(f"{syl}_repeats.pdf")

    plt.show()


if __name__ == "__main__":
    # this script plots transition matrix and diagrams
    #
    # INPUT:
    # sys.argv[1] = yaml file of the bird, example: example_yaml.yaml
    # sys.argv[2] = analysis catch or all files: input: catch, all
    #
    # OUTPUT:
    # figures
    syllable = "l+"
    path = "/home/acfw/Projects/data/gr88gr06_tweetynet/training_data"
    labels = [f"{syllable}"]
    bout_chunk = "ll"

    plt.rcParams["svg.fonttype"] = (
        "none"  # this is so that svg figures save text as text and not the single letters
    )
    plt.rcParams["font.size"] = 20
    plt.rcParams["axes.linewidth"] = 2

    main(path, syllable, labels, bout_chunk)
