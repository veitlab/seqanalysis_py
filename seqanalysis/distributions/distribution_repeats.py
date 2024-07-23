import re
import glob
import pathlib
import numpy as np
import seqanalysis.util.helper_functions as hf
import matplotlib.pyplot as plt
from IPython import embed


def get_data(path, intro_notes, bout_chunk):
    file_list = glob.glob(path)

    seqs = hf.get_labels(file_list, intro_notes)
    bouts, _ = hf.get_bouts(seqs, bout_chunk)

    return bouts


def get_catch_data(path, intro_notes, bout_chunk, intro_replacement):
    # assembles the files in the path
    file_list = []
    batch_catch_files = list(path.glob("**/[!syll*][!check]*.catch"))

    for i in range(len(batch_catch_files)):
        with open(batch_catch_files[i], "r") as file:
            line_list = file.readlines()
            file_list.extend(
                [
                    str(batch_catch_files[i].parent) + "/" + item.rstrip() + ".not.mat"
                    for item in line_list
                ]
            )
    seqs = hf.get_labels(file_list, intro_notes, intro_replacement)
    bouts, _ = hf.get_bouts(seqs, bout_chunk)

    return bouts


def mult_days(
    folder_path, syl, labels, bout_chunk, intro_recplacement, path_non_training
):
    folder_path = pathlib.Path(folder_path)
    days = [day for day in folder_path.iterdir() if day.is_dir()]
    path_non_training = pathlib.Path(path_non_training)
    days = sorted(days)
    days.append(path_non_training)
    fig1, ax1 = plt.subplots(1, 1, figsize=(8, 5), sharex=True)
    color_training_days = ["#220901", "#621708", "#941b0c", "#bc3908", "#f6aa1c"]
    # color_training_days = ["#0b132b", "#1c2541", "#3a506b", "#5bc0be", "#6fffe9"]
    for i, day in enumerate(days):
        files = list(day.glob("**/[!syll*][!check]*.not.mat"))
        intro_notes = [
            "x",
            "j",
            "a",
            "d",
            "c",
            "e",
        ]
        seqs = hf.get_labels(
            files,
            intro_notes,
            intro_recplacement,
        )
        if day.name == "first_screening":
            bouts, _ = hf.get_bouts(seqs, bout_chunk)
        else:
            bouts = get_catch_data(day, intro_notes, bout_chunk, intro_recplacement)
        all_b = re.findall(syl, bouts)
        folder_lengths = [len(item) - 0 for item in all_b]

        max_steps = np.max(folder_lengths)
        steps = np.linspace(0, max_steps + 1, max_steps + 2)

        fig1.subplots_adjust(
            left=0.13, bottom=0.17, top=0.87, right=0.99, wspace=0.5, hspace=0.7
        )

        counts, bins = np.histogram(folder_lengths, steps)
        if day.name == "first_screening":
            name = "baseline"
        else:
            name = day.name

        ax1.plot(
            bins[:-1],
            counts / np.sum(counts),
            label=name,
            linewidth=4,
            color=color_training_days[i] if day.name != "first_screening" else "red",
        )

    ax1.set_title("Repeats of syllable l+")
    ax1.set_xlabel("Repeat number", fontsize=25)
    ax1.set_ylabel("Rel. frequency", fontsize=25)
    plt.rcParams["font.size"] = 20
    plt.rcParams["axes.linewidth"] = 2
    ax1.spines["top"].set_visible(False)
    ax1.spines["right"].set_visible(False)
    ax1.tick_params(width=2)
    ax1.tick_params(axis="both", which="major", labelsize=20)
    # ax1.plot(bins[1:-1], counts[1:]/np.sum(counts[1:]), label=day, linewidth=4)

    plt.legend(fontsize=20)
    fig1.savefig(
        f"/home/alexander/projects/masterthesis/figures/results/{syl}_repeats_days.pdf"
    )

    plt.show()


def main(folder_path, syl, labels, bout_chunk, intro_recplacement):
    folder_path = pathlib.Path(folder_path)
    files = list(folder_path.glob("**/[!syll*][!check]*.not.mat"))

    seqs = hf.get_labels(
        files,
        [
            "x",
            "j",
            "a",
            "d",
            "c",
            "e",
        ],
        intro_recplacement,
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
    path = "/home/alexander/projects/data/gr88gr06/training/"
    labels = [f"{syllable}"]
    bout_chunk = "ll"
    intro_recplacement = "x"
    path_non_training = "/home/alexander/projects/data/gr88gr06/first_screening/"

    plt.rcParams["svg.fonttype"] = (
        "none"  # this is so that svg figures save text as text and not the single letters
    )
    plt.rcParams["font.size"] = 20
    plt.rcParams["axes.linewidth"] = 2
    mult_days(path, syllable, labels, bout_chunk, intro_recplacement, path_non_training)
    # main(path, syllable, labels, bout_chunk, intro_recplacement)
