import re
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from seqanalysis.util.helper_functions import get_catch_data


def plot_sensitivity_bar_plot(syllable, root_dir, folders, days, condition, colors):
    if len(colors) != len(condition) or colors is None:
        # Apply a fancy colormap to the figure
        cmap_jet = plt.get_cmap("jet")
        colors = cmap_jet(np.linspace(0, 1, len(condition)))
    plt.rcParams.update({"font.size": 15})
    sensitivity_matrix = np.full([len(syllable), len(folders)], np.nan)
    intro_replacement = ["x", "j", "a", "d", "c", "e"]
    for folder_idx, folder in enumerate(folders):
        bouts = get_catch_data(root_dir + folder, ["_"], "d", intro_replacement)
        for syl_idx, syl in enumerate(syllable):
            all_syl = re.findall(syl, bouts)
            sensitivity_matrix[syl_idx, folder_idx] = len(all_syl)

    sensitivity_matrix_percent = np.round(
        sensitivity_matrix / np.sum(sensitivity_matrix, axis=0) * 100, 2
    )
    df = pd.DataFrame(sensitivity_matrix_percent.T, columns=condition)
    df.insert(0, "days", days, True)
    fig, ax = plt.subplots(1, 1, figsize=(21 / 2.54, 12 / 2.54))
    fig.subplots_adjust(
        left=0.13, bottom=0.18, top=0.99, right=0.74, wspace=0.5, hspace=0.7
    )

    df.plot(x="days", kind="bar", stacked=True, color=colors, ax=ax)

    ax.set_xlabel("training day", fontsize=20)
    ax.set_ylabel("probability [%]", fontsize=20)
    plt.legend(bbox_to_anchor=(1.05, 0.5), fontsize=20)
    fig.savefig(
        "/home/alexander/projects/masterthesis/figures/results/sensitivity_plot.pdf"
    )
    plt.show()


if __name__ == "__main__":
    syllable = ["bnmml{,2}w+", "bnmml{,2}(?!w+)"]
    file_path_glob = "/home/alexander/projects/data/gr88gr06/"

    folders = [
        "training/240529/batch.notcatch",
        "training/240530/batch.notcatch",
        "training/240531/batch.notcatch",
        "training/240601/batch.notcatch",
        "training/240602/batch.notcatch",
    ]

    labels = ["T1", "T2", "T3", "T4", "T5"]
    plot_label = ["WN", "*"]

    colors_2 = ["#35787f", "#b4e57b"]
    plot_sensitivity_bar_plot(
        syllable, file_path_glob, folders, labels, plot_label, colors_2
    )
