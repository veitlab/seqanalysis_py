import re
import sys
import glob
import numpy as np
import seqanalysis.util.helper_functions as hf
import matplotlib.pyplot as plt
from IPython import embed

# --------------------------------------------------------------------------------------------------------------
syllables = ["g+", "g+Og*", "g+Mg*", "g+Ng*"]
playback_syl = ["O", "N", "M"]
path = "C:/Users/veitlab/Documents/pk02gr02/"
folders = ["exp1_target_g_with_gci", "baseline/male"]

colours = ["blue", "red", "black"]
# --------------------------------------------------------------------------------------------------------------

folder_lengths = []
for folder_idx in range(len(folders)):
    real_path = path + folders[folder_idx] + "/*/*cbin.not.mat"
    file_list = glob.glob(real_path)

    seqs = hf.get_labels(file_list, ["E", "S"])
    bouts, _ = hf.get_bouts(seqs, "g")

    fig1, ax1 = plt.subplots(1, 1, figsize=(21 / 2.54, 12 / 2.54))
    fig1.subplots_adjust(left=0.13, bottom=0.17, top=0.99, right=0.99)
    for PB_syl_idx in range(len(playback_syl)):
        [unique_g_PB_g, counts_g_PB_g] = np.unique(
            re.findall("Pg{,10}", re.sub("g+" + playback_syl[PB_syl_idx], "P", bouts)),
            return_counts=True,
        )
        x_axis = [len(i) for i in unique_g_PB_g]
        ax1.plot(
            x_axis,
            counts_g_PB_g / np.sum(counts_g_PB_g),
            "o-",
            color=colours[PB_syl_idx],
            label=playback_syl[PB_syl_idx] + "g",
        )
        ax1.set_xticks(x_axis, labels=unique_g_PB_g)
        ax1.spines["top"].set_visible(False)
        ax1.spines["right"].set_visible(False)
        ax1.tick_params(width=2)
        ax1.tick_params(axis="both", which="major", labelsize=15)
        ax1.set_xlabel("Repeat number", fontsize=20)
        ax1.set_ylabel("Rel. frequency", fontsize=20)
        plt.rcParams["font.size"] = 20
        plt.rcParams["axes.linewidth"] = 2

    plt.legend(loc=1)
    plt.show()

    gPBg_list = []
    for find_syl_idx in range(len(syllables)):
        g_playback_g = re.findall(syllables[find_syl_idx], bouts)
        gPBg_list.append(g_playback_g)

        embed()
        quit()
