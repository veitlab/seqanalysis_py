import re
import sys
import glob
import yaml
import numpy as np
import matplotlib.pyplot as plt

import seqanalysis.util.plot_functions as pf
import seqanalysis.util.helper_functions as hf

from IPython import embed

plt.rcParams["svg.fonttype"] = (
    "none"  # this is so that svg figures save text as text and not the single letters
)
plt.rcParams["font.size"] = 18
plt.rcParams["axes.linewidth"] = 2
plt.rcParams["xtick.major.width"] = 2
plt.rcParams["ytick.major.width"] = 2
plt.rcParams["axes.spines.right"] = False
plt.rcParams["axes.spines.top"] = False
colors_rep = ["#cd001a", "#ef6a00", "#f2cd00", "#79c300", "#1961ae", "#61007d"]
colors_9 = [
    "#ff595e",
    "#ff924c",
    "#ffca3a",
    "#c5ca30",
    "#8ac926",
    "#52a675",
    "#1982c4",
    "#4267ac",
    "#6a4c93",
]
colors_2 = ["#e9d758", "#297373"]


file_path = "Z:/Gulnurs_Project/pk02gr02/exp1_target_g_with_gci/230728/*.cbin.not.mat"
syllable = [
    "(?<=[^gMNO])gg+[N]g*(?<![^gMNO])",
    "(?<=[^gMNO])gg+[O]g*(?<![^gMNO])",
    "(?<=[^gMNO])ggg+(?<![^gMNO])",
]
labels = ["PB i", "PB c", "no PB"]

file_list = glob.glob(file_path)
seqs = hf.get_labels(file_list, ["E", "S"])
bouts, noise = hf.get_bouts(seqs, "g")

playback_list = []
for playback_idx in range(len(syllable)):
    all_instances_of_playback = re.findall(syllable[playback_idx], bouts)
    lengths = [len(item) for item in all_instances_of_playback]
    playback_list.append(lengths)

max_steps = np.max(np.hstack(playback_list))
steps = np.linspace(0, max_steps + 1, max_steps + 2)
# ----------------------------------------------------------------------------------------------------------------------
fig1, ax1 = plt.subplots(1, 1, figsize=(21 / 2.54, 12 / 2.54), sharex=True)
fig1.subplots_adjust(
    left=0.13, bottom=0.17, top=0.99, right=0.99, wspace=0.5, hspace=0.7
)
for plot_idx in range(len(playback_list)):
    counts, bins = np.histogram(playback_list[plot_idx], steps)
    ax1.plot(
        bins[1:-1],
        counts[1:] / np.sum(counts[1:]),
        label=labels[plot_idx],
        linewidth=4,
        color=colors_rep[plot_idx + 2],
    )
ax1.set_xlabel("Repeat number")
ax1.set_ylabel("Rel. frequency")
plt.legend()
fig1.savefig(f"figures\\pk02gr02_dist_0to10_repeats.svg")
fig1.savefig(f"figures\\pk02gr02_dist_0to10_repeats.png")

# ----------------------------------------------------------------------------------------------------------------------
fig1, ax1 = plt.subplots(1, 1, figsize=(21 / 2.54, 12 / 2.54), sharex=True)
fig1.subplots_adjust(
    left=0.13, bottom=0.17, top=0.99, right=0.99, wspace=0.5, hspace=0.7
)
for plot_idx in range(len(playback_list)):
    counts, bins = np.histogram(playback_list[plot_idx], steps)
    ax1.plot(
        bins[5:10],
        counts[5:10] / np.sum(counts[1:]),
        label=labels[plot_idx],
        linewidth=4,
        color=colors_rep[plot_idx + 2],
    )
ax1.set_xlabel("Repeat number")
ax1.set_ylabel("Rel. frequency")
plt.legend()
fig1.savefig(f"figures\\pk02gr02_dist_5to9_repeats.svg")
fig1.savefig(f"figures\\pk02gr02_dist_5to9_repeats.png")
plt.show()

#######################################################################################################################
file_path = "D:/Birds/or05pk04/Gulnur/exp*/*/*.cbin.not.mat"
syllable = [
    "(?<=[^bRT])b{1}[RT]b*(?<![^bRT])",
    "(?<=[^bRT])b{2}b*(?<![^bRT])",
    "(?<=[^bRT])b{2}[RT]b*(?<![^bRT])",
    "(?<=[^bRT])b{3}b*(?<![^bRT])",
    "(?<=[^bRT])b{3}[RT]b*(?<![^bRT])",
    "(?<=[^bRT])b{4}b*(?<![^bRT])",
]
labels = ["bW", "bb", "bbW", "bbb", "bbbW", "bbbb"]

file_list = glob.glob(file_path)
seqs = hf.get_labels(file_list, ["E", "S", "k", "f"])
bouts, noise = hf.get_bouts(seqs, "b")

playback_list = []
for playback_idx in range(len(syllable)):
    all_instances_of_playback = re.findall(syllable[playback_idx], bouts)
    lengths = [len(item) for item in all_instances_of_playback]
    playback_list.append(lengths)

max_steps = np.max(np.hstack(playback_list))
steps = np.linspace(0, max_steps + 1, max_steps + 2)
# ----------------------------------------------------------------------------------------------------------------------
all_counts = {}
fig1 = plt.figure(figsize=(21 / 2.54, 12 / 2.54))
for plot_idx, data_idx in zip([0, 1, 2], [0, 2, 4]):
    ax = fig1.add_subplot(3, 1, 1 + plot_idx)
    counts1, bins1 = np.histogram(playback_list[data_idx], steps)
    ax.plot(
        bins1[1:-1],
        counts1[1:] / np.sum(counts1[1:]),
        label=labels[data_idx],
        linewidth=4,
        color=colors_2[0],
    )
    counts2, bins2 = np.histogram(playback_list[data_idx + 1], steps)
    ax.plot(
        bins2[1:-1],
        counts2[1:] / np.sum(counts2[1:]),
        label=labels[data_idx + 1],
        linewidth=4,
        color=colors_2[1],
    )
    plt.legend()
    all_counts[labels[data_idx]] = counts1
    all_counts[labels[data_idx + 1]] = counts2
    if plot_idx < 2:
        plt.tick_params(axis="x", which="both", bottom=False, labelbottom=False)
ax.spines["bottom"].set_visible(True)
fig1.savefig(f"figures\\or05pk04_repeat_dist.svg")
fig1.savefig(f"figures\\or05pk04_repeat_dist.png")
plt.show()

np.save("figures/or05pk04_repeat_dist.npy", all_counts, allow_pickle=True)
embed()
quit()
