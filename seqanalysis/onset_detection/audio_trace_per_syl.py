import sys
import glob
import numpy as np
import vocalpy as voc
import matplotlib.pyplot as plt

from IPython import embed
from seqanalysis.util.filter_functions import envelop, bandpass_filter, threshold_crossings, threshold_crossing_times


def main(root_path, target_syl, threshold):
    # read wav/notmat data and import segments
    audio_paths = glob.glob(root_path + '\*.wav')
    notmat_paths = glob.glob(root_path + '\*.wav.not.mat')

    sounds = [voc.Sound.read(audio_path) for audio_path in audio_paths]
    annots = [voc.Annotation.read(notmat_path, format='notmat') for notmat_path in notmat_paths]

    # get onsts of the klicks
    onset_klick = []
    for file_idx in range(len(sounds)):
        labels = annots[file_idx].data.seq.labels
        offset = annots[file_idx].data.seq.offsets_s[labels == target_syl]
        onset = annots[file_idx].data.seq.onsets_s[labels == target_syl]

        sr = sounds[file_idx].samplerate
        y_data = sounds[file_idx].data[0]
        env, envrate = envelop(bandpass_filter(y_data, sr, 500.0, 10000.0), sr, 300.0)
        env_time = np.arange(len(env)) / envrate
        # plt.plot(env_time, env)
        #
        # for on, off in zip(onset, offset):
        #     plt.plot(env_time[(env_time >= on) & (env_time <= off)], env[(env_time >= on) & (env_time <= off)])
        # plt.show()

        for on, off in zip(onset, offset):
            env_syl = env[(env_time >= on) & (env_time <= off)]
            env_time_syl = env_time[(env_time >= on) & (env_time <= off)]
            plt.plot(env_time_syl - env_time_syl[0], env_syl)
            up_cross, down_cross = threshold_crossings(env_syl, threshold)
            up_t, down_t = threshold_crossing_times(env_time_syl - env_time_syl[0], env_syl, threshold, up_cross, down_cross)
            onset_klick.extend(up_t)
            print(file_idx)
            print(up_t)
            plt.plot(up_t, [threshold] * len(up_t), 'o', 'k')

    plt.show()
    embed()
    quit()
    fig1, ax1 = plt.subplots(1, 1, figsize=(21 / 2.54, 12 / 2.54), sharex=True)
    fig1.subplots_adjust(left=0.13, bottom=0.17, top=0.99, right=0.99, wspace=0.5, hspace=0.7)

    counts, bins = np.histogram(onset_klick)
    ax1.stairs(counts, bins, linewidth=4, color='pink')
    # ax1.hist(bins[:-1], bins, weights=counts)

    ax1.set_xlabel("clic onset time", fontsize=25)
    ax1.set_ylabel("Rel. frequency", fontsize=25)
    plt.rcParams["font.size"] = 25
    plt.rcParams["axes.linewidth"] = 2
    ax1.spines["top"].set_visible(False)
    ax1.spines["right"].set_visible(False)
    ax1.tick_params(width=2)
    ax1.tick_params(axis="both", which="major", labelsize=20)

    # fig1.savefig(f"{syl}_repeats.pdf")

    plt.show()



if __name__ == "__main__":
    main(sys.argv[1], sys.argv[2], float(sys.argv[3]))
