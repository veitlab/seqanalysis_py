import sys
import glob
import numpy as np
import vocalpy as voc
import matplotlib.pyplot as plt

from IPython import embed
from pathlib import Path
from seqanalysis.util.filter_functions import envelop, bandpass_filter


def main(root_path):
    # read wav/notmat data and import segments
    audio_paths = glob.glob(root_path + '\*.wav')
    notmat_paths = glob.glob(root_path + '\*.wav.not.mat')

    sounds = [voc.Sound.read(audio_path) for audio_path in audio_paths]
    annots = [voc.Annotation.read(notmat_path, format='notmat') for notmat_path in notmat_paths]

    # get onsts of the klicks
    for file_idx in range(len(sounds)):
        labels = annots[file_idx].data.seq.labels
        offset = annots[file_idx].data.seq.offsets_s
        onset = annots[file_idx].data.seq.onsets_s

        sr = sounds[file_idx].samplerate
        y_data = sounds[file_idx].data[0]
        env, envrate = envelop(bandpass_filter(y_data, sr, 500.0, 10000.0), sr, 300.0)
        env_time = np.arange(len(env)) / envrate

        plots_per_figure = 5
        if file_idx % plots_per_figure == 0:
            fig, axes = plt.subplots(plots_per_figure, 1, figsize=(10, 20), sharex='all', sharey='all', frameon=False)

        axes[file_idx % plots_per_figure].plot(env_time, env, 'k')
        axes[file_idx % plots_per_figure].plot(onset, [0.01]*len(onset), '.', color='red')
        axes[file_idx % plots_per_figure].plot(offset, [0.01]*len(offset),'.', color='blue')

        axes[file_idx % plots_per_figure].tick_params(tick1On=False, tick2On=False, label1On=False, label2On=False)
        axes[file_idx % plots_per_figure].axis('off')
        axes[file_idx % plots_per_figure].text(-0.1, 0.5, Path(annots[0].path).name, fontsize='small')

    plt.show()


if __name__ == "__main__":
    main(sys.argv[1])
