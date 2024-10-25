import re
import sys
import glob
import numpy as np
import pandas as pd
import vocalpy as voc
import matplotlib.pyplot as plt

from IPython import embed
from seqanalysis.util.get_data_transition_diagram import get_labels, get_bouts
from seqanalysis.util.filter_functions import highpass_filter, envelop, bandpass_filter, threshold_crossings, threshold_crossing_times


def get_catch_data(path, intro_notes, bout_chunk):
    list = glob.glob(path)

    file_list = []
    for i in range(len(list)):
        with open(list[i], 'r') as file:
            line_list = file.readlines()
            file_list.extend([list[i].rstrip('batch.notcatch') + item.rstrip() + '.not.mat' for item in line_list])

    seqs = get_labels(file_list, intro_notes)
    bouts, _ = get_bouts(seqs, bout_chunk)

    return bouts

def main(root_path):

    # read wav/notmat data and import segments
    audio_paths = glob.glob(root_path + '*.wav')
    notmat_paths = glob.glob(root_path + '*.wav.not.mat')

    sounds = [voc.Sound.read(audio_path) for audio_path in audio_paths]
    annots = [voc.Annotation.read(notmat_path, format='notmat') for notmat_path in notmat_paths]

    for file_idx, file in enumerate(sounds):
        labels = annots[0].data.seq.labels
        offset = annots[file_idx].data.seq.offsets_s[labels=='2']
        onset = annots[file_idx].data.seq.onsets_s[labels=='2']

        sr = sounds[file_idx].samplerate
        y_data = sounds[file_idx].data
        env, envrate = envelop(bandpass_filter(y_data, sr, 500.0, 10000.0), sr, 300.0)
        env_time = np.arange(len(env)) / envrate

        for on, off in zip(onset, offset):
            plt.plot(env[(env_time >= on) & (env_time <= off)])

    plt.show()

    embed()
    quit()
    
    up_cross, down_cross = threshold_crossings(env, threshold)
    up_t, down_t = threshold_crossing_times(env_time, env, threshold, up_cross, down_cross)


if __name__ == "__main__":

    main(sys.argv[1])