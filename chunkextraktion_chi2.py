import re
import sys
import glob
import yaml
import numpy as np
import scipy.io as sio
import scipy.stats as stats
import plot_functions as pf
import helper_functions as hf
import matplotlib.pyplot as plt

from IPython import embed


def get_labels(mat_list, notes):
    seqs = []
    for matidx in mat_list:
        mat = sio.loadmat(matidx)

        labels = mat['labels'][0]
        labels = 'S' + labels + 'E'

        if len(notes) > 0:
            try:
                labels = replace_intro_notes(labels, notes)
            except:
                print(matidx)

        seqs.append(labels)
    seqs = np.array(seqs)

    return seqs


def replace_intro_notes(s, intro_notes):
    unique_labels = sorted(list(set(s)))
    for i in range(len(intro_notes)):
        if intro_notes[i] in s:
            unique_labels.remove(intro_notes[i])

    motiv_start = []
    for i in range(len(unique_labels)):
        motiv_start.append(s.find(unique_labels[i]))

    temp = list(s)
    for i in range(1, np.min(motiv_start)):
        temp[i] = 'i'

    s = ''.join(temp)

    return s


def chi2_dist(dist1, dist2, alpha):
    '''

    :param dist1: observed distribution
    :param dist2: expected distribution
    :param alpha: without it alpha at 0.01

    :return:

    '''

    df = len(dist1) - 1  # freiheitsgrade
    obs_dist = np.array([dist1, dist2])

    sum_col = np.sum(obs_dist, axis=0)
    sum_row = np.sum(obs_dist, axis=1)
    sum_total = np.sum(obs_dist)
    exp_dist = np.zeros(np.shape(obs_dist))

    for idx_col in range(np.shape(obs_dist)[1]):
        for idx_row in range(np.shape(obs_dist)[0]):
            exp_dist[idx_row, idx_col] = sum_row[idx_row]*sum_col[idx_col]/sum_total

    diff = obs_dist - exp_dist
    norm_diff = (diff**2)/exp_dist
    chi2 = np.sum(norm_diff)
    p = 1- stats.chi2.cdf(chi2, df)

    if p < alpha:
        h = 1  # null hypothesis
    else:
        h = 0

    return h, np.round(chi2, 2), np.round(p, 4), df, exp_dist


def get_data(cfg):
    # file_list = glob.glob(cfg['paths']['folder_path'])

    # seqs = get_labels(file_list, ['S', 'E', 'm', 'n', 'b'])
    # seqs = ''.join(set(seqs))
    seqs = 'YiiiiiiiiiiiidabjjgekhcmmndabjjgekhcnkhcnkhcmmndabjjgekhcndabjjgekhcndabjjgekhcnkhcmmndabjjgekhcmmnndabjYiiiiiiiiiiiiiidabjjgekhcnkhcmmndabjjgekhcnkhcmmndabjjgekmmmmnndabjjgekhcmmnmmndabjjgekhcmmndabmmmYiiiiiiiiiiiiiiiidabjjgekhcnkhcmmndabjjgekhcnkhcmmndabjjgekhcmmmmnndabjjggkhcmmYiiiiiiiiiiiiiiidabjjgekhcmmnndabjjgekhcmmnndabjjgekhcmmndabjjgekhcnkhcmmndabjjgekhcnkhcmmnndabjjgekhcmmmmnndabjjgekhcmYiiiiiiiiiiiiiidabjjgekhcmmnndabjjgekhcmmndabjjgekhcnkhcmmndabjjgekhcnkhcndabjjgekhcmmmnndabjjjekhcmmYiiiiiiiiiiiiiidabjjgekhcmmnnndabjjgekhcnkhcmmnndabjjgekhcnkhcnnhbmmnndabjjgekhcnkhcmmnndabjjgekhcmmnndabjjgekhcmmmnndabjjgekhcmYiiiiiiiiiiiiiiiiidabjjgekhcnkhcmmmmndabjjgekhcmmmndabjjgekhcnkhcmmndabjjgeYiiiiiiiiiiidabjjgekhcnkhcmmndabjjgemYiiiiiiiiiiiiidabjjgekhcnkhcmmndabjjgeYiiiiiiiiiiiiidabjjgekhcnkhcmmnndabjjgekhcnkhcnkhcmmndabjjgekhcmmmnnndabjjgekhcmmmmnndabjjgekhcmnkhcnkhcmmndabjjgekhcmmYiiiiiiiiiiiiidabjjgekhcnkhcmmmmmndabjjgekhcmmmndabjjgekhcndabjjgekhcmmndabjjgekhcmmndabjjgekhcnkhcmmndabjjgekhcndabjjgemYiiiiiiiiiiiiiiiidabjjgekhcmmndabjjgekhcnkhcmmndabjjgekhcmmmmmndabjjgekhcmmmnndabjjgekhcmmmmnndabjjgekhcmYiiiiiiiiiiiiiiiiiiidabjjgekhcnkhcmmndabjjgekhcnkhcmmnndabjjgekhcnnhbndabjjgekhcnkhcmmnndabjjgekhcnkhcmmndabjjgekhcmmmmnndabjjggnhbmmmmnndabjjgeYiiiiiiiiiiiiiidabjjgekhcnkhcmmnndabjjgekhcnkhcmmndabjjgekhcnkhcmmndabjjgekhcnkhcmmmnndabjjgekhcmmmmndabjjgekhcmmYiiiiiiiiiiiiiiiidabjjgekhcndabjjgekhcnkhcmmndabjjgekhcnkhcmmndabjjgekhcmmYiiiiiiiiiiiidabjjgekhcnkhcmmndabjjgekhcnkhcmmndabjjgekhcnkhcmmndabjjgekhcnkhcmmndabjjgekhcmmnndabjjgekhcmmmmmnndabjjgekhcmmmmnndabjjgeYiiiiiiiiiiidabjjgekhcnkhcmmndabjjgekhcnkhcmmnndabjjgekhcnkhcmmndabjjgekhcnkhcmmnndabjjgekhcmmmndabjjgekhcmmYiiiiiiiiiiiiiidabjjgekhcnkhcmmndabjjgekhcnkhcmmndabjjgekhcnkhcmmndabjjgekhcnkhcmmnndabjjgekhcmmmndabjjgekhcmmndabmYiiiiiiiiiiiiiiiidabjjgekhcndabjjgekhcnkhcndabjjgekhcnkhcmmndabjjgekhcndabjjgeYiiiiiiiiiiiiiiiidabjjgekhcmmnndabjjgekhcnkhcmmndabjjggkhcnkhcmmnndabjjgekhcnkhcmmmnndabjjgekhcnkhcmmmmnndabjjgekhcmYiiiiiiiiiiiiidabjjgekhcnkhcmmndabjjgekhcnkhcmmmnndabjjgeYiiiiiiiiiiiiiiiidabjjgekhcmmndabjjgekhcmmndabjjgekhcmmnndabjjgekhcnkhcmmndabjjgekhcmYiiiiiiiiiiiiidabjjgekhcnkhcndabjjgekhcnkhcmmndabjjgekhcmmnndabjjgekhcmYiiiiiiiiiiiiiiiiiidabjjgekhcndabjjgekhcmmndabjjgekhcmmmYiiiiiiiiiiiiiiiidabjjgekhcmmmmmndabjjgemmmmmndabjjgekhcmmmmmnndabjjgekhcmYiiiiiiiiiiiiiiiidabjjgekhcmmndabjjgekhcmmnndabjjgekhcmmmmmndabjjgekhcnkhcmmmmndabjjgeYiiiiiiiiiiiiiidabjjgekhcnkhcmmndabjjgekhcmmmmmnndabjjgekhcnkhcmmndabmYiiiiiiiiiiiiidabjjgekhcnkhcmmndabjjgekhcmmmmYiiiiiiiiiiiidabjjgekhcnkhcmmndabjjgekhcnkhcmmmnndabjjgekhcmYiiiiiiiiiiiiiiiiiiidabjjgekhcmmndabjjgekhcnkhcmmndabjjgekhcnkhcmmnnndabjjgekhcmmmnndabjjgekhcnkhcmmYiiiiiiiiiiiiiiidabjjgemhbmmYiiiiiiiiiiidabjjgekhcnkhcmmndabjjgekhcmmmnndabjjgekhcnkhcmmndabjjgekhcm'
    # unique syl
    seqs = hf.replace_chunks(seqs, ['i+'])
    unq = sorted(list(set(seqs)))

    # prob of all syls to all possible following syl
    tm_sf, tmp_sf = hf.get_transition_matrix(seqs, unq)
    sylM_sf = np.array([['11'] * len(unq)] * len(unq))
    for idx1 in range(len(unq)):
        for idx2 in range(len(unq)):
            sylM_sf[idx1, idx2] = unq[idx1] + unq[idx2]

    # prob of all syls to all possible following syl and depending on before syl
    tm_bsf, tmp_bsf = hf.get_transition_matrix_befor_following_syl(seqs, unq)
    sylM_bsf = np.array([[['111'] * len(unq)] * len(unq)] * len(unq))
    for idx1 in range(len(unq)):
        for idx2 in range(len(unq)):
            for idx3 in range(len(unq)):
                sylM_bsf[idx2, idx1, idx3] = unq[idx2] + unq[idx1] + unq[idx3]

    # Chi2 Test -------------------------------------------------------------------------------------------------------

    # remove incredibly small branches which are likely mislabels
    tm_sf2 = tm_sf
    tm_sf2[tmp_sf[:] < 0.01] = 0

    # remove all lines with only one entry point
    nzeros = [len(np.nonzero(x)[0]) > 1 for x in tm_sf2]
    tm_sf2[~np.array(nzeros)] = 0

    # save tm_bsf
    test_tm_bsf = tm_bsf
    sumcol = np.sum(np.sum(tm_bsf, axis=2), axis=0)

    # chi2 test
    h = np.full((len(unq), len(unq)), np.nan)
    chi = np.full((len(unq), len(unq)), np.nan)
    p = np.full((len(unq), len(unq)), np.nan)

    # go through all cells and delete whole if sum within the list is less
    # than 1% of the total times this syllable is observed
    for idx1 in range(len(test_tm_bsf)):
        for idx2 in range(len(test_tm_bsf)):
            if np.sum(test_tm_bsf[idx1, idx2]) / sumcol[idx2] <= 0.01:
                test_tm_bsf[idx1, idx2] = 0

    #
    test_sylM_bsf = np.array([[['   '] * len(unq)] * len(unq)] * len(unq))
    for idx1 in range(len(unq)):
        existingpostsyls = np.argwhere(tm_sf2[idx1]/np.sum(tm_sf2[idx1]) > 0.01)
        if len(np.argwhere(tm_sf2[idx1] > 0.01)) > 0:
            for idx2 in range(len(unq)):
                numstates = len(np.unique(np.nonzero(test_tm_bsf[:, idx1])[0]))
                if np.sum(test_tm_bsf[idx2, idx1]) > 0:
                    if np.count_nonzero(test_tm_bsf[idx2, idx1][existingpostsyls]) > 1:
                        test_sylM_bsf[idx2, idx1][existingpostsyls] = sylM_bsf[idx2, idx1][existingpostsyls]
                        h_chi2, chi_chi2, p_chi2, _, _ = chi2_dist(np.hstack(tm_sf2[idx1][existingpostsyls]),
                                                                   np.hstack(test_tm_bsf[idx2, idx1][existingpostsyls]),
                                                                   0.01 / numstates)
                        h[idx2, idx1] = h_chi2
                        chi[idx2, idx1] = chi_chi2
                        p[idx2, idx1] = p_chi2

                    elif np.count_nonzero(test_tm_bsf[idx2, idx1][existingpostsyls]) == 1:
                        h[idx2, idx1] = 1


    embed()
    quit()
    return cfg


# def make_first_plots(cfg):
#     bouts = cfg['data']['bouts_rep']
#     tm, _ = hf.get_transition_matrix(bouts, cfg['labels']['unique_labels'])
#
#     # ---------------------------------------------------------------------------------------------------------------
#     k = np.where(sum(tm) / sum(sum(tm)) * 100 <= 0)
#     tmd = np.delete(tm, k, axis=1)
#     tmd = np.delete(tmd, k, axis=0)
#
#     tmpd = (tmd.T / np.sum(tmd, axis=1)).T
#     tmpd = hf.get_node_matrix(tmpd, cfg['constants']['edge_threshold'])
#     'Plot Transition Matrix and Transition Diagram'
#     node_size = np.round(np.sum(tmd, axis=1) / np.min(np.sum(tmd, axis=1)), 2) * 100
#     pf.plot_transition_diagram(tmpd, np.delete(cfg['labels']['unique_labels'], k),
#                                node_size,
#                                cfg['paths']['save_path']+cfg['title_figures']+'_graph_simple.jpg',
#                                cfg['title_figures'])
#     plt.show()


def main(yaml_file):
    with open(yaml_file) as f:
        cfg = yaml.load(f, Loader=yaml.FullLoader)
        f.close()

    cfg = get_data(cfg)

    # make_first_plots(cfg)
    print('\n Unique labels of Chunks: {}\n'.format(sorted(list(set(cfg['data']['chunk_bouts'])))))

    # with open(yaml_file, 'w') as f:
    #     yaml.dump(cfg, f)
    #     # print(yaml.dump(cfg))
    #     f.close()
    # embed()
    # quit()


if __name__ == '__main__':
    #     This function tried to automate extraction of chunks based on chi sq
    # analysis. First, we try to determine if a syllable's next transition
    # depends on the syllable that comes before it using chi sq analysis. If it
    # does, we relable the syllable and use it as a different 'state' of the
    # same syllable. Then using my chunk exrtraction function, we create chunks
    # with one-in-one-out branches, >80% transition prob as middle nodes.
    # Detailed notes are in the description of the function. This function
    # gives you chunks and plots old transition diagram and new transition
    # diagram automatically. Input= seqnew, a clean sequence string with intro
    # notes replaced and start symbol 'Y' already present if you dont want a
    # plot write 0 for plotting you need seqforchunks and chunks2replace for
    # chunk consistency analysis and labelidx also newseq gives you the seq w/
    # chunks replaced divprobtoplot2 = final transitionprob matrix and
    # patterncell2 gives you the patterncell for the diagraph labels2 = what
    # are the replaced labels numnewsyls = number of additional states
    #
    # INPUT:
    # sys.argv[1] = yaml file of the bird, example: example_yaml.yaml
    # sys.argv[2] = analysis catch or all files: input: catch, all
    #
    # OUTPUT:
    # figures

    main(sys.argv[1])
