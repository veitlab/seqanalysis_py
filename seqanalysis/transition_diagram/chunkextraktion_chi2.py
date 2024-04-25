import re
import sys
import glob
import yaml
import string
import numpy as np
import scipy.io as sio
import scipy.stats as stats
import seqanalysis.util.plot_functions as pf
import seqanalysis.util.helper_functions as hf
from seqanalysis.util.logging import config_logging
import matplotlib.pyplot as plt

log = config_logging()

from IPython import embed


def chi2_dist(dist1, dist2, alpha):
    """

    :param dist1: observed distribution
    :param dist2: expected distribution
    :param alpha: without it alpha at 0.01

    :return:

    """

    # NOTE: the degrees of freedoms should be (number of rows - 1)(number cols -1)
    # but the cols are always 2 so its is (len(dist1)-1)*1

    df = len(dist1) - 1  # freiheitsgrade
    obs_dist = np.array([dist1, dist2])

    sum_col = np.sum(obs_dist, axis=0)
    sum_row = np.sum(obs_dist, axis=1)
    sum_total = np.sum(obs_dist)
    exp_dist = np.zeros(np.shape(obs_dist))

    for idx_col in range(np.shape(obs_dist)[1]):
        for idx_row in range(np.shape(obs_dist)[0]):
            exp_dist[idx_row, idx_col] = sum_row[idx_row] * sum_col[idx_col] / sum_total

    diff = obs_dist - exp_dist
    norm_diff = (diff**2) / exp_dist
    chi2 = np.sum(norm_diff)
    p = 1 - stats.chi2.cdf(chi2, df)

    if p < alpha:
        h = 1  # reject the null hypthesis
    else:
        h = 0

    return h, np.round(chi2, 2), np.round(p, 4), df, exp_dist


def make_states_of_syl_chi2(seqs):
    # unique syl
    unq = sorted(list(set(seqs)))

    embed()
    exit()
    # prob of all syls to all possible following syl
    tm_sf, tmp_sf = hf.get_transition_matrix(seqs, unq)
    sylM_sf = np.zeros((len(unq), len(unq)), "U2")

    for idx1, syly in enumerate(unq):
        for idx2, sylx in enumerate(unq):
            sylM_sf[idx1, idx2] = syly + sylx

    # prob of all syls to all possible following syl and depending on before syl
    tm_bsf, _ = hf.get_transition_matrix_befor_following_syl(seqs, unq)
    sylM_bsf = np.zeros((len(unq), len(unq), len(unq)), "U3")
    for idx1, syly in enumerate(unq):
        for idx2, sylx in enumerate(unq):
            for idx3, sylz in enumerate(unq):
                # BUG: this was idx2, idx1, idx3
                sylM_bsf[idx1, idx2, idx3] = syly + sylx + sylz

    # clean up matrix and remove small branches and single lines ------------------------------------------------------

    # remove incredibly small branches which are likely mislabels
    tm_sf2 = tm_sf
    tm_sf2[tmp_sf < 0.01] = 0

    # remove all lines with only one entry point
    nonzeros_log = [len(np.nonzero(x)[0]) > 1 for x in tm_sf2]
    tm_sf2[~np.array(nonzeros_log)] = 0

    # remove Y transitions because Y doesn't depend on anything
    coly = unq.index("Y")
    tm_sf2[:, coly] = 0

    # save tm_bsf
    test_tm_bsf = tm_bsf
    # WARNING: is this the correct way to sum the columns? shoulnd it be axis=1 istead of axis=0?
    sumcol = np.sum(np.sum(tm_bsf, axis=2), axis=0)

    # NOTE: chi2 Test
    h = np.full((len(unq), len(unq)), np.nan)
    chi = np.full((len(unq), len(unq)), np.nan)
    p = np.full((len(unq), len(unq)), np.nan)

    # go through all cells and delete whole if sum within the list is less
    # than 1% of the total times this syllable is observed
    for idx1 in range(len(test_tm_bsf)):
        for idx2 in range(len(test_tm_bsf)):
            # BUG:this should not be the sum of but every element in the row
            if np.sum(test_tm_bsf[idx1, idx2]) / sumcol[idx2] <= 0.01:
                test_tm_bsf[idx1, idx2] = 0

    # BUG: the Code is not deleting transitions that come from the intro syll in the tm_bsf

    # NOTE: suggested code
    # for idx1 in range(len(unq)):
    #     for idx2 in range(len(unq)):
    #         for idx3 in range(len(unq)):
    #             if test_tm_bsf[idx2, idx2, idx3] / sumcol[idx2] <= 0.01:
    #                 test_tm_bsf[idx1, idx2, idx3] = 0
    # and prettier but dont know if this is better or even works
    # for idx1 in range(len(unq)):
    #   test_tm_bsf[:, idx1][(test_tm_bsf[:, idx1]/sumcol[idx1])<= 0.01] = 0

    test_sylM_bsf = np.zeros((len(unq), len(unq), len(unq)), "U3")
    for idx1 in range(len(unq)):
        # NOTE: Adding check if the sum of the row is 0 to avoid division by zero
        if np.sum(tm_sf2[idx1]) == 0:
            continue
        existingpostsyls = np.argwhere(tm_sf2[idx1] / np.sum(tm_sf2[idx1]) > 0.01)
        numstates = len(np.unique(np.nonzero(test_tm_bsf[:, idx1])[0]))
        # NOTE: this is should produce the same result for entering the if statement?
        # if len(existingpostsyls.shape[0]) > 0:
        if len(np.argwhere(tm_sf2[idx1] > 0.01)) > 0:
            for idx2 in range(len(unq)):
                if np.sum(test_tm_bsf[idx2, idx1]) > 0:
                    if np.count_nonzero(test_tm_bsf[idx2, idx1][existingpostsyls]) > 1:
                        h_chi2, chi_chi2, p_chi2, _, _ = chi2_dist(
                            np.hstack(tm_sf2[idx1][existingpostsyls]),
                            np.hstack(test_tm_bsf[idx2, idx1][existingpostsyls]),
                            0.01 / numstates,
                        )
                        # res = stats.chisquare(
                        #     np.array(
                        #         [
                        #             np.hstack(tm_sf2[idx1][existingpostsyls]),
                        #             np.hstack(
                        #                 test_tm_bsf[idx2, idx1][existingpostsyls]
                        #             ),
                        #         ]
                        #     ),
                        #     axis=None,
                        # )
                        # log.info(f"{res}\n")
                        # log.info(h_chi2)
                        if h_chi2 == 1:
                            test_sylM_bsf[idx2, idx1][existingpostsyls] = sylM_bsf[
                                idx2, idx1
                            ][existingpostsyls]

                        h[idx2, idx1] = h_chi2
                        chi[idx2, idx1] = chi_chi2
                        p[idx2, idx1] = p_chi2
                        log.info(f"{h_chi2}\n{p_chi2}")

                    elif (
                        np.count_nonzero(test_tm_bsf[idx2, idx1][existingpostsyls]) == 1
                    ):
                        h[idx2, idx1] = 1

    # keep only the transitions that are h = 1
    test_sylM_sf = sylM_sf
    test_sylM_sf[h != 1] = ""

    # now relabeling of original sequence -----------------------------------------------------------------------------
    # find letters not present in seqs

    embed()
    exit()
    all_letters = string.ascii_letters + string.digits
    available_letters = all_letters[::-1]
    for unq_idx in unq:
        available_letters = re.sub(unq_idx, "", available_letters)

    # collecting indices for repleacing and replacing together
    states_log = h == 1
    states_sum = sum(states_log)
    to_relabel = np.delete(test_sylM_sf, np.where(states_sum == 0), axis=1)

    relabel_table = np.full((int(sum(states_sum)), 3), "  ")
    relabel_pos = []

    av_letter_counter = 0
    for states_sum_idx in range(sum(states_sum != 0)):
        to_relabel_row = np.where(to_relabel[:, states_sum_idx] != "  ")[0]
        row_counter = 0
        for to_relabel_row_idx in to_relabel_row:
            state_seq = to_relabel[to_relabel_row_idx, states_sum_idx]
            relabel_pos.append(
                [match.end(0) - 1 for match in re.finditer(state_seq, seqs)]
            )

            relabel_table[av_letter_counter, 0] = np.array(
                available_letters[av_letter_counter], dtype="S1"
            )
            relabel_table[av_letter_counter, 1] = state_seq[1] + str(row_counter)
            relabel_table[av_letter_counter, 2] = state_seq

            row_counter += 1
            av_letter_counter += 1

    # collecting indices for replacing and replacing together
    relabel_tm_bsf = np.delete(test_tm_bsf, np.where(states_sum == 0), axis=1)

    which_to_merge = np.full(relabel_tm_bsf[:, :, 0].shape, np.nan).tolist()
    names_to_merge = np.full(relabel_tm_bsf[:, :, 0].shape, np.nan).tolist()

    for ip in range(np.shape(which_to_merge)[1]):
        # this is to keep track of the actual position in the matrix
        f_row, f_column = np.where(
            relabel_tm_bsf[:, ip, :]
            / np.sum(relabel_tm_bsf[:, ip, :], axis=1, keepdims=True)
            > 0.05
        )

        # to get which branches have the same transitions and compare them, if the state of the syllable is the same
        # or different depending on the previous syllable
        # value_idx = np.unique(np.nonzero(relabel_tm_bsf[:,ip,:])[0])
        # percent_of_iter = relabel_tm_bsf[value_idx, ip, :]/np.sum(relabel_tm_bsf[value_idx, ip, :], axis=1, keepdims=True)

        mapping = {}
        for key, value in zip(f_row, f_column):
            if key in mapping:
                mapping[key].append(value)
            else:
                mapping[key] = [value]

        out = []
        for key1, value1 in mapping.items():
            for key2, value2 in mapping.items():
                if key1 < key2 and value1 == value2:
                    out.append(value1)

        embed()
        quit()
        for ir in range(len(out)):
            index = []
            for k in np.unique(f_row):
                if mapping[k] == out[ir]:
                    ind = mapping[k] == out[ir]
                    index.append(k)

            # if len(out[ir]) > 1: #then do all the testing

        # ToDo: this for loop will not work because i'm working from the inside out of the loops in the matlab
        #  function (line 178-216)

    return seqs

    # for ir in range(len(out)):
    #     index = []
    #     for k in range(len(f)):
    #         if out[ir] in f[k]:
    #             index.append(k)
    #     if len(out[ir]) > 1:
    #         pairwise_combinations = np.array(list(combinations(index, 2)))
    #         numstates2 = len(pairwise_combinations) // 2
    #         for id in range(pairwise_combinations.shape[0]):
    #             data1 = remcountstg[pairwise_combinations[id,0],ip]
    #             data2 = remcountstg[pairwise_combinations[id,1],ip]
    #             h_testing_in_pairs = 0 if chi2.pvalue(chi2.stat(data1, data2), 1) > 0.01/numstates2 else 1
    #             if h_testing_in_pairs == 0:
    #                 if 'which_to_merge' not in locals():
    #                     which_to_merge = []
    #                 which_to_merge.append([pairwise_combinations[id,0],ip,pairwise_combinations[id,1],ip])
    #     elif len(out[ir]) == 1:
    #         if 'which_to_merge' not in locals():
    #             which_to_merge = []
    #         which_to_merge.extend([[i,ip] for i in index])
    # if 'which_to_merge' in locals():
    #     which_to_merge = np.unique(which_to_merge, axis=0)
    #     if 'names_to_merge' not in locals():
    #         names_to_merge = []
    #     names_to_merge.append([torelabel[row,col] for row, col in which_to_merge])


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

    sequence = cfg["data"]["bouts"]
    chunk_seqs = make_states_of_syl_chi2(sequence)


if __name__ == "__main__":
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
    #
    # OUTPUT:

    # ToDo: this function is not finished yet

    seqs = "YimcBljljmcBgmjkcljljmcjmcBtjljljmcBgmjkcYimcBjjljljmcBgmjkcljljljmcBgmjkcljljljmcBgmjkcljljljmcBgmjkcljljmcBgmjkcljljljmcBgmjkcljljljmcBgmjkcljljjmcBgYimcBjljmcBgmjkcljljljmcBgmjkcljmcBgmjkcljljljmcBgmjkcljljljjcbbbbbbbgmjkcjljljmcjljljmcBgYimcBfljljmcBgmjkcljljljmfjjjmcBgmjkcljmcBgmjkcljljljmjmcBgmjkcljljljmcBgmjkcljljljmcBgmjkcljljljmcBgYimcBjljljmcBgmjkcljljYimcBjljljmcBgmjkcljljljmcBgmjkcljljljmcBgmjkcljljljmcBgmjkcljljmcBgmjkcljljljmcBYimcBjljmcBgmjkcljljljmcBgmjkcljljljmcBgmjkcljljljmcBtjljmjmcBgmjkcljljljljmcBgmjkcljljljmcBgYicbbbbbbjljmcBgmjkcljmcBgmjkcljljmjmcBgmjkcljljmcBgmjkcljljljjcbbbbgljkcljljljjcbbbbbbgmjkcYlflflalflfljmcBjljljljmjmcBgmjkcljljljmcBgmjkcljljljmcBgmjkcljljljmcBgmjkcljmcBgmjkcljljljmcBgjljljmcBgmjkcljljljmcBgYiabbbbbbbbjljljmcBgkjkcljkjljljljmcBfjljljmjmcBgmjkcljljljmcBgmjkcljljljljjmcBgmjkcljljljljmcBgmjkcljljljmcBgjljljjcbbbbbgYimjmcBjljljljmcBgmjkcljljljmcBgmjkcljljkmcBjljmjmcBgmjkcljljljmcBgmjkcljljljmcBgmjkcljljljmcBgmjkclYimcBjljljmcBgmjkcljljljljmcBgmjkcljljljmcBgmjkcljljmcBgmjkcljljmcBtjljmjmcBgmjkcljljljjjcbbbgjljljmcBgjljjmcBgYimcBjljljmcBgmjkcljljljmcBgmjkcljljljmcBgmjkcljljmcBgmjkcljlfljmcBgmjkcljljmcBgljkcljljmcBYimcBjljmcBgmjkcljmcBgmjkcljljfjmcBgmjkcljmcBgmjkcljljljljmcBgYimcBflfmcBgmjkcjljljmcBgmjkcljmcBgkjkcljljljljmcBgmjkcljljljljmcBgmjkcljljljmcBgmjkclYimcBjljljmcBgmjkcljljlflfmcBgljkcljlfljljmcBgmjkcljljmcBgmjkcljlfljmcBjljmccgmjkcljlfljljljmcBgmjkcljljlfjmcBYimcBjljljmcBgmjkcljljlfljmcBgmjkcljmcBgmjkcljljljljmcBYimcBjljljjcbbbgmjkcljljljmcBgmjkcljljljmcBgmjkcljmcBgmjkcljljljmcBgmjkcljljmcBgjljmcBgjljljmcBgmjkcljljljmcBYimcBjljmcBgmjkcljljljmjmcBgmYimcBjljljmcBgmjkcljlfljmcBgmjkcljljljmcBjljmjmccgmjkcljmcBfjljljmcBgmjkcljljljmcBjjljljmjmcBgmjkcljljljmcBgjljmcBgYimjmcBjlljmcBgmjkcljljljmcBgmjkcljljmcBgmjkcljljmcBgYimcBjljljacbbbbgmjkcljljljmcBgmjkcljljmYimcBjjljljmcBgmjkcljmcBgmjkcljmcBgmjkcljmcBgmjkcljljmcBgmjkcljljljmjmcBgmjkcljljmcBYimcBjljljmcBgmjkcljljljmcBjjmcBgmjkcljljljjcbbbbjljljmcBgmjkcljljljmcBgmjkcljljljmcBgmjkcljmcBgYimcBjljljmcBgmjkcljljljmjmcBgmjkcljmcBgmjkcljmcBgmjkcljlfjljmcBgmjkcljmcBgmjkcljljljmcBgYimcBjljljmcBgmjkcljmcBgmjkcljljljmcBgmjkcljljljmcjmcBgmjkcljljljmcBgjljljmcBgYimfmjmcBjljljmcBgmjkcljljmcBjljmjmcBgmjkcljljljmclgmjkcljmcBgmjkcljljmcBgmjkcljljmcBgmjkcljlfljmcBgmjkcljljljlffjljmcBgYicbbbbajljmcBfljmcBgmjkcljljljljmcBgmjkcljmcBgmjkcljljljmcBgmjkcljljmcBgmjkcljljljmcBgmjkcljljljmcBtjljljmcBgYimcBjljljmcBgmjkcljmcBgmjkcljljljmcBgmjkcljmcBgmjkcljljmcBjljljmjmcBgmjkcljljljmcBjljljmcBgYimcBjljljmcBgmjkcljljljmcBgmjkcljljljljmjmcBgmjkcljljmcBgmjkcljljmcBgmjkcljljmcBgjljljmcBgmjkcjljljmcBjjmcBgmjkcljljljmcBgYimcBjljljmcBgmjkcljljljljmcBgmjkcljljljmcBtjljmjmcBgmjkclYimcBjljmcBjljmjmcBgmjkcljmcBjljljmjmcBgmjkcljljmcBgmjkcljljljmcBgmjkcljljljmcBgmjkcljljljmcBgmjkcljljljmcBjljljmcBgjljljljmcBgYimfljmcBjljjcbbbgmjmcljljmcBgmjkcljljljjcbbbbbbgmjkcljljljmcBgmjkcljljljjcbbbbbgmjkcljljljljljmcBfjljmcBgmjkcljmcBgmjkcljljljjmcBgYicbbbbbajljljmcBgmjkcljljljmjmcBgmjkcljmcBgmjkcljljljljmcBgmjkcljljljmcBgmjkcljljljmcBgmjkcljljljljmcBgjljljmcBgjljmcBgmjkcljljljmcBgYimcBjljkjmjmcBgmjkcljljmcBgmjkcljmcBgmjkcljljljmcBgmjkcljljmcBgmjkcljmcBgmjkcljljmcBgYimcBjljljljmcBgmjkcljljljmcBtjljmjmcBgmjkcljljljmcBgmjmcljkjljljljmcBgmjkcljljljmcBgmjkcjljljmcBgjljYimcBjljljjcbbgmjkcljljljmcBgmjkcljljljmcBgYimjjcbbbbjljljmcBgmjkcljkjljmjmcBgmjkcljljljmcBgmjkcljljljmcBgmjkcljljljmcBtjljljmcBgmjkcljljmcjmcBgmjkcljljljmcBgYikjmcBjljljjcbgmjkcljljljmjmcBgmjkcljljljmcBgmjkcljljljmcBgmjkcljljljmcBgmjkcljljmcBgmjljljmcBgYimcBjljljljmcBjljljmcBgmjkcljmcBgmjkcljljljmcjljljmcBgmjkcljlfljmcBgmjkcljljljmcBtjljljmjmcBgYimcBjjljljmcBgmjkcljljljmcBgmjkcljljljmcBgmjkcljmcBgmjkcljljljmcBgmjkcljljljjcbbbbgYimcBjljjcbbbbgmjkcljljljmcBgmjkcljmcBgmjkcljmcBgmjkcljljljljmcBgmjkcljljmcBgmjkcljljljljmcBgYimcBjljljmcBgmjkcljmcBgmjkcljljmcBgmjmcljljljmjmcBgmjkcljljljmcBgmjkcljljmcBgmjkcljljljjcbbbbbgmjkcjljljjcbbbbbgYimcBjjljmcBgmjkcljljljmcBgmjkcljljmcBjljmjmcBgmjkcljljljmcBgmjkcljljljmcBgmjkcljljmcBgjljljmcBgYimcBjljljmcBgmjkcljljljmjmcBgmjkcljmcBgmjkcljljljmcBgmjkcljljljmcBgmjkcljljmcBfYimcBjljljmcBgmjkcljljljmcBgmjkcljljljmcBgmjkcljljljmcBgmjkcljljljmcBgjljljmcBgmjkcljmcBgmjkcljljmcBfjljljmcBjljmcBgYicbbbbbbbjljljjcbbbgmjkcljmcBgmjkcljljljjcbbbbbgmjkcljljljmcBjjljmcBgmjkcljljljmcBgmjkcljljmjmcBgmjkcljljljmcBgjljljmYimjmcBjljljjcbbgmjkcjljljmcBgmjkcljmcBgmjkcljljljmcBYimcBjljljmcBgmjkcljljmcBgmjkcljlfljmcBgmjkcljljljmcBgmjkcljljljmcBgjljljmcBgmjkcljlfmcBgmjkcljljljmcBgjljljmcBgmjkcljljljmcBgYimcBjljmcBjlfljmcBgmjkcljljljkjmcBgmjkcljmcBjljmjmcBgmjkcljljljmcBjljmcBgmjkcljljljbcbbbbgYimcBjljjcbbbjljmjmcBgmjkcljljmcBgmjkcljljljmcBgmjkcljljljmcBgmjkcljljmcBgmjkcljljljljljjcbbbbbgmjkcljljjljljmcBgjljljmcBgYimcBjljljmcBjljmjmcBgmjkcljljmjmcBgmjkcljmcBgmjkcljljljjcbbbbbgmjkcljljljmcBgmjkcljljljmcBgmjkcljljljmcBgjljYimcBjljmcBgmjkcljljljmjmcBgmjkcljljljmcBgmjkcljljljmcBgmjkcljljlYimcBfljmcBgmjkcljljljmcBgmjkcljljljmcBgmjkcljljljbcbbbgmjkcljljlfmjmcBgmYimcBjljljmcBgmjkcljljljljmcBcjljmjmcBgmjkcljljmcBgmYimcBljljmjmcBgmjkcljljljmcBgmjkcljljljmcBgmjkcljljljmcBgmjkcljljljjcbbbbbgmjkcljljljmcBgmjkcljljljmcBgmjkjljljmcBgYimcBjljljmcBgjljmjmcBgmjkcljmcBgmjkcljmcBgmjkcljljmcBgmjkcljljljmcBgmjkcljljljmcBgmjkcljljljjcbbbgjljljljmcBgYimcBjljljjcbbgmjkcljljljmjmcBgmjkcljljmcBgmjkclYimcBjjljljmcBgmjkcljljljmcBgmjkcljmcBgmjkcljljljjcbbbbbgmjkcljlfljmcBfjljljmcBgmjkcljlfljmcBgmjkcljmcBgjljmcBgYfflfljmcjljljmcBjljjcbbgmjkcljmcBgmjkcjljljjcbbbbgmjkcljljljjcbbbgmjkcljljmjmcBgmjkcljljljmcBgjljljmcBgYimcBjljljmcBgmjkcljljljljmcBgmjkcljljmcBgmjkcljljljmcBgmjkcljljljmcBgmjkcjljljmcBgYimcBjljmcBgmjkcljljmcBgmjkcljljljmcBjljmcBgmjkcljljljmcBgmjkcljljljmjljmcBgmjkcljljljmcBYicbbbbbbljljmcBgmjkcljmcBgmjkcljljmcBgmjkcljljljmcBgmjkcljljljmcBgmjkcljljljmcBgmjkcljljljjcbbbbbgmjkcljljmcBjljmcBgYimcBjljjcbbbbbjljljmjmcBgmjkcljljmcBgmjkcljmcBgmjkcljljljmcBgmjkcljljljmcBgmjkcljljljmcBgmjkcljljljjcbbbbYicbbbbbbjjljmcBgmjkcljljmcBgmjkcljljmjmcBgmjkcljljljmcBjjljljmjmcBgmjkcljljljmjmcBgmYimcBfljmcBgmjkcljljmjmcBgmjkcljmcBYicbbbbbbjjljmcBgmjkcljljmcBgmjkcljmcBgmjkcljljljmcBgjljljmcBgmjljljljmcBgmjkcljljljmcBjljmcBjljmjmcjljljmcBYimcBjljljmcBgmjkcljljljmcBjljmjmcBgmjkcljljljmcBgmjkcljljljmYimcBjljmcBgmjkcljljljmjmcBgmjkcljmcBgmjkcljljljmcBtjljljmcBgmjkcljljljmcBjljmjmcBgmjkcljljljmcBgYimjjljmcBjjljljmcBgmjkcljljljmjmcBfjljmclgmjkcljmcBgjljljmcBgmjkcljljljmcBgYimcjljmcBfjljmjmcBgmjkcljljljmcBjljmjmcBgmjkcljmcBgmjkcljljljmcBgmjkcljljljmcBgljmcBgjljljmcBgmjkcljljljmcBgmjkcYimcBjljljmcBgmjkcljljljmcBgmjkcljljljmcBgmjkcljljmcBgmjkcljljljmcBgmjkcljljljjjmcBtYicbbbbbbjjljmcBjljljmjmcBgmjkcljljmcBgmjkcljljljmcBgmjkcljmcBgmjkcljljljmcBgmjkcljljljmcBjljljmcBjmcBgjljljmjkcljljljmcBgmjkcljljljjcbbYimcBgmjljljmcBgmjkcljmcBgmjkcljljmcBgmjkcljljljmcBcjljljmcBgmjkcljljljmcBgYimcBjljljmcBgmjkcljljljljjcbbbbgmjkcljljmcBjljljljmjmcBgmjkcljljmcBgmjkcljmcBgmjljmjmcBgmjkcljljmcBgYimcBjljmcBgmjkcljljljmcBgmjkcljljmcBgmjkcljljlYimcBljmcBfljmcBgmjkcljmcBgmjkcljljljmcBgmjkcljljmcBgmjkcljljljmcBgjljljYimcBjljljmcBgmjkcljmcBgmjkcljljjmjjmcBgmjkcljljljmcBgmjkcljljljmcBgmjkcYimcBjljljmcBgmjkcljljljljmjkcljljljmcBgmjkcljljljmcBgmjkcljmcBgmjkcljljljljljmcBgmjkcljljYimcBjljljmcBgmjkcljljmcBgmjkcljmcBgmjkcljljljmcBgmYimcBjljljjcbbbgmjkcljljljlfmjmcBgmjkcljmcBgmjkcljmcBgmjkcljljljmcBgmjkcljljljmcBgmjkcljljljmjmcBgjljljmcBgYimcBjljljmcBgmjkcljkjljljmcBgmjkcljljmjmcBgmjkcljljjmcBgmjkcljljljmcBYimcBjljmcBgmjkcljljmjmcBgmjkcljljljmjmcBgmjkcljljljmcBgmjkcljljmcBgmjkcljljljmcBjljmcBgmjkcljljljmcBgYimcBjljljljmcBgmjkcljljmjmcBjljmjmcBgmjkcljljljmcBgmjkcljljljmcBgmjkcljmcBgmjkcljljljjcbbbbbgjljljmcBgYimcBjljljmcBgmjkcljljljljmcBjljmjmcBgmjkcljmcBgmjkcljljljmcBgmjkcljljljmcBgmjkcljljljmYimcBjljljmcBjjljmcBgmjkcljljmjmcBjljljmcBgmjkcljljmcBgmjkcljljmcBgjljljmjmcYimcBjljmcBjgmjkcljmcBgmjkcljmcBjlfljljmcBgmjkcljljmcBgmjkcljljmcBgYimcBjljljmcBgmjkcljljljlfmcBgmjkcYimcBjljmcBgmjkcljljljmcBgmjkcljljmcBgmjkcljljljmcBgYimcBjjljmcBgmjkcljljljmcBgmjkcljljmcBgmjkcljljljmcBgmjkcljljljmcBgmjYjlffljmcBjfYimcBjljljmcBgmjkcljmcBfjljmjmcBgmjmclYimcBjljljmcBgmjkcljljljmjmcBgmjkcljlfljmcBgmjkcljljjjmcBgYimcBjljmcBgmjkcljmcBgmjkcljljljmjmcBgYimcBjljljmcBgmjmcBgmjkcljljmcBgmjkcljljmcBgjljljmcBYimjmcBjljjmcBgmjkcljljljmcBgmjkcljljljmcBtjljljmjmcBgYimcBjlYimcBjaYimcljmcBjljljmcBgmjkcljljljmcBgmjkcljmcBgmjkcljljljmcBljjljljmcBgYimcBjljljmcBjljmcjljljmcBjljmjmcBgmjkcljljljmcBgYimcBjljljmcBgmjkcljljljmcBjljmjmcBgmjkcljljljmcBgmjkcljmcBjljlYimcBjljljmcBgmjkcljljljlfmjmcBgmjkcljljljmcBgmjkcljljljljmcBgYimcBjljljmcBgmjkcljljljmcBgmjkcljljmcBgmjkcljjljljlfjmcBgmjkcljljljljmcBgmjkcjljljlYimcBjljmcBgmjkcjljljmcBjljmcBgmjkcljmcBgmjkcljljljmcBgYimcBjljljmcBgmjkcljljljmcBgmjkcljmcBjljljljljmcBgmjkcljljmcBgmjmcBjjljljmcBgYimcBjljmcBgmjkcljmjmcBgmjkcljljmcBgmjkcljmcBgmjkcljljljmcBgmjkcljljljmcBgYimcBjljmcBgmjmcjljmjmcBgmjkcljmcBgmjkcljljljmcBgmjkcljmcBgmjkcljljljmcBgmjkclYimcBjljljmcBgmjkcljljkjmcBgmjkcljmcBgmjkcljljljmcBgmjkcljljmcBgmjkcljlfljmcBgYimcBjjlljljljljmcBjljmjmcBgmjkcljljmcBgmjkcljljljmcBgmjkcljljljmcBgmjljljljmcBjljmcBgYimcBjljljmcBgmjkcljlfljmcBgmjkcljljljmcBjljljmjmcBgmjkcljljljmcBgmjkcllfljmcBgmjljljmjmcBgYimcBjljljmjmcBgmjYimcBjljljmcBgmjkcljmcBjljmjmcBgmjkcljljljmcBgmjkcljljljmcBgmjkcljljljmYimcBjljljmcBgmjmcBgmjkcljljmcBgmjkcljmcBgmjmcBgmjkcljljljmcBgmjkcljljljmcBljljmjmcBgmjkcljljljmcBgYimjmcBjljljmcBgmjkcljmcBgmjkcljljljmjmcBgmjkcljljljmcBgmjkcljljljmcjmcBgmjkcljljljmcBjljljmcBgYimcBjljmcBgmjkcljljljmjmcBgmjkcljljljmcBfjmcBgmjkcljljmcBgmjkcYimcBjljljmcBgmjkcljljmjmcBgmjkcljljljmcBgmjlYimcBjljljmcBgmjkcljljljmjmcBgmjkcljljljmcBgmjkcljmcBgmjkcljljljmcBgjljljmcBgYikjmcBmjljljmcBjljmjmcBgmjkcljljljmcBgmjkcljljmcBgmjkcljljljmjmcBgmjkcljljmcBgYimcBjljjljmcBgmjkcljljljmcBgmjYimcBjljljmcBgmjkcljljmcBgmjkcljljljmcBgmjkcljljmcBfjljljmcBgmjkcljmcBgmjkcljmcBgjljmjmcBtjjljljmcBfYimcBgmjkcljljljmcBYimcBjljljmjmcBgmjkcljljljmcBgmjkcljmcBgmjkcljljmcBgmjkcljmcBgmjkcljljljmcBjljljmcBgjljljmcBgmjkcljljljmcBgYimcBjljljmcBgmjkcljljljmjmcBgmjkcljljmcBgmjkcljmcBgmjkcljljljmcBgjljmcBgjljljmcBgmjkcljljljmcBgYimcBjljljmcBgmjkcljljljljmcBgmjkcljljljmcBjljljmjmcBgmjkcljmcBgmjkcljljmcBgYimcBjljmcBgmjkcljmcBgmjkcljljljmcBgmjkcljljljmcBgmjkcljljmcBjljljmcBgmjkcljljjmcBgmjkcljljljmcBgYimcBjljljmcBgmjkcljljljmcBgmjkcljmcBgmjkcljljljmcBgmjkcljljljmjmcBjjljjmcBgmjkcljljljmcBgYimcBjljljmcBgmjkcljljmjmcjmcBgYimcBjljljmjmcBgmjkcljljljmcBgmjkcljljljmcBgmjkcljljljmcBgmjkcljljljmcBjljljmcBgmjkcljljljmcBgmjkclYimcBjljmcBgmjkcljljljjmcBgmjkcljmcBgmjkcljljljmcBYimcBjljmcBjljmjmcBgmjkcljljljmcBgmjkcljljljmcBjljmjmcBgmjkcljljljmcBgmjkcljljljljmcBgmjkcljljljmcBfjljljmcBjYimcBjljljmjmcBgmjkcljljljljmcBgmjkcljmcBgmjkcljljljmclgmjkcljljljmcBgmjkcljljljmcBgmjkcljljljljmcBYimcBjljljljmjmcBgmjkcljljljmcBgmjkcljljljmcBYimcBjljljmcBgmjkcljljljmcBtjljmjmcBgmjkcljljljmcBgmjkcljljljmcBgmjkcljljljYimcBjljjljmcBjljmjmcBgmjkcljljljmcBgmjkcljmcBgmjkcljmcBgmjkcljljljmcBYimcBjljjljmcBgmjkcljljljmcBgmjkcljmcBgmjkcljljljljmcBgmjkcljljljmcBjjljmcBgYimcBjljljmcBgmjkcljljljljmjmcBgmjkcljljljmcBgmjkcljljljmcBjljljmjmcBgmjkcljmcBgmjkcljljYimjmcBjljmcBgmjkcljljljmcBgmjkcljmcBgmjkcljljljmcBjljljmcBgmjkcljljljmcBgmjkcljljljmcBgjljmjmcBgjjljljmcBgYimcBjljljmcBgmjkcljmcBfjljmjmcBgmjkcljljljmcBgmjkcljljljmcBjgmjkcljYimcBjljljmcBgmjkcljljljmcBjljmjmcBgmjkcljmcBgmjkcljljljmcBYimcBjljmcBgmjkcljljljmcBgmjkcljljmcBgmjkcljljljljmcBgYlflfljmcBjljljljmcBgmjkcljmcBgYimcBjljmjmcBgmjkcljljljmjmcBgmjkcljljmcBgmjkcljljljmcBgmjkcljljljmcBgmjkcljljljmcBgjljmcBgmjkclfjljljmcBgmjkcljmcjljljmcBgmjkcljljljmcBgYimcBjljmjmcBgmjkcljljljmjmcBgmjkcljljljmcBgmjkcljljljmcBgmjkcljljljmcBgmjkcljljmcBgmjkcljlfljmcBgmjkcljljljljmcBgmjkcljljljmcBYimcjljmcBjljljmcBgmjkcljmcBgmjkcljljljmcBjljjjmcBgmjkcljljljljmcBgmjkcljljljmcBgmjkcljljljmcBgmjkcljljljmcBgmjkcljljmcBgmjkcjljljmcBgYimjljmcBjljmjkcjlfljmcBgmjkcljljljmcBgmjkcljljljmjmcBlYimcBjljljmcBgmjkcljljmjmcBgmjkcljljljmcBgmjkcljljljmcBgmjkcljljljmcBgmjkcljljmcBgmjkcljljljmcBgmjkcljljljljmcBgjljljjjmcBgmjkclYimcBjljljmcBgmjkcljljljmjmcBgmjkcljljljmcBfjljmjmcBgmjkcljljlmcBgmjkclYimcBlfljmcBgmjkcljljljmcBgmjkcljljmcBgmjkcljljmcBgmjkcljljljmjmcBjjljljmcBgmjkcljljljmcBgmjkcljljljmcBjYimcBjljljmcBgmjkcljmcBgmjkcljmcBgmjkcljmcjljljmjmcBgmjkcljljmcBgmjkcljlfmcBgmjkcljljljmcBgmjkcljljljljmcBgmjljljmcBgYimcBjljljmcBgmjljmjmcBgmjkcljljljmcBgmjkcljljmcBgmjkcljljljmcBgmjkcljljljljmcBgjljljmcBgmjkcljljljmcBgYimcBjljmcBgmjkcljljljmcBgmjkcljljljmcBtjljmcBgmjkcljljljmjmcBgmjkcljljljljmcBgmjkcljljljkjmcBgmjkcljljljljmcBgYimcBjljljmcBgjljljmcBgmjkcljljljmjmcBgmjkcljljljljmcBgmjkcljljljmcBgYimcBjljljmcBgmjkcljljljljmjmcBgmjkcljljljmcBgmjkcljljljmcBgmjkcljljljljljmcBgmjkcljljmcBgmjkcljmcBgmjkclYimcBjljljmcBgmjkcljljljkjmcBgmjkcljljmcBgmjkcljljljmcBgmjkcljljljmcBgmjkcjljljmcBjljmcBgmjkcljljljmcBgmYimcBjjljljmcBgmjkcljljljmcBYimcBjgmjkcljljljkjmcBgmjkcljljljmcBgmjkcljljljljmcBgmjkcljljmcBgmjkcljljljmcBjljmcBgmjkcljljljmcBgmjkcljljljmcBgjljljmcBgmjkclYimcBjljmcBgmjkcljljmjmcBgmjkcljljmcBgmjkcljljljljljmcBgmjkcljljljmcBgmjkcljljljkcjmcBgmjkcljljljmcBjgYimcBfljmcBjljljmcBgmjkcljljljljmcBgmjkcljljljmcBgmjkcljljljmcBgmjkcljljljmcBgmjkcljljljljYimcBjljljmcBgmjkcljljmjmcBgmjkcljljljmcBgmjkcljljljljmcBgmjkcljljljmcBgmjkcljljljmcBjljmcBgmjkcljljmcBgmjkcjljljmcBjjjjmcBgYimcBjljljmcBgmjkcljljljmjmcBgmjkcljljljljmcBgmjkcljljljmcBgmjkcljljljmcBgmjkcljljljmcBgmjkcjljljmcBgjljljmcBgjmcBgmYimcBjljljmcBgmjkcljljljmcBgmjkcljljljmcBtjljkjmcBgmjkcljljljmcBgmjkcljljljmcBgjljljmcBgmjkcljljljmcBgmjkcljljmcBgmjkcljljljljmcBgYimcBjljljmcBgmjkcljljljmcBgmjkcljljljljmcBjljmjmcBgmjkcljljljmcBgmjkcljljljmcBjYimcBjljljmcBgmjkcljljljmjjmcBgmjkcljljlfljmcBjljmjmcBgmjkcljljljmcBgmjkcljljljljmcBgmjkcljljljljljmcBgmjkcljljmcBYimcBjljljmcBgmjkcljljljljmcBgmjkcljljljmcBgmjkcljljljmcBgmjkcljljljmcBgmjkcljljmcBgYimcBjljmjmcBgmjkcljljljljmjmcBgmjkcljmcBgmjkcljljmjkcljljljmcBjljljmjmcBgmjkcljljljmcBgjljljmcBgmjkcljljljmcBgmjkcljljljmcBjjmcBgmjmcBgYimcBjjljmcBgmjkcljljmjmcBgmjkcljljljmcBgmjkcljljmcBgmjkcljlfljmjmcBgmjkcljljljmcBgmjkcljljljjljmcBgYimcBjjljljmcBgmjkcljljljljmcBjlYimcBjljljmcBgmjkcljljljljmcBgmjkcljljljljmcBgmjkcljljljmcBgmjkcljljljmcjkcljljmcBYimcBjljljmcBtjljmjmcBgmjkcljljljljljmjmcBgmjkcljljljmcBjljmjmcBgmjkcljljljmcBgjYimcBjljljmcBgmjkcljljljljmcBgmjkcljljljmcBtjljmjmcBgmjkcljljljmcBgmjkcljljljljmcBgmjkcljljljmcBgjljljljljmcBjYimcBjfljmcBgmjkcljljljmcBgmjkcljljljljmcBgmjkcljljljmcBgjljljmcBgmjkcljljljmcBjljjjmcBgmjkcljmcBfgjljljmcBgmjkcljljljmcBgjljfjmcBgmjkcljljljmcBgYimcBjjjljljmcBgmjkcjljljljmjmcBgmjmcljljljljmcBgmjkcljljljmcBgmjkcljljljmcBgmjkcljljljmcBgjYimcBjjljljmcBgmjkcljljljljmcBgmjkcljljljmcBgmjkcjjljljljmcBgmjkcljmcBgjljljmjmcBgmjkcljljljmcBgYimcBlfljmcBjljljmcBgmjkcljljljmcBjljljmjmcBgmjkcljljmcBgmjkcljmcBgmjkcljljljmcBjljmjmcBgmjkcljljljljmcBgYimcBjlfljmcBgmjkcljljljljmjmcBgmjkcljljljljmcBgmjkcljljljmcBgmjkcljljmcBgmjkcjljljmcBgmjkcljlYimcBjljmcBgmjmcBgmjmcBgmjkcljljljmcBgmjkcljljljmjmcBtjljljmcBjljljmcBgmjkcjljljmcBgmjkcljljljmcBgjljljmcBgmjkcljljljmcBgjljljmcBgmjkcljljmcjjYimcBjljljljmcBgmjkcljlfljljmcBgmjkcljljljmcBgmjkcljljmcBgmjkcljmcBgjjljljmcBgmjkcljljljmcBgmjkcljljljmcBgYimcljmcBjljmcBgmjkcljmcBgmjkcljljljmjmcBgmjkcljljmcBgmjkcljljljmcBYimcBjljljmcBgmjkcljmcBgmjkcjljljmcBgmjkcljljljmjmcBgjljljmcBgmjkcljljljmcBgjljljmcBgYimcBjljljmcBgmjkcljjljljljljmcBgmjkcljmcBgmjkcljljljmcBgmjkcljljljmcBgmjkcljljljmcBgmjkcljljljmcBgmjljljmcBgmjkcjljljmcBgjljmcBjjljljmcBYimjmcBjljljmcBgmjkcljljljjmcBgmjkcljljljmcBgmjkcljmcBgmjkcljljljmcBgmjkcljljljmcBjjljljmjmcBgjljmcBgmjkcljljljmcBgYimcBfljmcBgmjkcljljljmcBjljmjmcBgmjkcljljljljmcBgmjkcljljljmcBgmjkcljljljmcBgmjkcljljljmcBgmjkcljljljmcBgYimcBjljljmcBgmjkcljljljmcBgmjkcljljljmcBjljmjmcBgmjkcljljljmcBgmjkcljlfljmcBtjkcljmcBgYimcBjjljljmcBjljmjmcBgmjkcljljljmcBfjmcBgmjkcljljlfljmcBgmjkjljljmcBgmjkcljljljljmcBgmjkcljljljmcBgmjkcljljljmcBgmjkcljljljmcBgYimcBgmjkcljljljmcBgmjkcljljljmcBgmjkcljljljmcBgmjkcljljljljmcBgmjkcljljljmcBgmjkcljljljljmcBgmjkcljljljmcBgmjkcljljljmcBgmjkcljljljmcBgYimcBjljljmcBjljljkjmcBgmjkcljljljljmcBgmjkcljljljmcBgmjkcljljljmcBgmjkcljljljmcBgmjkcljljmcBgmjkcjljljmcBgjljljmcBgYimcBjljljmcBgmjkcljljljmjmcBgmjkcljljljljmcBgmjkcljljljmcBgmjkcljljljYimcBjljljmcBgmjkcljljmjmcBgmjkcljmcBgYimcBjljljmcBgmjkcljljmjmcBgmjkcljljljmcBgmjkcljljmcBgmjkcljljljmcBjljljmcBgYimjljmcBjljljmcBgmjkcljljmjmcBgmjkcljljljmcBgmjljmjmcBYimcBjljljmcBgmjkcljljljljmcBgmjkcljljljmcBgmjkcljmcBjjljmjmjkcljljljmcBgjljljljmcBgjljljmcBYimcBjljljmcBgmjkcljljmjmcBgmjkcljljljmcBgmjkcljljljmcBjljljmcBgmjkcljljljmcBgmjkcljljljmcBgjljljmcBgjljljmcBgYimcBjljljmcBgmjkcljljljmcBgmjkcljljmcBgmjkcljljljmcBgmjkcljljljmcBgjljljmjmcBgmjkcljmcBgmjkcljljljmcBgjljmcBgmjkcljljljmcBgjljljYimcBjljljmcBjljljmcBgmjkcljmcBgmjkcljljljmcBgmjkcljljljmcBYimcBjljljljmcBgmjkcljYimcBjljljmcBgmjkcljmcBgmjkcljljljmcBYimcBjljljmcBgmjkcljljkjmcBgmjkcljlfljmcBgmjkcljmcBjljljmcBgmjkcljljljmcBgmjkcljljljljmcYflfljmcBjljljmcBgmjkcljlfljmcBjljmjmcBgmjkcljljljljmcBjljljmjmcBgmjkcljljljmcBgmjkcljljljlfmcBgmjkcljljljmcBgYimcBjljmcBgmjkcjljljkjmcBgmjkcljljljljmcBgmjkcljljljmcBtjljljmjmcBgmjkcljljmcBgmjkcljljljljmcBjlYimcBjljljmcBgmjkcljlfmcBgmjkcljljljmcBjljmjmcBgmjkcljljljmcBgmjkcljljljmcBgmjkcljmcljljljmcBgmjkcljljmcBgjljljmcBgYimcBjljmcBgmjkcljljljmcBgmjkcljmcBgmjkcljljljmcBgmjkcljlfljmcBgmjkcljljljljmcBjljljmcBgmjkclYimcBjjljmcBgmjkcljljmcBgmjkcljmcBgmjkcljljljmcBgmjkcljljljljljljljmcBgmjkcljljljmcBgmjkcjljljmcBjljljmcBgjljljmcBgYimcBjljljmcBgmjkcljljljmcBgmjkcljmcBgmjkcljljljmjmcBgmjkcljljljkjmcBgmjkcljljmcBgmjkcljljljmcBgYimcBjljljmcBgmjkcljljljljmcBgmjkcljljljmcBgmjkcljljmcBgmjkcljljljmcBmcjljmjmcBYimcBjljmjmcBgmjkcljljljljmcBgmjkcljljmcBgmjkcljljljmcBgmjkcljljmcBgmjkcljljljmcBgYimcBjljmcBjljljmcBgYimcBjjljljmcBgmjkcljljljmcBgmjkcljljmjmcBgmjkcljljmcBgmjkcljljljmcBgmjkcljljmcBfYimcBjljljmcBgmjkcljmcBgmjkcljljljmcBtjljljljmcBgmjkcljljljmcBgYimcBjljljmcBgmjkcljljmjmcBgmjkcljljljmcBtjljmjmcBgmjkcljljljmcjmcBgmjkcljljljmcBgmjkcljljmcBgmjkcljljljmcBtjljljjmcBgmjkcYimcBjjljmcBgmjkcljljljmjmcBgmjkcljljmcBjljkjmcBgmjkcljljljmcBgmjkcljljljmcBjljljmjmcBgmjkcljljljmcBgmYimcBjljljmcBgmjkcljljljmcBgmjkcljmcBjljmjmcBgmjkcljljljmcBgmjkcljljljljmcBgmjkcljljljmcBgmjkcljmcBgjljljjkjmcBgYimcBjljmcBgmjkcljljljkjmcBgmjkcljljmcBgmjkcljmcBgmjkcljljljmcBjljljmjmcBgmjkcljljljmcBgYimcBjljljmcBgmjkcljljljmcBgmjkcljmcBgmjkcljljljljmcBgmjkcljmcBgmjkcljljljmcBgmjkcljljljmcBgmjkcljljljmcBgmjjljljjkjmcBgYimcBjljljmcBgmjkcljljljmjmcBgmjkcljljljmcBgmjkcljljljmcBgmjkcljljljmcBgmjkcljljmcBgmjkcljljlYimcBjjmcBjjmjmcBgmjkcljljmjmcBgmjYimcBjljmcBgmjkcljljljljmcBgmjkcljmcBgmjkcljljljmcBgmjkcljljljljmcBgmjkcljmcBgmjkcljljljmcBgjYimcBjljmcBgmjmcBgmjkcljljljmcBgmjkcljmcBgmjkcljljljmcBgmjkcljljljmcBgmjkcljljljmcBYimcBjlfljmcBgmjkcljmcBgmjkcljljljljmcBgmjkcljljljmcBgmjkcljljljljmcBgmYimcBjljmcBgmjkcljmcBfjljljljmcBgmjkcljljljmcBgmjkcljljljmcBjljmjmcBgmjkcljljljmcBgmjkcljljljmcBjjljljmcBYimcBjljljmcBgmjkcljljljmcBgmjkcljmcBgmjkcljljljmcBgmjkcljljljljmcBgmjkcljljljmcBgYimcBjljmcBgmjkcljmcBgmjkcljljljmcBjljljmcBjjljljljmcBgmjkcljljljmjmcBgmjkcljljljmcBgmYimcBfljmcBgmjkcljmcBjljmjmcBgmjkcljljljmcBfjljmjmYimcBjljljmcBjljmcBgmjkcljljljmcBgmjkcljmcBgmjkcljlfljmcBjljljmcBjljmjmcBgmjkcljljljmcBgmjkcljljljmcBgmjkcljljljmcBgYimcBjljmcBgmjkcljmcBgmjkcljljmcBgmjkcljmcBgmjkcljljljmcBgmYimcBjljmcBgmjkcljmcBgmjkcljmcBgmjkcljljljmcBgmjkcljljmcBljljljmcjmcBgmjkcljljljmcBgjljljmcBgmjkcljljljljmcBtjYimcBjljmcBgmjkcljljmjmcBgmjkcljljmcBjljljmjmcBgmjkYimcBjljljmcBgmjkcljljljmcBgmjkcljljljljmcBgmjkcljljljmcBgmjkcljmcBgYimcBjljmcBgmjkcljljljmjmcBgmjkcljljljljmcBgmjkcljljljmcBgmjkcljljljmcBgmjkcljljljmcBgmjkcljljljmcBgmjkclYimcBjljljmcBgmjkcljljljmcBgmjkcljljmcBgmjkcljmcBgmjkcljlYimcBjljljmcBgmjkcljljlfljmcBjgmjkcljYimcBjljljmcBgmjkcljljljmjmcBgYimcBjljljmcBgmjkcljljljmcBgmjkcljljljmcBgmjkcljljljmcBgmjkcljmcBgmjkcljljlfljmcBjljmcBgmjkcljljljmcBgYimcBjljmcBgmjkcljljljmcBgmjkcljljljmcBgmjkcljljljmcBgmjkcljljljmcBgmjkcljljljljmcBgmjljljmcBgjljljmcBgmjkcljljljmcBgYimcBjljljmcBgmjkcljljljmcBgmjkcljljljmcBgmjkcljlYimcBjljljmcBgmjkcljljljmcBgmjkcljljljljmcjljljmcBgmjkcljljljmcBgmjkcjljljmcBgmjkcljljljmcBgmjkcljljljljljmcBgYimcBjljljmcBgmjkcljljljmcBgmjkcljljljmcBgmjkcljljljmcBjljljmcBgmjkcljmcBgmjkcjljmcBgjljljmcBgmjkcljljljmcBjljmcBgYimcBjfljmcBjljljmcBgmjkcljljljljmcBgmjkcljljYimcBjljljmcBgmjkcljljljljmcBgmjkcljljljmcBjljmjmcBgmjkcljmcBgmjkcljljljmcBgmjkcljljljmcBgmjkcljljljljmcBgjljljljmcBgmjkcljljljmcjYimcBjljmcBgmjkcljljljmjmcBgmjkcljljljmcBgmjkcljljljmcBgmjkcljljljmcBgmjkcljljljmcBgmjkcljljljmcBgmjkcljljljmcBgjljjYimcBjljmcBgmjkcljljljmcBgmjkcljljljmcBgmjkcljljmcBgmjkcljljljljmcBgmjkcljljljmcBfYimcBjljljmcBgmjkcljljljmjmcBgmjkcljljljmcBjljljljmjmcBgmjkcljljljmcBgmjkcljljljmcBgmjkcljYimcBjljljmcBgmjkcljljljmcBgmjkcljljljmcBgmjkcljljljmcBgmjkcljljljmcBgmjkcljljljmcBjlfljmjmcBgYimcBjljljmcBgmjkcljljmjmcBjljmjmcBgmjkcljljljmcBjljljmcBgmjkcljljljkYimcBjjljmcBgmjkcljljljmcBgmjkcljljljmcBgmjkcljljljmcBtjljmjmcBgmjkcljljjmcBgmjkcljljljmcBgYimcBjljljmcBgmjkcljljljmcBgmjkcljljljmcBgmjmcBgmjkcljljljmcBgmjkcljljljljmcBgmjkcljljljljljmcBgmjkcljljljmcBYimcBjjjljmcBjjljljmcBgmjkcljmcBgmjkcljljljljmcBgmjkcljljljmcBgmjkcljljljmcBgmjkcljlYimcBjljmcBgmjkcljljljmcBgmjkcljljljmcBgmjkcljljljmcBgmjkcljljljmcBgmjkcljljljmcBgYimcBjljmcBgmjkcljljljmcBjljmcBjljljmcBgmjkcljljljmcBgmjkcljljmcBjljmcBgmjkcljmcBgmjkcljljljmcBgmjkcljljljmcBgYimcBjljljmcBgmjkcljljljljmcBgmjkcljkjljmcBgmjljljmjmcBgmjkcljljljmcBljlYimflfljmcBjljmcBjljljljmcBgmjkcljljljmjmcBgmjkcljljljmcBfjljmcBgmjkcljmcBgmjkcljljljmcBgmjkcYimcBjljmcBjljmcBgmjkcljljljmcBgmjkcljljmcBgmjkcljmcBgmjkcljljljmcBjljmjmcBgmjkcljmcBgmjkcljljljljljljmcBYimcBjljljmcBjljljmcBgmjkcljljljmcBgmjkcljmcBgmjkcljljljmcBgmjkcljljmcBgmjkcljljljmcBgYikjjljmcBjljljmjmcBgmjkcljmcBgmjkcljljljmcBgmjkcljljljmcBjljmjmcBgmjkcljljljmcBgmjkcljljljmcBtjljmcBgmYimcBjljljmcBgmjkcljljljmjmcBgmjkcljljljmcBgmjkcljljmcBgmjkcljljljmcBgmjmcBgmjkcljljljmcBYimcBjljljljmcBgmjkcljljljmjmcBgmjkcljljljmcBgmjkcljmcBgmjkcljljljljmcBgjljljmcBgYimcBjljmcBgmjjljmcBgmjkcljljljmcBgmjkcljljljmcBgmjkcljljljljmcBgmjkcljljljmcBgmjmjmcBgmjkcljmcBgmjkcljljljmcBgjljjmcBtjljljmcBYimcBjljljmcBljljmjmcBgmjkcljljljmcBgmjkcljljljYimcBjljljmcBgmjkcljljljmcBYimcBtjljmcBgmjkcljljljmcBgmjkcljljljmcBgmjkcljljljmcBgmjkcljmcBgmjkcljljljmcBgmjkcljljljmcBYimcBgmjkcljljljmcBgmjkcljmcBgmjkcljmcBgmjkcljljljmcBgmjkcljkjmjmcBgmjkcljljljmcBgmYimcBjljljmcBgmjljljmcBgmjkcljljmcBgjmcBgmjkcjljljmcBgmjkcljljljmcBYimcBgmjkcljljkjmcBgmjkcljljljljmcBgmjkcljljljmcBgmjkcljljljljmcBjljmcBgmjkcljmcBgYimcBjljmcBjjljljmcBgmjkcljmcBgmjkcljljljmcBjljmjmclgmjkcljljljmcBjljljmcBgjmcBgmjkcljljmcBgYimcBgmjjljmcBjljmcBgmjkcljmcBgmjkcljljljmcBtjYimcBjljljmcBgmjkcljljljmcBgjljljmcBgmjkcljmcBgmjkljljljmcBgmjkcljljljmcBgmjmcBgmjkcljljljmcBgmjkcljljljmcBgmjkcljljljmcBgmjljmcBgYimcBjljmcBgmjkcljmcBgmjkcljljmcBgmjkcljljljmcBgmjkcljljmjmcBgmjkcljljljmcBgmjkcljljljljmcBgYimcBjfljmcBgYimcBjljmcBgmjkcljmcBgmjkcljmcBgmjkcljljljljmcBgmjkcljljljmcBgmjkcljljljmcBgYimjmcBjljljmcBjljmjmcBgmjkcljljljmcBgmjkcljljljmcBgmjkcljljljmcBjljmcBgmjkcljljmcBgmjkcljljljljmcBgYimcBtjljmcBjljmjmcBgmjkcljljljmcBgmjkcljljljmcBgmjkcljmcBgmjkcljljljmcBYimcBjljljmcBgmjljljmcBgmjkcljljmjmcBgmjkcljmcBgmjkcljljljljmcBgmjkcljljljljmcBgmjkcljljljmcBgmjkcljljljljmcBgmjkcljmcjljljljmcBgYimcBjljmcBgmjkcljljljmjmcBgmjkcljljljmcBgmjkcljmcBgmjkcljljljmcBgmjmcBgmjkcljljljmcBgYimcBjljmcBgmjkcljmcBgmjkcljljljmcBgmjkcljljjmcBgmjkcljljljmcBgmjkcljljljmcBgYimcBgmjmjljmcjljmcBjljmjmcBgmjkcljljljmcBgmjkcljmcjYimcBjljljmcBgmjkcljljljmcBgmjkYimcBjljljmclfljmcBgmjkcljmcBgmjkcljljmcBgmjkcljljljmcBgmjkcljljljmcBjljmjmcBgmjkcljljljmcBgYimcBjljmcBjljljmcBgmjkcljljljljmcBgmjkcljljljmcBjljmjmcBgmjjljljljmcBgmjkcljmcBgjljljljmjmcBgmjkcljmcBgYimcBjljljmcBgmjkcljmcBgmjkcljmcBgmjkcljljljmcBtjljmjmcBgmjkcljljljmcBgmjkcljljljljmcBtjjljljmcBgjljmcBjjljljljmcBYimcBjljmcBjljmcBgmjkcljljljmcBgmjkcljmcBgmjkcljljljmcBgmjkcljljljmcBgmjkcljljljmcBgYimcBjljmcBgmjkcjljljmcBgjljljmcBgmjkcljljljmcBgmjkcljljmcBgmjkcljljljmcBgmjkcljljjljmcBgmjkcljljljmcBgmjkcljljljljmjmcBgYimcBjfljmcBfjljljmcBgmjkcljljljkjmcBgmjkcljljljmcBgmjkcljljljmcBgmjkcljljljjmcBtjYimcBfljmcBjljljmjmcBgmjkcljmcBYimcBjljljmcBgmjkcljljljmcBjljmjkcljmcBgmjkcljljljljmjmcBgmjkcljljljmcjYimcBjljljmcBgmjkcljljmcBgmjkcljljljmcBgmjkcljljljmcBgmjkcljljljljmcBgmjkcljljljmcBgmjkcljljljmcBgjljljmcBgYimjmcBjljjmcBjljljmcBgmjkcljljmcBgmjkcljmcBgmjYimcBjljljmcBgmjkcljljljmcBgmjkcljljljmcBgmjkcljmcBgmjkcljljljljljmcBgmjljmcBgjljljmcBgjljljmcBgmjkcljljmcBgjljmcBgYimcBjljmcBgmjkcjljljmcBgmjkcljljljmcBgmjkcjljljmcBgmjkcljljljmcBgmjkcljlfljmcBgmjkcljljljljmcBgmjkcljljljmcBgmjkcljljljmcjjljljmcBjljmcBYimcBjljljmcBgmjkcljljljljmcBgmjkcljljljljmcBgmjkcljljljmcBgmjljmcBgmjkcljljljljYimcBjljljljmcBgmjkcljljmcBgmjkcljljljmcBgmjkcljmcBgmjkcljljljmcBjljmcBgmjmcBgmjkcljljmcBgYimjljmcBjljljmcBgmjkcljmcBgmjkcljljljmcBgmjkcljmcBgmjkcljljljmcBgmjljljmcBgjljljmcBgmjkcljljljmcBgjljkjkjmcjmcBgYimcBgmjkcljljljmcBgmjkcljmcBgmjkcljljljljmcBjjljljmcBgmjkcljljjmcBfjljljmcBgmjkcljljljmcBgmY"

    new_seqs = make_states_of_syl_chi2(seqs)
