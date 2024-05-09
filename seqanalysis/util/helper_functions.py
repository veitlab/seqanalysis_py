import glob
import re
import string

import numpy as np
import scipy.io as sio
from IPython import embed

from seqanalysis.util.logging import config_logging

log = config_logging()


def get_data(path, intro_notes, bout_chunk):
    file_list = glob.glob(path)

    seqs = get_labels(file_list, intro_notes)
    bouts, _ = get_bouts(seqs, bout_chunk)

    return bouts


def get_labels(mat_list, notes):
    """
    Extracts sequence labels from a list of .mat files.

    Parameters:
    - mat_list (list): List of .mat file paths.
    - notes (str): Additional notes to replace in the labels.

    Returns:
    - seqs (numpy array): Array of sequence labels.
    """
    seqs = []
    for matidx in mat_list:
        mat = sio.loadmat(matidx)
        log.debug(f"Processing file: {matidx}")

        labels = mat["labels"][0]
        labels = "_" + labels

        if len(notes) > 0:
            try:
                labels = replace_intro_notes(labels, notes)
                log.debug(f"Intro notes replaced in file: {matidx}")
            except ValueError:
                log.error(f"Intro notes not found in file: {matidx} ")

        seqs.append(labels)
    seqs = np.array(seqs)

    return seqs


def replace_intro_notes(s, intro_notes):
    """
    Replaces introductory notes in a sequence.

    Parameters:
    - s (str): Input sequence.
    - intro_notes (list): List of introductory notes to be replaced.

    Returns:
    - s (str): Sequence with replaced introductory notes.
    """
    unique_labels = sorted(list(set(s)))
    for i, intro_note in enumerate(intro_notes):
        if intro_note in s:
            unique_labels.remove(intro_note)

    motiv_start = []
    for unique_label in unique_labels:
        motiv_start.append(s.find(unique_label))

    temp = list(s)
    for i in range(1, np.min(motiv_start)):
        temp[i] = "x"

    s = "".join(temp)

    return s


def replace_chunks(s, chunks):
    """
    Replaces chunks in a sequence.

    Parameters:
    - s (str): Input sequence.
    - chunks (list): List of chunks to be replaced.

    Returns:
    - s (str): Sequence with replaced chunks.
    """
    asci_letters = list(string.ascii_uppercase)
    ch = []
    for i, chunk in enumerate(chunks):
        log.info(f"Replacing chunk: {chunk}, with {asci_letters[i]}")
        ch.append((chunk, asci_letters[i]))
        s = re.sub(chunk, asci_letters[i], s)
    return s, ch


def get_syl_dur():
    """
    Retrieves syllable durations from a list of .mat files.

    Returns:
    - durs (numpy array): Array of syllable durations.
    """
    mat_list = glob.glob("*cbin.not.mat")
    durs = []
    for matidx in mat_list:
        mat = sio.loadmat(matidx)

        dur = mat["offsets"] - mat["onsets"]
        durs.append(dur)

    durs = np.array(durs)

    return durs


def get_gap_dur():
    """
    Retrieves gap durations from a list of .mat files.

    Returns:
    - durs (numpy array): Array of gap durations.
    """
    mat_list = glob.glob("*cbin.not.mat")
    durs = []
    for matidx in mat_list:
        mat = sio.loadmat(matidx)

        dur = mat["onsets"][1:] - mat["offsets"][:-1]
        durs.append(dur)

    durs = np.array(durs)

    return durs


def get_bouts(seqs, bout_string):
    """
    Extracts bouts and noise from a list of sequences.

    Parameters:
    - seqs (numpy array): Array of sequences.
    - bout_string (str): String to identify bouts.

    Returns:
    - bouts (str): Concatenated bouts.
    - noise (str): Concatenated non-bout sequences.
    """
    bouts = ""
    noise = ""
    for seqsidx in range(len(seqs)):
        if seqs[seqsidx].find(bout_string) >= 0:
            bouts = bouts + seqs[seqsidx]
        elif seqs[seqsidx].find(bout_string) < 0:
            log.debug(f"Sequence {seqsidx} is not a bout")
            noise = noise + seqs[seqsidx]

    return bouts, noise


def get_transition_matrix(bout, unique_labels):
    """
    Computes transition matrix and probability matrix for a given bout.

    Parameters:
    - bout (str): Bout sequence.
    - unique_labels (list): List of unique labels.

    Returns:
    - transM (numpy array): Transition matrix.
    - transM_prob (numpy array): Transition probability matrix.
    """

    transM = np.zeros((len(unique_labels), len(unique_labels)))
    transM_prob = np.zeros((len(unique_labels), len(unique_labels)))

    alphabet = {letter: index for index, letter in enumerate(unique_labels)}
    numbers = [alphabet[character] for character in bout if character in alphabet]

    for idx in range(len(numbers) - 1):
        transM[numbers[idx], numbers[idx + 1]] += 1

    # Normalize transition matrix
    transM_prob = (transM.T / np.sum(transM, axis=1)).T
    transM = transM.astype(int)

    return transM, transM_prob


def get_transition_matrix_befor_following_syl(bout, unique_lables):
    """
    Computes transition matrix and probability matrix for preceding and following syllables in a given bout.

    Parameters:
    - bout (str): Bout sequence.
    - unique_lables (list): List of unique labels.

    Returns:
    - transM_bsf (numpy array): Transition matrix.
    - transM_prob_bsf (numpy array): Transition probability matrix.
    """
    transM_bsf = np.zeros((len(unique_lables), len(unique_lables), len(unique_lables)))
    transM_prob_bsf = np.zeros(
        (len(unique_lables), len(unique_lables), len(unique_lables))
    )
    alphabet = {letter: index for index, letter in enumerate(unique_lables)}
    numbers = [alphabet[character] for character in bout if character in alphabet]

    for idx in range(1, len(numbers) - 1, 1):
        transM_bsf[numbers[idx - 1], numbers[idx], numbers[idx + 1]] += 1
        transM_prob_bsf[numbers[idx - 1], numbers[idx], numbers[idx + 1]] += 1

    # BUG: Transition Maxtrix is not updated
    # transM_prob_bsf = (transM_prob_bsf.T / np.sum(transM_bsf, axis=1)).T
    transM_bsf = transM_bsf.astype(int)

    return transM_bsf, transM_prob_bsf


def get_node_positions(source_target_list):
    """
    Computes node positions based on a source-target list.

    Parameters:
    - source_target_list (list): List of source-target pairs.

    Returns:
    - pos (numpy array): Array of node positions.
    """
    xpos = np.array([int(string[0]) for string in source_target_list])
    ypos = np.zeros(len(xpos))
    for i in range(len(np.unique(xpos))):
        ypos[xpos == i] = [*range(len(xpos[xpos == i]))]

    pos = np.column_stack((xpos, np.array(ypos)))

    return pos


def get_node_matrix(matrix, edge_thres):
    """
    Calculates transition probabilities in percent and sets edges below a threshold to zero.

    Parameters:
    - matrix (numpy array): Transition probability matrix.
    - edge_thres (int): Threshold for edges to be set to zero.

    Returns:
    - matrix (numpy array): Updated matrix with probabilities in percent.
    """
    matrix = np.around(matrix, 2) * 100
    matrix = matrix.astype(int)
    matrix[matrix < edge_thres] = 0

    return matrix


def findMaximalNonBranchingPaths(Graph):
    """
    Finds maximal non-branching paths in a graph.

    Parameters:
    - Graph (networkx.DiGraph): Directed graph.

    Returns:
    - paths (list): List of maximal non-branching paths.
    """
    paths = []
    nodes1in1out = set()  # 1-in-1-out nodes
    nExplored = set()  # 1-in-1-out nodes which were explored
    for v in Graph.adj.keys():
        if not (1 == Graph.in_degree[v] and 1 == Graph.out_degree[v]):
            if Graph.out_degree[v] > 0:
                for w in Graph.adj[v].keys():
                    nbPath = [v, w]  # NonBranchingPath

        else:
            nodes1in1out.add(v)

    for v in nodes1in1out:
        if v not in nExplored:
            w = v
            nbPath = []
            while w in nodes1in1out:
                nbPath.append(w)
                if w == v and len(nbPath) > 1:
                    paths.append(nbPath)
                    for node in nbPath:
                        nExplored.add(node)
                    break
                w = Graph.adj[w][0]
    return paths
