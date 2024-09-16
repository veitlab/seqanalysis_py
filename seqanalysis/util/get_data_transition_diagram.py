import re
import string

import numpy as np
import scipy.io as sio
from IPython import embed

from seqanalysis.util.logging import config_logging

log = config_logging()


def get_labels(mat_list, notes, intro_replacement):
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
        try:
            mat = sio.loadmat(str(matidx))

            # the following if statment are because nils labels are not consistently only one type of data,
            # following comments are the data type in matlab
            if type(mat['labels']) == np.ndarray and type(mat['labels'][0]) == np.str_: # if 1x37 char
                labels = mat['labels']
                labels = '_' + ''.join(str(x) for x in labels)
            elif type(mat['labels'][0]) == np.ndarray and type(mat['labels'][0][0]) == np.int32: # if int
                labels = mat['labels'][0]
                labels = '_' + ''.join(str(x) for x in labels)
            elif type(mat['labels'][0]) == np.str_: # if 37 char
                labels = mat['labels'][0]
                labels = '_' + labels

        except OSError:
            log.error(f"File not found: {matidx}")
            continue

        # log.debug(f"Processing file: {matidx}")

        if len(notes) > 0:
            try:
                labels = replace_intro_notes(labels, notes, intro_replacement)
                log.debug(f"Intro notes replaced in file: {matidx}")
            except ValueError:
                log.error(f"Intro notes not found in file: {matidx} ")

        seqs.append(labels)
    seqs = np.array(seqs)

    return seqs


def replace_intro_notes(s, intro_notes, replacement):
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
        if str(intro_note) in s:
            unique_labels.remove(str(intro_note))

    motiv_start = []
    for unique_label in unique_labels:
        motiv_start.append(s.find(unique_label))

    temp = list(s)
    for i in range(1, np.min(motiv_start)):
        temp[i] = replacement[0]

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
        ch.append((str(chunk), asci_letters[i]))
        s = re.sub(str(chunk), asci_letters[i], s)
    return s, ch


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

