import glob
import numpy as np
import scipy.io as sio
import re


def get_labels(mat_list, notes):
    seqs = []
    for matidx in mat_list:
        mat = sio.loadmat(matidx)

        labels = mat['labels'][0]
        labels = 'S' + labels + 'E'

        if len(notes) > 0:
            labels = replace_intro_notes(labels, notes)

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


def replace_chunks(s, chunks):
    # ToDo: was wenn diechunks mit dem selben buchstaben beginnen?

    for i in range(len(chunks)):
        s = re.sub(chunks[i], chunks[i][0].upper(), s)

    return s


def get_syl_dur():
    mat_list = glob.glob('*cbin.not.mat')
    durs = []
    for matidx in mat_list:
        mat = sio.loadmat(matidx)

        dur = mat['offsets'] - mat['onsets']
        durs.append(dur)

    durs = np.array(durs)

    return durs


def get_gap_dur():
    mat_list = glob.glob('*cbin.not.mat')
    durs = []
    for matidx in mat_list:
        mat = sio.loadmat(matidx)

        dur = mat['onsets'][1:] - mat['offsets'][:-1]
        durs.append(dur)

    durs = np.array(durs)

    return durs


def get_bouts(seqs, bout_string):
    bouts = ''
    noise = ''
    for seqsidx in range(len(seqs)):
        if seqs[seqsidx].find(bout_string) >= 0:
            bouts = bouts + seqs[seqsidx]
        elif seqs[seqsidx].find(bout_string) < 0:
            noise = noise + seqs[seqsidx]

    # bouts = np.array(bouts)
    # noise = np.array(noise)
    return bouts, noise


def get_transition_matrix(bout, unique_labels):
    transM = np.zeros((len(unique_labels), len(unique_labels)))
    transM_prob = np.zeros((len(unique_labels), len(unique_labels)))

    alphabet = {letter: index for index, letter in enumerate(unique_labels)}
    numbers = [alphabet[character] for character in bout if character in alphabet]

    for idx in range(len(numbers) - 1):
        transM[numbers[idx], numbers[idx + 1]] += 1
        transM_prob[numbers[idx], numbers[idx + 1]] += 1

    transM_prob = (transM_prob.T / np.sum(transM, axis=1)).T
    transM = transM.astype(int)

    return transM, transM_prob


def get_node_positions(source_target_list):
    xpos = np.array([int(string[0]) for string in source_target_list])
    ypos = np.zeros(len(xpos))
    for i in range(len(np.unique(xpos))):
        ypos[xpos == i] = [*range(len(xpos[xpos == i]))]

    pos = np.column_stack((xpos, np.array(ypos)))

    return pos


def get_node_matrix(matrix, edge_thres):
    '''
    calculates the transition probabilities int percent values and sets every edge below the edge_threshold to zero
    :param
            matrix: numpy array; matrix of the transition probabilities
            edge_thres: int; threshold where edges below go to zero
    :return:
            matrix: numpy array; same size as input matrix, with probabilities in percent between 0-100
    '''

    matrix = np.around(matrix, 2) * 100
    matrix = matrix.astype(int)
    matrix[matrix < edge_thres] = 0

    return matrix

def findMaximalNonBranchingPaths(Graph):
    paths = []
    nodes1in1out = set()  # 1-in-1-out nodes
    nExplored = set()  # 1-in-1-out nodes which were explored
    for v in Graph.adj.keys():
        if not (1 == Graph.in_degree[v] and 1 == Graph.out_degree[v]):
            if Graph.out_degree[v] > 0:
                for w in Graph.adj[v].keys():
                    nbPath = [v, w]  # NonBrachingPath

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
