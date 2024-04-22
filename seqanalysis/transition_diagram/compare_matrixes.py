import re
import sys
import glob
import yaml
import numpy as np
import seaborn as sns
import matplotlib.pyplot as plt
import seqanalysis.util.plot_functions as pf
import seqanalysis.util.helper_functions as hf

from IPython import embed


def get_matrix(bird):
    bouts = bird["data"]["chunk_bouts"]
    tm, _ = hf.get_transition_matrix(bouts, bird["labels"]["unique_labels"])

    # Filter out nodes with low occurrence
    k = np.where(np.sum(tm, axis=0) / np.sum(tm) * 100 <= 0.01)
    tmd = np.delete(tm, k, axis=1)
    tmd = np.delete(tmd, k, axis=0)
    print(np.delete(bird["labels"]["unique_labels"], k))

    # Normalize transition matrix and create node matrix
    tmpd = (tmd.T / np.sum(tmd, axis=1)).T
    tmpd = hf.get_node_matrix(tmpd, bird["constants"]["edge_threshold"])

    return tmd, tmpd


def load_bird_data(yaml_file):
    with open(yaml_file) as f:
        bird = yaml.load(f, Loader=yaml.FullLoader)
    return bird


def main(yaml_file1, yaml_file2):
    bird1 = load_bird_data(yaml_file1)
    bird2 = load_bird_data(yaml_file2)

    matrix1, matrix1p = get_matrix(bird1)
    matrix2, matrix2p = get_matrix(bird2)

    embed()

    # labels = ['E', 'I', 'B', 'C', 'D', 'e', 'a', 'F', 'K']
    labels = ["E", "I", "B", "C", "D", "e", "a", "F"]
    matrix = matrix2p - matrix1p
    matrix = matrix[:-1, :-1]

    fig, ax = plt.subplots()
    colormap1 = sns.diverging_palette(160, 40, as_cmap=True)
    hm = sns.heatmap(
        matrix,
        ax=ax,
        annot=True,
        vmin=-100,
        vmax=100,
        fmt="d",
        cmap=colormap1,
        xticklabels=labels,
        yticklabels=labels,
    )
    ax.set_yticklabels(hm.get_yticklabels(), rotation=0)
    ax.tick_params(left=False, bottom=False)
    sns.despine(top=False, right=False, left=False, bottom=False)
    plt.show()


if __name__ == "__main__":
    # this script compares two matrixes and plots the comparison matrix
    #
    # INPUT:
    # sys.argv[1] = yaml file 1 of the bird, example: example_yaml.yaml
    # sys.argv[2] = yaml file 1 of the bird, example: example_yaml.yaml
    #
    # OUTPUT:
    # figures

    main(sys.argv[1], sys.argv[2])
