import glob
import pathlib
import re
import sys

import matplotlib.pyplot as plt
import numpy as np
import yaml
from IPython import embed

import seqanalysis.util.helper_functions as hf
import seqanalysis.util.plot_functions as pf
from seqanalysis.util.logging import config_logging

log = config_logging()


def get_catch_data(cfg):
    # assembles the files in the path
    file_list = []
    files_path = pathlib.Path(cfg["paths"]["catch_path"])
    if not files_path.exists():
        log.error(f"Path {files_path} does not exist")
        FileNotFoundError(f"Path {files_path} does not exist")
    files = list(files_path.glob("*/*.not.mat"))
    for i in range(len(files)):
        with open(files[i] + cfg["paths"]["catch_file"], "r") as file:
            line_list = file.readlines()
            file_list.extend(
                [list[i] + item.rstrip() + ".not.mat" for item in line_list]
            )

    seqs = hf.get_labels(file_list, cfg["labels"]["intro_notes"])
    cfg["data"]["bouts"], cfg["data"]["noise"] = hf.get_bouts(
        seqs, cfg["labels"]["bout_chunk"]
    )
    if cfg["labels"]["double_syl"] != [None]:
        for i in range(len(cfg["labels"]["double_syl"])):
            if i == 0:
                cfg["data"]["bouts_rep"] = re.sub(
                    cfg["labels"]["double_syl"][i],
                    cfg["labels"]["double_syl_rep"][i],
                    cfg["data"]["bouts"],
                )
            else:
                cfg["data"]["bouts_rep"] = re.sub(
                    cfg["labels"]["double_syl"][i],
                    cfg["labels"]["double_syl_rep"][i],
                    cfg["data"]["bouts_rep"],
                )
    else:
        cfg["data"]["bouts_rep"] = cfg["data"]["bouts"]

    cfg["data"]["chunk_bouts"] = hf.replace_chunks(
        cfg["data"]["bouts_rep"], cfg["labels"]["chunks"]
    )

    return cfg


def get_data(cfg):
    file_list = pathlib.Path((cfg["paths"]["folder_path"]))
    if not file_list.exists():
        log.error(f"Path {file_list} does not exist")
        FileNotFoundError(f"Path {file_list} does not exist")
    file_list = list(file_list.glob("*/*.not.mat"))
    if not file_list:
        log.error(f"No files found in {file_list}")
        FileNotFoundError(f"No files found in {file_list}")
    log.info(f"Files found: {len(file_list)}")

    seqs = hf.get_labels(file_list, cfg["labels"]["intro_notes"])
    cfg["data"]["bouts"], cfg["data"]["noise"] = hf.get_bouts(
        seqs, cfg["labels"]["bout_chunk"]
    )
    if cfg["labels"]["double_syl"]:
        log.info("Replacing double syllables")
        for i in range(len(cfg["labels"]["double_syl"])):
            if i == 0:
                cfg["data"]["bouts_rep"] = re.sub(
                    cfg["labels"]["double_syl"][i],
                    cfg["labels"]["double_syl_rep"][i],
                    cfg["data"]["bouts"],
                )
            else:
                cfg["data"]["bouts_rep"] = re.sub(
                    cfg["labels"]["double_syl"][i],
                    cfg["labels"]["double_syl_rep"][i],
                    cfg["data"]["bouts_rep"],
                )
    else:
        cfg["data"]["bouts_rep"] = cfg["data"]["bouts"]

    log.info("Replacing chunks")
    cfg["data"]["chunk_bouts"] = hf.replace_chunks(
        cfg["data"]["bouts_rep"], cfg["labels"]["chunks"]
    )

    return cfg


def make_first_plots(cfg):
    bouts = cfg["data"]["bouts_rep"]
    # unique_labels in bouts
    unique_labels = sorted(list(set(bouts)))
    log.info(f"Unique labels of Chunks from bouts_rep: {unique_labels}\n")
    tm, _ = hf.get_transition_matrix(bouts, unique_labels)

    # U2 is the size of the strings in zeros
    label_matrix = np.zeros((len(unique_labels), len(unique_labels)), "U2")
    for i, labely in enumerate(unique_labels):
        for j, labelx in enumerate(unique_labels):
            labelylabelx = str(labely + labelx)
            label_matrix[i, j] = labelylabelx

    # NOTE: Sort tm by the transitions with the highest probability
    tm_prob = np.around((tm.T / np.sum(tm, axis=0)).T, 2) * 100
    tm_sorted = np.zeros(tm.shape)
    label_matrix_sorted = np.zeros((len(unique_labels), len(unique_labels)), "U2")
    _multiple_index = [0]
    # NOTE: Add the first element befor entering the loop
    tm_sorted[0] = tm[0]
    label_matrix_sorted[0] = label_matrix[0]
    for i in range(1, tm.shape[0]):
        for sort in np.argsort(tm_sorted[i - 1])[::-1]:
            log.info(sort)
            if sort not in _multiple_index:
                log.info(_multiple_index)
                tm_sorted[i] = tm[sort]
                label_matrix_sorted[i] = label_matrix[sort]
                _multiple_index.append(sort)
                log.info(_multiple_index)
                break
            else:
                continue

    k = np.where(
        np.sum(tm, axis=0) / np.sum(tm) * 100 <= cfg["constants"]["node_threshold"]
    )
    tmd = np.delete(tm, k, axis=1)
    tmd = np.delete(tmd, k, axis=0)

    tmpd = (tmd.T / np.sum(tmd, axis=0)).T
    tmpd = hf.get_node_matrix(tmpd, cfg["constants"]["edge_threshold"])
    # "Plot Transition Matrix and Transition Diagram"
    node_size = np.round(np.sum(tmd, axis=1) / np.min(np.sum(tmd, axis=1)), 2) * 100
    # get them into the right order

    pf.plot_transition_diagram(
        tmpd,
        np.delete(unique_labels, k),
        node_size,
        cfg["constants"]["edge_width"],
        cfg["paths"]["save_path"] + cfg["title_figures"] + "_graph_simple.pdf",
        cfg["title_figures"],
    )

    pf.plot_transition_matrix(
        tmpd,
        np.delete(unique_labels, k),
        cfg["paths"]["save_path"] + cfg["title_figures"] + "_matrix_simple.pdf",
        cfg["title_figures"],
    )
    plt.show()


def make_final_plots(cfg):
    bouts = cfg["data"]["chunk_bouts"]
    tm, tmp = hf.get_transition_matrix(bouts, cfg["labels"]["unique_labels"])
    print(tmp)

    # Filter out nodes with low occurrence
    k = np.where(sum(tm) / sum(sum(tm)) * 100 <= cfg["constants"]["node_threshold"])
    tmd = np.delete(tm, k, axis=1)
    tmd = np.delete(tmd, k, axis=0)
    print(np.delete(cfg["labels"]["unique_labels"], k))

    # Normalize transition matrix and create node matrix
    tmpd = (tmd.T / np.sum(tmd, axis=1)).T
    tmpd = hf.get_node_matrix(tmpd, cfg["constants"]["edge_threshold"])
    node_size = (
        np.round(np.sum(tmd, axis=1) / np.min(np.sum(tmd, axis=1)), 2)
        * cfg["constants"]["node_size"]
    )
    "Plot Transition Matrix and Transition Diagram"
    pf.plot_transition_matrix(
        tmpd,
        cfg["labels"]["node_labels"],
        cfg["paths"]["save_path"] + cfg["title_figures"] + "_matrix.pdf",
        cfg["title_figures"],
    )
    pf.plot_transition_diagram(
        tmpd,
        cfg["labels"]["node_labels"],
        node_size,
        cfg["constants"]["edge_width"],
        cfg["paths"]["save_path"] + cfg["title_figures"] + "_graph.pdf",
        cfg["title_figures"],
    )
    plt.show()


def main(yaml_file, analyse_files):
    with open(yaml_file) as f:
        cfg = yaml.load(f, Loader=yaml.FullLoader)
        f.close()

    if not cfg["data"]["bouts"]:
        log.info("No bouts found in yaml file")
        if analyse_files == "all":
            cfg = get_data(cfg)
        elif analyse_files == "catch":
            cfg = get_catch_data(cfg)

        make_first_plots(cfg)
        log.info(
            f"Unique labels of Chunks: {sorted(list(set(cfg["data"]["chunk_bouts"])))}\n"
        )
        make_final_plots(cfg)

    else:
        log.info(
            f"Unique labels of Chunks: {sorted(list(set(cfg["data"]["chunk_bouts"])))}\n"
        )
        # print('{} Unique labels of Chunks: '.format(sorted(list(set(cfg['data']['chunk_bouts'])))))
        make_final_plots(cfg)

    with open(yaml_file, "w") as f:
        yaml.dump(cfg, f)
        # print(yaml.dump(cfg))
        f.close()


if __name__ == "__main__":
    # this script plots transition matrix and diagrams
    #
    # INPUT:
    # sys.argv[1] = yaml file of the bird, example: example_yaml.yaml
    # sys.argv[2] = analysis catch or all files: input: catch, all
    #
    # OUTPUT:
    # figures

    main(sys.argv[1], sys.argv[2])
