import re
import sys
import glob
import yaml
import numpy as np
import matplotlib.pyplot as plt

import seqanalysis.util.plot_functions as pf
import seqanalysis.util.helper_functions as hf


def get_catch_data(cfg):
    # assembles the files in the path
    file_list = []
    list = glob.glob(cfg["paths"]["catch_path"])
    for i in range(len(list)):
        with open(list[i] + cfg["paths"]["catch_file"], "r") as file:
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
    file_list = glob.glob(cfg["paths"]["folder_path"])

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
        cfg['data']['bouts_rep'] = cfg['data']['bouts']

    cfg["data"]["chunk_bouts"] = hf.replace_chunks(
        cfg["data"]["bouts_rep"], cfg["labels"]["chunks"]
    )

    return cfg


def make_first_plots(cfg):
    bouts = cfg["data"]["bouts_rep"]
    tm, _ = hf.get_transition_matrix(bouts, cfg["labels"]["unique_labels"])

    # ---------------------------------------------------------------------------------------------------------------
    k = np.where(sum(tm) / sum(sum(tm)) * 100 <= 0)
    tmd = np.delete(tm, k, axis=1)
    tmd = np.delete(tmd, k, axis=0)

    tmpd = (tmd.T / np.sum(tmd, axis=1)).T
    tmpd = hf.get_node_matrix(tmpd, cfg["constants"]["edge_threshold"])
    "Plot Transition Matrix and Transition Diagram"
    node_size = np.round(np.sum(tmd, axis=1) / np.min(np.sum(tmd, axis=1)), 2) * 100
    pf.plot_transition_diagram(
        tmpd,
        np.delete(cfg["labels"]["unique_labels"], k),
        node_size,
        cfg["constants"]["edge_width"],
        cfg["paths"]["save_path"] + cfg["title_figures"] + "_graph_simple.pdf",
        cfg["title_figures"],
    )

    pf.plot_transition_matrix(
        tmpd,
        np.delete(cfg["labels"]["unique_labels"], k),
        cfg["paths"]["save_path"] + cfg["title_figures"] + "_matrix_simple.pdf",
        cfg["title_figures"],
    )
    plt.show()


def make_final_plots(cfg):
    bouts = cfg["data"]["chunk_bouts"]
    tm, tmp = hf.get_transition_matrix(bouts, cfg["labels"]["unique_labels"])
    print(tmp)

    # Filter out nodes with low occurrence
    k = np.where(sum(tm) / sum(sum(tm)) * 100 <= cfg['constants']['node_threshold'])
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
=======
    tmpd = hf.get_node_matrix(tmpd, cfg['constants']['edge_threshold'])

    # Plot Transition Matrix and Transition Diagram
    node_size = np.round(np.sum(tmd, axis=1) / np.min(np.sum(tmd, axis=1)), 2) * cfg['constants']['node_size']

    pf.plot_transition_matrix(tmpd,
                              cfg['labels']['node_labels'],
                              cfg['paths']['save_path'] + cfg['title_figures'] + '_matrix.pdf',
                              cfg['title_figures'])
    pf.plot_transition_diagram(tmpd,
                               cfg['labels']['node_labels'],
                               node_size,
                               cfg['constants']['edge_width'],
                               cfg['paths']['save_path'] + cfg['title_figures'] + '_graph.pdf',
                               cfg['title_figures'])
>>>>>>> main
    plt.show()


def main(yaml_file, analyse_files):
    with open(yaml_file) as f:
        cfg = yaml.load(f, Loader=yaml.FullLoader)
        f.close()

    if not cfg["data"]["bouts"]:
        if analyse_files == "all":
            cfg = get_data(cfg)
        elif analyse_files == "catch":
            cfg = get_catch_data(cfg)

        make_first_plots(cfg)
        print(
            "\n Unique labels of Chunks: {}\n".format(
                sorted(list(set(cfg["data"]["chunk_bouts"])))
            )
        )

    else:
        # print('{} Unique labels of Chunks: '.format(sorted(list(set(cfg['data']['chunk_bouts'])))))
        make_final_plots(cfg)

    with open(yaml_file, "w") as f:
        yaml.dump(cfg, f)
        # print(yaml.dump(cfg))
        f.close()
    # embed()
    # quit()


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
