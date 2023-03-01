import re
import sys
import glob
import yaml
import numpy as np
import plot_functions as pf
import helper_functions as hf
import matplotlib.pyplot as plt

from IPython import embed


def get_data(cfg):
    file_list = glob.glob(cfg['paths']['folder_path'] + '/*/*.not.mat')

    seqs = hf.get_labels(file_list, cfg['labels']['intro_notes'])
    cfg['data']['bouts'], cfg['data']['noise'] = hf.get_bouts(seqs, cfg['labels']['bout_chunk'])

    for i in range(len(cfg['labels']['double_syl'])):
        cfg['data']['bouts_rep'] = re.sub(cfg['labels']['double_syl'][i], cfg['labels']['double_syl_rep'][i], cfg['data']['bouts'])

    cfg['data']['chunk_bouts'] = hf.replace_chunks(cfg['data']['bouts_rep'], cfg['labels']['chunks'])

    return cfg


def make_first_plots(cfg):
    bouts = cfg['data']['bouts_rep']
    tm, _ = hf.get_transition_matrix(bouts, cfg['labels']['unique_labels'])

    # ---------------------------------------------------------------------------------------------------------------
    k = np.where(sum(tm) / sum(sum(tm)) * 100 <= 0)
    tmd = np.delete(tm, k, axis=1)
    tmd = np.delete(tmd, k, axis=0)

    tmpd = (tmd.T / np.sum(tmd, axis=1)).T
    tmpd = hf.get_node_matrix(tmpd, cfg['constants']['edge_threshold'])
    'Plot Transition Matrix and Transition Diagram'
    node_size = np.round(np.sum(tmd, axis=1) / np.min(np.sum(tmd, axis=1)), 2) * 100
    pf.plot_transition_diagram(tmpd, np.delete(cfg['labels']['unique_labels'], k),
                               node_size,
                               cfg['paths']['save_path']+cfg['title_figures']+'_graph_simple.jpg',
                               cfg['title_figures'])
    plt.show()


def make_final_plots(cfg):
    bouts = cfg['data']['chunk_bouts']
    tm, tmp = hf.get_transition_matrix(bouts, cfg['labels']['unique_labels'])

    # ---------------------------------------------------------------------------------------------------------------
    k = np.where(sum(tm) / sum(sum(tm)) * 100 <= cfg['constants']['node_threshold'])
    tmd = np.delete(tm, k, axis=1)
    tmd = np.delete(tmd, k, axis=0)

    tmpd = (tmd.T / np.sum(tmd, axis=1)).T
    tmpd = hf.get_node_matrix(tmpd, cfg['constants']['edge_threshold'])
    'Plot Transition Matrix and Transition Diagram'
    pf.plot_transition_matrix(tmpd,
                              cfg['labels']['node_labels'],
                              cfg['paths']['save_path'] + cfg['title_figures'] + '_matrix.jpg',
                              cfg['title_figures'])
    pf.plot_transition_diagram(tmpd,
                               cfg['labels']['node_labels'],
                               np.round(np.sum(tmd, axis=1) / np.min(np.sum(tmd, axis=1)), 2) * 500,
                               cfg['paths']['save_path']+cfg['title_figures']+'_graph.jpg',
                               cfg['title_figures'])
    plt.show()


def main(yaml_file):

    with open(yaml_file) as f:
        cfg = yaml.load(f, Loader=yaml.FullLoader)
        f.close()

    if not cfg['data']['bouts']:
        cfg = get_data(cfg)
        make_first_plots(cfg)

    else:
        make_final_plots(cfg)

    with open(yaml_file, 'w') as f:
        yaml.dump(cfg, f)
        f.close()
    # embed()
    # quit()



if __name__ == '__main__':
    # this script plots transition matrix and diagrams
    #
    # INPUT:
    # sys.argv[1] = yaml file of the bird, example: example_yaml.yaml
    #
    # OUTPUT:
    # figures

    main(sys.argv[1])
