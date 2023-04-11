import re
import glob
import sys
import yaml
import helper_functions as hf

from IPython import embed


def get_data(cfg):

    l = glob.glob('D:/Birds/wh09pk39/otargeted_gtargeten/*/*.rec')
    catch_list = []
    templ_list = []
    for idx in range(len(l)):
        with open(l[idx], 'r') as file:
            lines = file.readlines()
            lines = [item.rstrip() for item in lines]
            # catch_list.append(lines[9][-1])
            templ_list.append(lines[-1])

    print(lines)
    embed()
    quit()



    file_list = glob.glob(cfg['paths']['folder_path'])

    seqs = hf.get_labels(file_list, cfg['labels']['intro_notes'])
    cfg['data']['bouts'], cfg['data']['noise'] = hf.get_bouts(seqs, cfg['labels']['bout_chunk'])

    for i in range(len(cfg['labels']['double_syl'])):
        if i == 0:
            cfg['data']['bouts_rep'] = re.sub(cfg['labels']['double_syl'][i],
                                              cfg['labels']['double_syl_rep'][i],
                                              cfg['data']['bouts'])
        else:
            cfg['data']['bouts_rep'] = re.sub(cfg['labels']['double_syl'][i],
                                              cfg['labels']['double_syl_rep'][i],
                                              cfg['data']['bouts_rep'])

    cfg['data']['chunk_bouts'] = hf.replace_chunks(cfg['data']['bouts_rep'], cfg['labels']['chunks'])

    return cfg


def main(yaml_file):
    with open(yaml_file) as f:
        cfg = yaml.load(f, Loader=yaml.FullLoader)
        f.close()
        cfg = get_data(cfg)

    with open(yaml_file, 'w') as f:
        yaml.dump(cfg, f)
        # print(yaml.dump(cfg))
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
