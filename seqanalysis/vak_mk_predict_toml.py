""" This skript makes predict toml files for an entire folder with dates in it

INPUT:
    * sys.argv[1] = toml file of the bird, example: example_toml.yaml
    * sys.argv[2] = path to file folder
    * sys.argv[3] = birdname

OUTPUT:
    * predict_toml files for every day in the folder
"""

import os
import sys
import toml
from IPython import embed

def main(toml_file, folder_path, birdname):

    old_file = toml.load(toml_file)
    old_data_dir = old_file['PREP']['data_dir']

    save_file_path = toml_file.rsplit('\\', 1)[0]
    date_list = os.listdir(folder_path)

    for i in range(len(date_list)):

        file_path = save_file_path+'/'+birdname+'_'+date_list[i]+'_predict.toml'
        with open(file_path, 'w') as file:
            new_file = old_file
            new_file['PREP']['data_dir'] = old_data_dir.rsplit('/', 1)[0] + '/' + date_list[i]
            new_file['PREDICT']['annot_csv_filename'] = birdname + '_' + date_list[i] + '_annot.csv'
            toml.dump(new_file, file)


if __name__ == '__main__':
    main(sys.argv[1], sys.argv[2], sys.argv[3])
