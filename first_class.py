import os
import re
import glob
import numpy as np
import plot_functions as pf
import helper_functions as hf
import matplotlib.pyplot as plt

from IPython import embed




class SongAnalysisParameters:
    def __init__(self):
        # constants ---------------------------------------------------------------------------------------------------
        self.threshold = 5

        # paths -------------------------------------------------------------------------------------------------------
        self.datapath = 'D:/Birds/bu01ye01'
        self.savepath = 'C:/Users/User/Documents/Doktorarbeit/data/bu01ye01/py_data/'
        self.folder = ['intan_wav']

        # intronotes, labels and chunks -------------------------------------------------------------------------------
        self.bout_chunk = 'he'  # chunk has to be in a bout to be analysed as a bout
        self.intronotes = ['S', 'E', 'k', 'c']  # all intro nodes are replaced by 'i',
                                                # 'S' and 'E' have to be in the list so that the code works
        self.labelsorder = [['E', 'S', 'I', 'P', 'p', 'j', 'B', 'b', 's', 'H', 'h', 'e',
                             'G', 'g', 'd', 'F', 'f', 'l', 'c']]  # unique labels; list of labels for each folder
        self.chunks = ['pj', 'bs', 'he', 'gd', 'he', 'fl', 'i+', 'c+', 'ES']

        # double syllables, replaces_double_syllables -----------------------------------------------------------------
        self.double_syl = ['bj']
        self.rep_double_syl = ['bs']

    def getdata_func(self):
        os.chdir(self.datapath)

        for folderidx in range(len(self.folder)):
            os.chdir(self.folder[folderidx])
            print(os.getcwd())
            days = [d for d in os.listdir(os.getcwd()) if os.path.isdir(d)]

            all_bouts = []

            for dayidx in range(1):
                os.chdir(days[dayidx])
                print(os.getcwd())
                seqs = hf.get_labels(self.intronotes)
                bouts, _ = hf.get_bouts(seqs, self.bout_chunk)
                all_bouts.append(bouts)
                os.chdir('..')


class PlotParameters():
    def __init__(self):
        # labels ---------------------------------------------------------------------------------------------------
        self.nodelabels = ['Start', 'i+', 'pj', 'p', 'bj', 'he', 'gd', 'fl']