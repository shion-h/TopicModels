#
# calcldavlb.py
#
# Copyright (c) 2017 Shion Hosoda
#
# This software is released under the MIT License.
# http://opensource.org/licenses/mit-license.php
#
#
# Includes pandas
#-----------------------------------------------------------------------------
# Copyright (c) 2012, PyData Development Team
# All rights reserved.
#
# Distributed under the terms of the BSD Simplified License.
#
# The full license is in the LICENSE file, distributed with this software.
#-----------------------------------------------------------------------------

import os
import sys
import random
import argparse
import subprocess
import numpy as np
import pandas as pd
from scipy.stats import dirichlet
from datetime import datetime as dt
from iomodule import create_command_string
from bowmodule import make_bag_of_words_num, make_frequency_matrix
from calcldappl import calculate_log_likelihood
lda_directory_path = os.path.split(os.path.abspath(os.path.dirname(__file__)))[0]
timestamp = dt.now().strftime('%m%d_%a_%H_%M_%S')


class VLBCalculator():
    def __init__(self, BOW_filename, k_min, k_max, kstep):
        self.BOW_filename = BOW_filename
        self.k_min = int(k_min)
        self.k_max = int(k_max)
        self.kstep = int(kstep)
        self.iteration = 100 # dummy

    def execute_estimation(self):
        try:
            os.mkdir(lda_directory_path + '/output/' + timestamp)
        except FileExistsError:
            exit()
        for topic_num in range(self.k_min, self.k_max, self.kstep):
            sys.stdout.write('estimating(' +str(topic_num)+ 'topics-LDA)...\r')
            output_dir = lda_directory_path + '/output/' + timestamp + '/K_'+str(topic_num)
            try:
                os.mkdir(output_dir)
            except FileExistsError:
                pass
            command = create_command_string(topic_num, self.iteration,
                                            lda_directory_path, output_dir,
                                            self.BOW_filename, algorythm_flag=2)
            subprocess.call(command.split(' '),
                            stdout=subprocess.DEVNULL, stderr=subprocess.STDOUT)

    def calculate_VLB_each_topic_num(self):
        BOW, _ = make_bag_of_words_num(self.BOW_filename)
        words_num = 0
        for document in BOW:
            words_num += len(document)
        for topic_num in range(self.k_min, self.k_max, self.kstep):
            VLB = np.loadtxt(lda_directory_path+'/output/'+
                    timestamp+'/K_'+str(topic_num)+'/variationalLowerBound.csv')
            print('{0},{1}'.format(topic_num, VLB))

    def run(self):
        self.execute_estimation()
        print('\nestimate done')
        self.calculate_VLB_each_topic_num()


if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument("BOW_filename", help="Bag Of Words filename")
    parser.add_argument("k_min", help="minimum number of topics for calculation of perplexity", type=int)
    parser.add_argument("k_max", help="max number of topics for calculation of perplexity", type=int)
    parser.add_argument("-s", "--kstep", help="step number of topics for calculation of perplexity", type=int, default=1)
    args = parser.parse_args()
    calculator = VLBCalculator(args.BOW_filename, args.k_min, args.k_max, args.kstep)
    calculator.run()

