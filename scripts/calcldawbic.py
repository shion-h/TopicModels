#
# calcldawbic.py
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
from multiprocessing import Pool
from iomodule import create_command_string
from bowmodule import make_bag_of_words_num, make_frequency_matrix
from calcldappl import calculate_log_likelihood
lda_directory_path = os.path.split(os.path.abspath(os.path.dirname(__file__)))[0]
timestamp = dt.now().strftime('%m%d_%a_%H_%M')


class WBICCalculator():
    def __init__(self, BOW_filename, k_min, k_max, kstep, iteration):
        self.BOW_filename = BOW_filename
        self.k_min = int(k_min)
        self.k_max = int(k_max)
        self.kstep = int(kstep)
        self.iteration = int(iteration)
        self.ID = ''

    def make_output_dir(self):
        count = 0
        while True:
            try:
                os.mkdir(lda_directory_path + '/output/' + timestamp + '_' + str(count))
                self.ID = '_' + str(count)
                break
            except FileExistsError:
                pass
            count += 1

    def execute_estimation(self,topic_num):
        output_dir = lda_directory_path + '/output/' + timestamp + self.ID + '/K_'+str(topic_num)
        try:
            os.mkdir(output_dir)
        except FileExistsError:
            pass
        command = create_command_string(topic_num, lda_directory_path, 
                                        output_dir, self.BOW_filename, 
                                        iteration=self.iteration)
        subprocess.call(command.split(' '),
                        stdout=subprocess.DEVNULL,
                        stderr=subprocess.STDOUT)

    def calculate_WBIC(self, topic_num):
        command = lda_directory_path + "/LDAWBIC " + \
                  self.BOW_filename + ' ' + lda_directory_path + \
                  '/output/' + timestamp+self.ID+'/K_' + str(topic_num) + \
                  '/alpha.csv ' + lda_directory_path + '/output/' + \
                  timestamp+self.ID+'/K_' + str(topic_num) + '/beta.csv'
        # print(command)
        subprocess.call(command.split(' '))

    def run(self):
        self.make_output_dir()
        pool = Pool(processes=4)
        _ = pool.map(self.execute_estimation,
                 range(self.k_min, self.k_max, self.kstep))
        _ = pool.map(self.calculate_WBIC,
                 range(self.k_min, self.k_max, self.kstep))


if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument("BOW_filename", help="Bag Of Words filename")
    parser.add_argument("k_min", help="minimum number of topics for calculation of wbic", type=int)
    parser.add_argument("k_max", help="max number of topics for calculation of wbic", type=int)
    parser.add_argument("-s", "--kstep", help="step number of topics for calculation of wbic", type=int, default=1)
    parser.add_argument("-i", "--iteration", help="number of iterations(>10)", default=100)
    args = parser.parse_args()
    calculator = WBICCalculator(args.BOW_filename, args.k_min, args.k_max, args.kstep, args.iteration)
    calculator.run()
