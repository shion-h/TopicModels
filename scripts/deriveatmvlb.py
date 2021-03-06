#
# deriveldavlb.py
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
from atmiomodule import create_command_string
from bowmodule import make_bag_of_words_num, make_frequency_matrix
from calcldappl import calculate_log_likelihood
lda_directory_path = os.path.split(os.path.abspath(os.path.dirname(__file__)))[0]
timestamp = dt.now().strftime('%m%d_%a_%H_%M')


class VLBDeriver():
    def __init__(self, BOW_filename, author_filename, k_min, k_max, kstep, conv_det, output_dir):
        self.BOW_filename = BOW_filename
        self.author_filename = author_filename
        self.k_min = int(k_min)
        self.k_max = int(k_max)
        self.kstep = int(kstep)
        self.conv_det = float(conv_det)
        self.output_dir = output_dir
        if len(output_dir) != 0:
            if output_dir[-1] == '/':
                self.output_dir = output_dir[:-1]

    def make_output_dir(self):
        count = 0
        while True:
            try:
                os.mkdir(lda_directory_path + '/output/' + timestamp + '_' + str(count))
                self.output_dir = lda_directory_path + '/output/' + timestamp +'_' + str(count)
                break
            except FileExistsError:
                pass
            count += 1

    def execute_estimation(self, topic_num):
        output_dir = self.output_dir + '/K_'+str(topic_num)
        try:
            os.mkdir(output_dir)
        except FileExistsError:
            pass
        command = create_command_string(topic_num,lda_directory_path,
                                        output_dir,self.BOW_filename,
                                        self.author_filename, conv_det=self.conv_det)
        subprocess.call(command.split(' '),
                        stdout=subprocess.DEVNULL, stderr=subprocess.STDOUT)

    def derive_VLB_each_topic_num(self):
        for topic_num in range(self.k_min, self.k_max, self.kstep):
            try:
                VLB = np.loadtxt(self.output_dir+'/K_'+str(topic_num)+'/variationalLowerBound.csv')
                print('{0},{1}'.format(topic_num, VLB))
            except FileNotFoundError:
                print('{0},{1}'.format(topic_num, 'Error'))

    def run(self):
        if not os.path.isdir(self.output_dir):
            self.make_output_dir()
        pool = Pool(processes=5)
        _ = pool.map(self.execute_estimation,
                 range(self.k_min, self.k_max, self.kstep))
        self.derive_VLB_each_topic_num()


if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument("BOW_filename", help="Bag Of Words filename")
    parser.add_argument("author_filename", help="Author filename")
    parser.add_argument("k_min", help="minimum number of topics for calculation of maximized variational lower bound", type=int)
    parser.add_argument("k_max", help="max number of topics for calculation of maximized variational lower bound", type=int)
    parser.add_argument("-s", "--kstep", help="step number of topics for calculation of maximized variational lower bound", type=int, default=1)
    parser.add_argument("-d", "--conv_det", help="convergence determination", type=float, default=0.001)
    parser.add_argument("-o", "--output_dir", help="output directory", type=str, default='')
    args = parser.parse_args()
    deriver = VLBDeriver(args.BOW_filename, args.author_filename, args.k_min, args.k_max, args.kstep, args.conv_det, args.output_dir)
    deriver.run()
