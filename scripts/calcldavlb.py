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
from multiprocessing import Pool
from iomodule import create_command_string
from bowmodule import make_bag_of_words_num, make_frequency_matrix
from calcldappl import calculate_log_likelihood
lda_directory_path = os.path.split(os.path.abspath(os.path.dirname(__file__)))[0]
timestamp = dt.now().strftime('%m%d_%a_%H_%M')


class VLBCalculator():
    def __init__(self, BOW_filename, k_min, k_max, kstep, conv_det):
        self.BOW_filename = BOW_filename
        self.k_min = int(k_min)
        self.k_max = int(k_max)
        self.kstep = int(kstep)
        self.conv_det = float(conv_det)
        self.iteration = 100 # dummy
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

    def execute_estimation(self, topic_num):
        output_dir = lda_directory_path + '/output/' + timestamp + self.ID + '/K_'+str(topic_num)
        try:
            os.mkdir(output_dir)
        except FileExistsError:
            pass
        command = create_command_string(topic_num,lda_directory_path,
                                        output_dir,self.BOW_filename,
                                        algorythm_flag=2, conv_det=self.conv_det)
        subprocess.call(command.split(' '),
                        stdout=subprocess.DEVNULL, stderr=subprocess.STDOUT)
        # p = subprocess.Popen(command.split(' '),
                        # stdout=subprocess.PIPE, stderr=subprocess.PIPE)
        # print(p.communicate())

    def derive_VLB_each_topic_num(self):
        BOW, _ = make_bag_of_words_num(self.BOW_filename)
        words_num = 0
        for document in BOW:
            words_num += len(document)
        vlb_list = []
        for topic_num in range(self.k_min, self.k_max, self.kstep):
            VLB = np.loadtxt(lda_directory_path+'/output/'+
                             timestamp+self.ID+'/K_'+str(topic_num)+'/variationalLowerBound.csv')
            print('{0},{1}'.format(topic_num, VLB))

    def run(self):
        self.make_output_dir()
        pool = Pool(processes=4)
        _ = pool.map(self.execute_estimation,
                 range(self.k_min, self.k_max, self.kstep))
        self.derive_VLB_each_topic_num()


if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument("BOW_filename", help="Bag Of Words filename")
    parser.add_argument("k_min", help="minimum number of topics for calculation of maximized variational lower bound", type=int)
    parser.add_argument("k_max", help="max number of topics for calculation of maximized variational lower bound", type=int)
    parser.add_argument("-s", "--kstep", help="step number of topics for calculation of maximized variational lower bound", type=int, default=1)
    parser.add_argument("-d", "--conv_det", help="convergence determination for vb algorythm", type=float, default=0.001)
    args = parser.parse_args()
    calculator = VLBCalculator(args.BOW_filename, args.k_min, args.k_max, args.kstep, args.conv_det)
    calculator.run()
