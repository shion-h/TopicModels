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
from msmodule import BaseLDATopicNumberSelector, lda_directory_path
from bowmodule import make_bag_of_words_num, make_frequency_matrix
lda_directory_path = os.path.split(os.path.abspath(os.path.dirname(__file__)))[0]
timestamp = dt.now().strftime('%m%d_%a_%H_%M')


class LDAWBICCalculator(BaseLDATopicNumberSelector):
    def calculate_WBIC(self, topic_num):
        # command = lda_directory_path + "/LDAWBIC " + \
        #           self.BOW_filename + ' ' + lda_directory_path + \
        #           '/output/' + timestamp+self.ID+'/K_' + str(topic_num) + \
        #           '/alpha.csv ' + lda_directory_path + '/output/' + \
        #           timestamp+self.ID+'/K_' + str(topic_num) + '/beta.csv'
        command = lda_directory_path + "/LDAWBIC " + \
                  self.BOW_filename + ' ' + \
                  self.output_dir + \
                  '/K_' + str(topic_num) + '/alpha.csv ' + \
                  self.output_dir + \
                  '/K_' + str(topic_num) + '/beta.csv'
        subprocess.call(command.split(' '))

    def run(self):
        self.make_output_dir()
        pool = Pool(processes=5)
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
    parser.add_argument("-d", "--conv_det", help="convergence determination", type=float, default=0.001)
    parser.add_argument("-o", "--output_dir", help="output directory", type=str, default='')
    args = parser.parse_args()
    calculator = LDAWBICCalculator(args.BOW_filename, args.k_min, args.k_max, args.kstep, args.conv_det, args.output_dir)
    calculator.run()
