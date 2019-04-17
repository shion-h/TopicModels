#
# derivesldavlb.py
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
import argparse
import numpy as np
from multiprocessing import Pool
from deriveldavlb import LDAVLBDeriver
from msmodule import lda_directory_path


class sLDAVLBDeriver(LDAVLBDeriver):
    def __init__(self, BOW_filename, label_filename, k_min, k_max, kstep, conv_det, output_dir):
        super().__init__(BOW_filename, k_min, k_max, kstep, conv_det, output_dir)
        self.label_filename = label_filename

    # override
    def create_command_string(self, topic_num):
        output_dir = self.output_dir + '/K_'+str(topic_num)
        try:
            os.mkdir(output_dir)
        except FileExistsError:
            pass
        command = lda_directory_path + '/sLDA ' + self.BOW_filename
        command += ' ' + self.label_filename
        command += (' -k '+str(topic_num))
        command += (' -d '+str(self.conv_det))
        command += ' -o ' + output_dir + '/'
        return command

    def run(self):
        pool = Pool(processes=5)
        _ = pool.map(self.execute_estimation,
                 range(self.k_min, self.k_max, self.kstep))
        self.derive_VLB_each_topic_num()


if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument("BOW_filename", help="Bag Of Words filename")
    parser.add_argument("label_filename", help="label filename")
    parser.add_argument("k_min", help="minimum number of topics for calculation of maximized variational lower bound", type=int)
    parser.add_argument("k_max", help="max number of topics for calculation of maximized variational lower bound", type=int)
    parser.add_argument("-s", "--kstep", help="step number of topics for calculation of maximized variational lower bound", type=int, default=1)
    parser.add_argument("-d", "--conv_det", help="convergence determination", type=float, default=0.001)
    parser.add_argument("-o", "--output_dir", help="output directory", type=str, default='')
    args = parser.parse_args()
    deriver = sLDAVLBDeriver(args.BOW_filename, args.label_filename, args.k_min, args.k_max, args.kstep, args.conv_det, args.output_dir)
    deriver.run()
