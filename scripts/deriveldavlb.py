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

import argparse
import numpy as np
from multiprocessing import Pool
from multiprocessing import cpu_count
from msmodule import BaseLDATopicNumberSelector


class LDAVLBDeriver(BaseLDATopicNumberSelector):
    def derive_VLB_each_topic_num(self):
        for topic_num in range(self.k_min, self.k_max, self.kstep):
            try:
                VLB = np.loadtxt(self.output_dir+'/K_'+str(topic_num)+'/variationalLowerBound.csv')
                print('{0},{1}'.format(topic_num, VLB))
            except FileNotFoundError:
                print('{0},{1}'.format(topic_num, 'Error'))

    def run(self):
        pool = Pool(processes=5)
        _ = pool.map(self.execute_estimation,
                 range(self.k_min, self.k_max, self.kstep))
        self.derive_VLB_each_topic_num()


if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument("BOW_filename", help="Bag Of Words filename")
    parser.add_argument("k_min", help="minimum number of topics for calculation of maximized variational lower bound", type=int)
    parser.add_argument("k_max", help="max number of topics for calculation of maximized variational lower bound", type=int)
    parser.add_argument("-s", "--kstep", help="step number of topics for calculation of maximized variational lower bound", type=int, default=1)
    parser.add_argument("-d", "--conv_det", help="convergence determination", type=float, default=0.001)
    parser.add_argument("-o", "--output_dir", help="output directory", type=str, default='')
    args = parser.parse_args()
    deriver = LDAVLBDeriver(args.BOW_filename, args.k_min, args.k_max, args.kstep, args.conv_det, args.output_dir)
    deriver.run()
