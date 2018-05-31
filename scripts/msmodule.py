#
# msmodule.py
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
import subprocess
import pandas as pd
from datetime import datetime as dt
lda_directory_path = os.path.split(os.path.abspath(os.path.dirname(__file__)))[0]
timestamp = dt.now().strftime('%m%d_%a_%H_%M')


class BaseLDATopicNumberSelector(object):
    def __init__(self, BOW_filename, k_min, k_max, kstep, conv_det, output_dir):
        self.BOW_filename = BOW_filename
        self.k_min = int(k_min)
        self.k_max = int(k_max)
        self.kstep = int(kstep)
        self.conv_det = float(conv_det)
        self.output_dir = output_dir
        if len(output_dir) != 0:
            if output_dir[-1] == '/':
                self.output_dir = output_dir[:-1]
        if not os.path.isdir(self.output_dir):
            self.make_output_dir()

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
        command = self.create_command_string(topic_num)
        subprocess.call(command.split(' '),
                        stdout=subprocess.DEVNULL, stderr=subprocess.STDOUT)

    def create_command_string(self, topic_num):
        output_dir = self.output_dir + '/K_'+str(topic_num)
        try:
            os.mkdir(output_dir)
        except FileExistsError:
            pass
        command = lda_directory_path + '/LDA ' + self.BOW_filename
        command += (' -k '+str(topic_num))
        command += (' -d '+str(self.conv_det))
        command += ' -o ' + output_dir + '/'
        return command


def parse_result(theta_estimated_filename, phi_estimated_filename, word_list_filename):
    theta = pd.read_csv(theta_estimated_filename, header=None)
    phi = pd.read_csv(phi_estimated_filename, header=None)
    phi.columns = phi.columns.astype('str')
    word_list = pd.read_csv(word_list_filename, header=None)
    word_list[0] = word_list[0].astype('str')
    phi.columns = word_list.iloc[:, 0]
    return theta, phi
