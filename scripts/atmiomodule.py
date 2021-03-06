#
# atmiomodule.py
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


import pandas as pd


def create_command_string(topic_num, lda_directory_path, output_dir, BOW_filename, author_filename,  conv_det=0.001):
    command = lda_directory_path + '/ATM ' + BOW_filename + ' ' + author_filename
    command += (' -k '+str(topic_num))
    command += (' -d '+str(conv_det))
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
