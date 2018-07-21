#
# calcldappl.py
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
from multiprocessing import Pool
from datetime import datetime as dt
from msmodule import BaseLDATopicNumberSelector, lda_directory_path
from bowmodule import make_bag_of_words_num, make_frequency_matrix


def calculate_lda_log_likelihood(theta, phi, BOW):
    log_likelihood = 0.0
    word_prob_df = theta.dot(phi)
    for d, document in enumerate(BOW):
        document_log_likelihood = 0.0
        for i, word in enumerate(document):
            document_log_likelihood += np.log(word_prob_df.iloc[d, word])
        log_likelihood += document_log_likelihood
    return log_likelihood


def calculate_lda_perplexity(theta, phi, test_BOW, test_words_num=0):
    if test_words_num == 0:
        for document in self.test_BOW:
            test_words_num += len(document)
    log_likelihood = calculate_lda_log_likelihood(theta, phi, test_BOW)
    perplexity = np.exp(-log_likelihood/test_words_num)
    return perplexity


class LDATestPerplexityCalculator(BaseLDATopicNumberSelector):
    def __init__(self, BOW_filename, k_min, k_max, kstep, rate_of_train_data, conv_det, output_dir):
        if (float(rate_of_train_data)>1) or (float(rate_of_train_data)<=0):
            exit()
        super().__init__(BOW_filename, k_min, k_max, kstep, conv_det, output_dir)
        self.BOW, self.word_list = make_bag_of_words_num(BOW_filename)
        self.rate_of_train_data = float(rate_of_train_data)
        self.make_data_dir()

    def make_data_dir(self):
        self.data_dir = self.output_dir

    def split_train_test(self):
        train_BOW = []
        test_BOW = []
        for d, document in enumerate(self.BOW):
            all_index = list(range(len(self.BOW[d])))
            all_index = list(range(len(self.BOW[d])))
            train_index = random.sample(all_index, int(len(self.BOW[d])*self.rate_of_train_data))
            test_index = list(set(all_index)-set(train_index))
            train_words = [self.BOW[d][i] for i in train_index]
            test_words = [self.BOW[d][i] for i in test_index]
            train_BOW.append(train_words)
            test_BOW.append(test_words)
        make_frequency_matrix(train_BOW, self.word_list,
                              self.data_dir + '/trainBOW')
        make_frequency_matrix(test_BOW, self.word_list,
                              self.data_dir + '/testBOW')
        self.test_BOW = test_BOW

    # override
    def create_command_string(self, topic_num):
        output_dir = self.output_dir + '/K_'+str(topic_num)
        try:
            os.mkdir(output_dir)
        except FileExistsError:
            pass
        command = lda_directory_path + '/LDA ' + self.data_dir + '/trainBOW'
        command += (' -k '+str(topic_num))
        command += (' -d '+str(self.conv_det))
        command += ' -o ' + output_dir + '/'
        return command

    def calculate_perplexity_each_topic_num(self):
        test_words_num = 0
        for document in self.test_BOW:
            test_words_num += len(document)
        for topic_num in range(self.k_min, self.k_max, self.kstep):
            theta, phi = parse_result(
                self.output_dir + '/K_'+str(topic_num) + '/theta.csv',
                self.output_dir + '/K_'+str(topic_num) + '/phi.csv',
                self.output_dir + '/K_'+str(topic_num) + '/wordList.csv'
            )
            perplexity = calculate_lda_perplexity(theta, phi, self.test_BOW, test_words_num)
            print('{0},{1}'.format(topic_num, perplexity))

    def run(self):
        self.split_train_test()
        pool = Pool(processes=5)
        _ = pool.map(self.execute_estimation,
                 range(self.k_min, self.k_max, self.kstep))
        self.calculate_perplexity_each_topic_num()


if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument("BOW_filename", help="Bag Of Words filename")
    parser.add_argument("k_min", help="minimum number of topics for calculation of perplexity", type=int)
    parser.add_argument("k_max", help="max number of topics for calculation of perplexity", type=int)
    parser.add_argument("-s", "--kstep", help="step number of topics for calculation of perplexity", type=int, default=1)
    parser.add_argument("-r", "--rate", help="rate of train data", default=0.7)
    parser.add_argument("-d", "--conv_det", help="convergence determination", type=float, default=0.001)
    parser.add_argument("-o", "--output_dir", help="output directory", type=str, default='')
    args = parser.parse_args()
    calculator = LDATestPerplexityCalculator(args.BOW_filename, args.k_min, args.k_max, args.kstep, args.rate, args.conv_det, args.output_dir)
    calculator.run()
