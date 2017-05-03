#
# calculate_perplexity.py
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
from datetime import datetime as dt
from iomodule import create_command_string, parse_result
from bowmodule import make_bag_of_words_num, make_frequency_matrix

lda_directory_path = os.path.split(os.path.abspath(os.path.dirname(__file__)))[0]
timestamp = dt.now().strftime('%m%d_%a_%H_%M_%S')


def calculate_perplexity(theta, phi, test_BOW, test_words_num=0):
    if test_words_num == 0:
        for document in self.test_BOW:
            test_words_num += len(document)
    log_likelihood = 0.0
    word_prob_df = theta.dot(phi)
    for d, document in enumerate(test_BOW):
        document_log_likelihood = 0.0
        for i, word in enumerate(document):
            document_log_likelihood += np.log(word_prob_df.iloc[d, word])
        log_likelihood += document_log_likelihood
    perplexity = np.exp(-log_likelihood/test_words_num)
    return perplexity


class PerplexityCalculator():
    def __init__(self, BOW_filename, k_min, k_max, kstep, rate_of_train_data, iteration):
        if (float(rate_of_train_data)>1) & (float(rate_of_train_data)<=0):
            exit()
        self.BOW, self.word_list = make_bag_of_words_num(BOW_filename)
        self.k_min = int(k_min)
        self.k_max = int(k_max)
        self.kstep = int(kstep)
        self.rate_of_train_data = float(rate_of_train_data)
        self.iteration = int(iteration)

    # TODO:訓練perplexityも計算
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
        try:
            os.mkdir(lda_directory_path + '/data/' + timestamp)
        except FileExistsError:
            pass
        make_frequency_matrix(train_BOW, self.word_list,
                              lda_directory_path + '/data/' + timestamp + '/trainBOW')
        make_frequency_matrix(test_BOW, self.word_list,
                              lda_directory_path+'/data/' + timestamp + '/testBOW')
        self.test_BOW = test_BOW

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
                                            lda_directory_path, output_dir)
            subprocess.call(command.split(' '),
                            stdout=subprocess.DEVNULL, stderr=subprocess.STDOUT)

    def calculate_perplexity_each_topic_num(self):
        test_words_num = 0
        for document in self.test_BOW:
            test_words_num += len(document)
        for topic_num in range(self.k_min, self.k_max, self.kstep):
            theta, phi = parse_result(
                lda_directory_path+'/output/'+timestamp+'/K_'+str(topic_num)+'/theta.csv',
                lda_directory_path+'/output/'+timestamp+'/K_'+str(topic_num)+'/phi.csv',
                lda_directory_path+'/output/'+timestamp+'/K_'+str(topic_num)+'/wordList.csv'
            )
            perplexity = calculate_perplexity(theta, phi, self.test_BOW, test_words_num)
            print('{0},{1}'.format(topic_num, perplexity))


    def run(self):
        self.split_train_test()
        print('split done')
        self.execute_estimation()
        print('\nestimate done')
        self.calculate_perplexity_each_topic_num()


if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument("BOW_filename", help="Bag Of Words filename")
    parser.add_argument("k_min", help="minimum number of topics for calculation of perplexity", type=int)
    parser.add_argument("k_max", help="max number of topics for calculation of perplexity", type=int)
    parser.add_argument("-s", "--kstep", help="step number of topics for calculation of perplexity", type=int, default=1)
    parser.add_argument("-r", "--rate", help="rate of train data", default=0.7)
    parser.add_argument("-i", "--iteration", help="number of iterations(>10)", default=100)
    args = parser.parse_args()
    calculator = PerplexityCalculator(args.BOW_filename, args.k_min, args.k_max, args.kstep, args.rate, args.iteration)
    calculator.run()
