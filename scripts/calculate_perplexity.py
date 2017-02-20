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

lda_directory_path=os.path.split(os.path.abspath(os.path.dirname(__file__)))[0]

def read_BOW_file(BOW_filename):
    f=open(BOW_filename, 'r')
    full_string=f.read()
    f.close()
    BOW=[]
    for line in full_string.split('\n'):
        document=[]
        for word in line.split(','):
            if word!='':
                document.append(word)
        if document!=[]:
            BOW.append(document)
    return BOW


def write_BOW_file(BOW, BOW_filename):
    f=open(BOW_filename, 'w')
    for d, line in enumerate(BOW):
        for i, word in enumerate(line):
            f.write(word)
            if i!=(len(line)-1):
                f.write(',')
        f.write('\n')


def create_command_string(topic_num, iteration):
    command=lda_directory_path+'/bin/LDA '+lda_directory_path+'/data/trainBOW'
    command+=(' -k '+str(topic_num))
    command+=(' -s '+str(iteration))
    command+=(' -b '+str(int(iteration*0.8)))
    command+=(' -i '+str(5))
    try:
        os.mkdir(lda_directory_path+'/output/K_'+str(topic_num))
    except FileExistsError:
        pass
    command+=' -o '+lda_directory_path+'/output/K_'+str(topic_num)+'/'
    return command


def parse_result(theta_estimated_filename, phi_estimated_filename, word_list_filename):
    theta=pd.read_csv(theta_estimated_filename, header=None)
    phi=pd.read_csv(phi_estimated_filename, header=None)
    phi.columns=phi.columns.astype('str')
    word_list=pd.read_csv(word_list_filename, header=None)
    word_list[0]=word_list[0].astype('str')
    phi.columns=word_list.iloc[:, 0]
    return theta, phi, word_list


class PerplexityCalculator():
    def __init__(self, BOW_filename, k_min, k_max, kstep, rate_of_train_data, iteration):
        if rate_of_train_data>1:
            exit()
        self.BOW=read_BOW_file(BOW_filename)
        self.k_min=int(k_min)
        self.k_max=int(k_max)
        self.kstep=int(kstep)
        self.rate_of_train_data=float(rate_of_train_data)
        self.iteration=int(iteration)

    # TODO:訓練perplexityも計算
    def split_train_test(self):
        train_BOW=[]
        test_BOW=[]
        for d, document in enumerate(self.BOW):
            all_index=list(range(len(self.BOW[d])))
            train_index=random.sample(all_index, int(len(self.BOW[d])*self.rate_of_train_data))
            test_index=list(set(all_index)-set(train_index))
            train_words=[self.BOW[d][i] for i in train_index]
            test_words=[self.BOW[d][i] for i in test_index]
            train_BOW.append(train_words)
            test_BOW.append(test_words)
        write_BOW_file(train_BOW, lda_directory_path+'/data/trainBOW')
        self.test_BOW=test_BOW

    def execute_estimation(self):
        for i in range(self.k_min, self.k_max, self.kstep):
            command=create_command_string(i, self.iteration)
            subprocess.call(command.split(' '), stdout=subprocess.DEVNULL, stderr=subprocess.STDOUT)

    def calculate_perplexity(self):
        test_words_num=0
        for document in self.test_BOW:
            test_words_num+=len(document)
        for topic_num in range(self.k_min, self.k_max, self.kstep):
            theta, phi, word_list=parse_result(
                lda_directory_path+'/output/K_'+str(topic_num)+'/theta.csv',
                lda_directory_path+'/output/K_'+str(topic_num)+'/phi.csv',
                lda_directory_path+'/output/K_'+str(topic_num)+'/wordList.csv'
            )
            log_likelihood=0.0
            word_prob_df=theta.dot(phi)
            for d, document in enumerate(self.test_BOW):
                document_perplexity=0.0
                for i, word in enumerate(document):
                    if not (word in phi.columns.tolist()):
                        continue
                    document_perplexity+=np.log(word_prob_df.ix[d, word])
                log_likelihood+=document_perplexity
            perplexity=np.exp(-log_likelihood/test_words_num)
            print('topic{0}:{1}'.format(topic_num, perplexity))


    def run(self):
        self.split_train_test()
        self.execute_estimation()
        self.calculate_perplexity()


if __name__=='__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument("BOW_filename", help="Bag Of Words filename")
    parser.add_argument("k_min", help="minimum number of topics for calculation of perplexity", type=int)
    parser.add_argument("k_max", help="max number of topics for calculation of perplexity", type=int)
    parser.add_argument("-s", "--kstep", help="step number of topics for calculation of perplexity", type=int, default=1)
    parser.add_argument("-r", "--rate", help="rate of train data", default=0.7)
    parser.add_argument("-i", "--iteration", help="number of iterations", default=100)
    args = parser.parse_args()
    calculator=PerplexityCalculator(args.BOW_filename, args.k_min, args.k_max, args.kstep, args.rate, args.iteration)
    calculator.run()
