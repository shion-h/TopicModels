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
from iomodule import create_command_string
from bowmodule import make_bag_of_words_num, make_frequency_matrix
from calcldappl import calculate_log_likelihood
lda_directory_path = os.path.split(os.path.abspath(os.path.dirname(__file__)))[0]
timestamp = dt.now().strftime('%m%d_%a_%H_%M_%S')


def calculate_log_prior_value(theta, phi, alpha, beta):
    theta = np.array(theta)
    phi = np.array(phi)
    log_prior_value = 0.0
    for doc in theta:
        log_prior_value += dirichlet.logpdf(x=doc, alpha=alpha)
    for topic in phi:
        log_prior_value += dirichlet.logpdf(x=topic, alpha=beta)
    return log_prior_value


def calculate_target_distribution_value(log_likelihood, log_prior_value, n):
    inverse_temperature = 1.0 / np.log(n)
    target_distribution_value = inverse_temperature * log_likelihood + log_prior_value
    return target_distribution_value


def sampling_param_from_dirichlet(document_num, topic_num, alpha, beta):
    theta = []
    for d in range(document_num):
        theta.append(np.random.dirichlet(alpha))
    phi = []
    for k in range(topic_num):
        phi.append(np.random.dirichlet(beta))
    theta = pd.DataFrame(theta)
    phi = pd.DataFrame(phi)
    return theta, phi


def run_WBIC_metropolis(BOW, alpha, beta, n=0):
    if n == 0:
        for document in BOW:
            n += len(document)
    iteration = 1000
    burn_in = 800
    sampling_interval = 1
    sampling_times = (iteration - burn_in) // sampling_interval
    acceptance_times = 0
    log_likelihood_average = 0
    theta, phi = sampling_param_from_dirichlet(len(BOW), len(alpha), alpha, beta)
    log_likelihood = calculate_log_likelihood(theta, phi, BOW)
    log_prior_value = calculate_log_prior_value(theta, phi, alpha, beta)
    target_dist_value = calculate_target_distribution_value(
                         log_likelihood, log_prior_value, n)
    for i in range(iteration):
        theta_candidate, phi_candidate = sampling_param_from_dirichlet(
                                     len(BOW), len(alpha), alpha, beta)
        log_likelihood_candidate = calculate_log_likelihood(
                            theta_candidate, phi_candidate, BOW)
        log_prior_value_candidate = calculate_log_prior_value(
                            theta_candidate, phi_candidate, alpha, beta)
        target_dist_value_candidate = calculate_target_distribution_value(
                            log_likelihood_candidate, log_prior_value_candidate, n)
        acceptance_prob = np.exp(target_dist_value_candidate - target_dist_value +
                                 log_prior_value - log_prior_value_candidate)
        # print(target_dist_value_candidate, target_dist_value, log_prior_value, log_prior_value_candidate)
        # print(acceptance_prob)
        if acceptance_prob > np.random.uniform():
            acceptance_times += 1
            theta = theta_candidate
            phi = phi_candidate
            log_likelihood = log_likelihood_candidate
            log_prior_value = log_prior_value_candidate
            target_dist_value = target_dist_value_candidate
        if (i > burn_in-1) & (i-burn_in+1)%sampling_interval == 0:
            log_likelihood_average += log_likelihood
    log_likelihood_average /= sampling_times
    print('acceptance_times: ' + str(acceptance_times))
    return log_likelihood_average


class WBICCalculator():
    def __init__(self, BOW_filename, k_min, k_max, kstep, iteration):
        self.BOW_filename = BOW_filename
        self.k_min = int(k_min)
        self.k_max = int(k_max)
        self.kstep = int(kstep)
        self.iteration = int(iteration)

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
                                            lda_directory_path, output_dir,
                                            self.BOW_filename)
            subprocess.call(command.split(' '),
                            stdout=subprocess.DEVNULL, stderr=subprocess.STDOUT)

    def calculate_WBIC_each_topic_num(self):
        BOW, _ = make_bag_of_words_num(self.BOW_filename)
        words_num = 0
        for document in BOW:
            words_num += len(document)
        for topic_num in range(self.k_min, self.k_max, self.kstep):
            alpha = np.loadtxt(lda_directory_path+'/output/'+
                    timestamp+'/K_'+str(topic_num)+'/alpha.csv')
            beta = np.loadtxt(lda_directory_path+'/output/'+
                    timestamp+'/K_'+str(topic_num)+'/beta.csv')
            WBIC = run_WBIC_metropolis(BOW, alpha, beta, words_num)
            print('{0},{1}'.format(topic_num, WBIC))

    def run(self):
        self.execute_estimation()
        print('\nestimate done')
        self.calculate_WBIC_each_topic_num()


if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument("BOW_filename", help="Bag Of Words filename")
    parser.add_argument("k_min", help="minimum number of topics for calculation of perplexity", type=int)
    parser.add_argument("k_max", help="max number of topics for calculation of perplexity", type=int)
    parser.add_argument("-s", "--kstep", help="step number of topics for calculation of perplexity", type=int, default=1)
    parser.add_argument("-i", "--iteration", help="number of iterations(>10)", default=100)
    args = parser.parse_args()
    calculator = WBICCalculator(args.BOW_filename, args.k_min, args.k_max, args.kstep, args.iteration)
    calculator.run()
