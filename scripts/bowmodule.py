#
# bow_module.py
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


def read_BOW_file(BOW_filename):
    f = open(BOW_filename, 'r')
    full_string = f.read()
    f.close()
    BOW = []
    for line in full_string.split('\n'):
        document = []
        for word in line.split(','):
            if word != '':
                document.append(word)
        if document != []:
            BOW.append(document)
    return BOW


def write_BOW_file(BOW, BOW_filename):
    f = open(BOW_filename, 'w')
    for d, line in enumerate(BOW):
        for i, word in enumerate(line):
            f.write(word)
            if i != (len(line)-1):
                f.write(',')
        f.write('\n')


def make_bag_of_words_num(input_filename):
    frequency_matrix = pd.read_csv(input_filename)
    frequency_matrix.index = frequency_matrix.iloc[:,0]
    del frequency_matrix[frequency_matrix.columns[0]]
    word_list = frequency_matrix.columns.tolist()
    BOW = []
    for d, document in enumerate(frequency_matrix.index.tolist()):
        bufList = []
        for v, word in enumerate(frequency_matrix.columns.tolist()):
            for i in range(frequency_matrix.loc[document, word]):
                bufList.append(v)
        BOW.append(bufList)
    return BOW, word_list


def make_frequency_matrix(BOW, word_list, output_filename):
    frequency_matrix = []
    for document in BOW:
        frequency_list = []
        for i, word in enumerate(word_list):
            frequency = document.count(i)
            frequency_list.append(frequency)
        frequency_matrix.append(frequency_list)
    fm_df = pd.DataFrame(frequency_matrix, columns=word_list)
    fm_df.to_csv(output_filename)
