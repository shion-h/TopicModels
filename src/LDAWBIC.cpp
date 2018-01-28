//
// LDAWBICMetropolisSampler.cpp
//
// Copyright (c) 2017 Shion Hosoda
//
// This software is released under the MIT License.
// http://opensource.org/licenses/mit-license.php
//

#include "include/LDAWBICMetropolisSampler.hpp"

using namespace std;

int main(int argc, char *argv[]){
    BOWFileParser parser(argv[1]);
    vector<double> alpha = readHyperParam(argv[2]);
    vector<double> beta = readHyperParam(argv[3]);
    unsigned int n;
    if(argc == 5){
        n = atoi(argv[4]);
    }
    cout<<alpha.size()<<','<<runWBICMetropolis(parser.getFrequencyMatrix(), parser.getDocVoca(), alpha, beta)<<endl;
    return 0;
}
