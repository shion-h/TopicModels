//
// LDAWBICMetropolisSampler.cpp
//
// Copyright (c) 2017 Shion Hosoda
//
// This software is released under the MIT License.
// http://opensource.org/licenses/mit-license.php
//

#include "include/LDAWBICMetropolisSampler.hpp"
int main(int argc, char *argv[]){
    BOWFileParser parser(argv[1]);
    parser.readBOWFile();
    parser.makeBagOfWordsNum();
    vector<vector<unsigned int> > BOW = parser.getBagOfWordsNum();
    vector<double> alpha = readHyperParam(argv[2]);
    vector<double> beta = readHyperParam(argv[3]);
    unsigned int n;
    if(argc == 5){
        n = atoi(argv[4]);
    }
    cout<<runWBICMetropolis(BOW, alpha, beta)<<endl;
    return 0;
}
