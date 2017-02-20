//
// CollapsedGibbsSamplerFromLDA.hpp
//
// Copyright (c) 2017 Shion Hosoda
//
// This software is released under the MIT License.
// http://opensource.org/licenses/mit-license.php
//

#include"GibbsSamplerFromLDA.hpp"

using namespace std;

class CollapsedGibbsSamplerFromLDA:public BaseGibbsSampler<LDA>{
public:
    CollapsedGibbsSamplerFromLDA(const vector<vector<unsigned int> > &bagOfWordsNum,const vector<string > &wordList,const unsigned int K,const unsigned int V,const vector<unsigned int> nThIterationSample);
    void zUpdate();
    virtual void runIteraion();
};
