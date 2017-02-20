//
// CollapsedGibbsSamplerFromATM.hpp
//
// Copyright (c) 2017 Shion Hosoda
//
// This software is released under the MIT License.
// http://opensource.org/licenses/mit-license.php
//

#include"GibbsSamplerFromATM.hpp"

using namespace std;

class CollapsedGibbsSamplerFromATM:public GibbsSamplerFromATM{
public:
    CollapsedGibbsSamplerFromATM(const vector<vector<unsigned int> > &bagOfWordsNum,const vector<vector<unsigned int> > &bagOfAuthorsNum,const vector<string> &wordList,const vector<string> &authorList,const unsigned int K,const unsigned int V,const unsigned int A,vector<unsigned int> nThIterationSample);
    void zUpdate();
    virtual void runIteraion();
};
