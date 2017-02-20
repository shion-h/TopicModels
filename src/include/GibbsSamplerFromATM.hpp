//
// GibbsSamplerFromATM.hpp
//
// Copyright (c) 2017 Shion Hosoda
//
// This software is released under the MIT License.
// http://opensource.org/licenses/mit-license.php
//

#include"GibbsSamplerFromLDA.hpp"

using namespace std;

struct ATM:public LDA{
    const vector<vector<unsigned int> > &bagOfAuthorsNum;
    const unsigned int A;
    vector<vector<unsigned int > > y;
    vector<unsigned int > wordCountsInAuthor;
    vector<vector<unsigned int > > topicCountsInAuthor;

    ATM(const vector<vector<unsigned int> > &bagOfWordsNum,const vector<vector<unsigned int> > &bagOfAuthorsNum,const unsigned int K,const unsigned int V,const unsigned int A);
    virtual void initializeParam();
    virtual void calculateEx();
    virtual void nUpdate();
};

class GibbsSamplerFromATM:public BaseGibbsSampler<ATM>{
    const vector<string> &_authorList;
public:
    GibbsSamplerFromATM(const vector<vector<unsigned int > > &bagOfWordsNum,const vector<vector<unsigned int> > &bagOfAuthorsNum,const vector<string> &wordList,const vector<string> &authorList,const unsigned int K,const unsigned int V,unsigned int A,vector<unsigned int> nThIterationSample);
    virtual ~GibbsSamplerFromATM()=0;
    virtual void yUpdate();
    virtual void hyperParamUpdate();
    virtual void printTopFactor(int numOfTopFactor)const;
    virtual void runIteraion()=0;
};

