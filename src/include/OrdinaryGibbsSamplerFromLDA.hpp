//
// OrdinaryGibbsSamplerFromLDA.hpp
//
// Copyright (c) 2017 Shion Hosoda
//
// This software is released under the MIT License.
// http://opensource.org/licenses/mit-license.php
//

#include"GibbsSamplerFromLDA.hpp"

class OrdinaryGibbsSamplerFromLDA:public BaseGibbsSampler<LDA>{
private:
    vector<vector<double> > theta,phi;
public:
    OrdinaryGibbsSamplerFromLDA(const vector<vector<unsigned int > > &bagOfWordsNum,const vector<string > &wordList,const unsigned int K,const unsigned int V,const vector<unsigned int> nThIterationSample);
    void zUpdate();
    vector<double> samplingDirichlet(vector<unsigned int> param,unsigned int dim,double random_value);
    void thetaUpdate();
    void phiUpdate();
    virtual void runIteraion();
};
