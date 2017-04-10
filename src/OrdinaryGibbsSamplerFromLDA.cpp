//
// OrdinaryGibbsSamplerFromLDA.cpp
//
// Copyright (c) 2017 Shion Hosoda
//
// This software is released under the MIT License.
// http://opensource.org/licenses/mit-license.php
//

#include"include/OrdinaryGibbsSamplerFromLDA.hpp"

using namespace std;

OrdinaryGibbsSamplerFromLDA::OrdinaryGibbsSamplerFromLDA(const vector<vector<unsigned int > > &bagOfWordsNum,const vector<string > &wordList,const unsigned int K,const unsigned int V,const vector<unsigned int> nThIterationSample):GibbsSamplerFromLDA(bagOfWordsNum,wordList,K,V,nThIterationSample){//{{{
}//}}}

void OrdinaryGibbsSamplerFromLDA::zUpdate(){//{{{
    random_device rnd;
    mt19937 mt(rnd());
    uniform_real_distribution<double> randN(0,1);
    for(int d=0;d<_model->z.size();d++){
        for(int i=0;i<_model->z[d].size();i++){
            double A=0;
            double cum_prob=0;
            for(int k=0;k<_model->K;k++){
                A+=_model->theta[d][k]*_model->phi[k][_model->bagOfWordsNum[d][i]];
            }
            double random_value=randN(mt);
            for(int k=0;k<_model->K;k++){
                cum_prob+=_model->theta[d][k]*_model->phi[k][_model->bagOfWordsNum[d][i]]/A;
                if(cum_prob>random_value || k==_model->K-1){
                    _model->z[d][i]=k;
                    break;
                }
            }
        }
    }
}//}}}

vector<double> OrdinaryGibbsSamplerFromLDA::samplingDirichlet(vector<unsigned int> param,unsigned int dim,double random_value){//{{{
    int i;
    random_device rnd;
    vector<double> sample;
    for(i=0;i<dim;i++){
        gamma_distribution<double> gamma(param[i],1);
        sample.push_back(gamma(rnd));
    }
    double sampleSum=accumulate(sample.begin(),sample.end(),0.0);
    for(i=0;i<sample.size();i++){
        sample[i]/=sampleSum;
    }

    return sample;
}//}}}

void OrdinaryGibbsSamplerFromLDA::thetaUpdate(){//{{{
    vector<unsigned int> param;

    random_device rnd;
    mt19937 mt(rnd());
    uniform_real_distribution<double> randN(0,1);
    for(int d=0;d<_model->z.size();d++){//D
        for(int k=0;k<_model->K;k++){
            param.push_back(_model->alpha[k]+_model->topicCountsInDoc[d][k]);
        }
        double random_value=randN(mt);
        _model->theta[d]=samplingDirichlet(param,_model->K,random_value);
        param.clear();
    }
}//}}}

void OrdinaryGibbsSamplerFromLDA::phiUpdate(){//{{{
    vector<unsigned int> param;

    random_device rnd;
    mt19937 mt(rnd());
    uniform_real_distribution<double> randN(0,1);
    for(int k=0;k<_model->K;k++){//_model->K
        for(int v=0;v<_model->V;v++){
            param.push_back(_model->beta[v]+_model->vocaCountsFromTopic[k][v]);
        }
        double random_value=randN(mt);
        _model->phi[k]=samplingDirichlet(param,_model->V,random_value);
        param.clear();
    }

}//}}}

void OrdinaryGibbsSamplerFromLDA::runIteraion(){//{{{
    _model->nUpdate();
    thetaUpdate();
    phiUpdate();
    zUpdate();
    _iteration++;
    for(int i=0;i<_nThIterationSample.size();i++){
        if(_iteration==_nThIterationSample[i]){
            samplingParameter();
            break;
        }
    }
    hyperParamUpdate();
}//}}}
