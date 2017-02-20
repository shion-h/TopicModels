//
// CollapsedGibbsSamplerFromATM.cpp
//
// Copyright (c) 2017 Shion Hosoda
//
// This software is released under the MIT License.
// http://opensource.org/licenses/mit-license.php
//

#include"include/CollapsedGibbsSamplerFromATM.hpp"

using namespace std;

CollapsedGibbsSamplerFromATM::CollapsedGibbsSamplerFromATM(const vector<vector<unsigned int> > &bagOfWordsNum,const vector<vector<unsigned int> > &bagOfAuthorsNum,const vector<string> &wordList,const vector<string> &authorList,const unsigned int K,const unsigned int V,const unsigned int A,vector<unsigned int> nThIterationSample)//{{{
    :GibbsSamplerFromATM(bagOfWordsNum,bagOfAuthorsNum,wordList,authorList,K,V,A,nThIterationSample){
}//}}}

void CollapsedGibbsSamplerFromATM::zUpdate(){//{{{
    random_device rnd;
    mt19937 mt(rnd());
    uniform_real_distribution<double> randN(0,1);

    for(int d=0;d<_model->z.size();d++){
        for(int i=0;i<_model->z[d].size();i++){
            double cumProb=0;
            vector<double> topicProb;
            for(int k=0;k<_model->K;k++){
                double nak,nkv,na,nk;
                na=_model->wordCountsInAuthor[_model->y[d][i]]-1;
                if(_model->z[d][i]==k){
                    nk=_model->topicCounts[k]-1;
                    nkv=_model->vocaCountsFromTopic[k][_model->bagOfWordsNum[d][i]]-1;
                    nak=_model->topicCountsInAuthor[_model->y[d][i]][k]-1;
                }else{
                    nk=_model->topicCounts[k];
                    nkv=_model->vocaCountsFromTopic[k][_model->bagOfWordsNum[d][i]];
                    nak=_model->topicCountsInAuthor[_model->y[d][i]][k];
                }
                topicProb.push_back(static_cast<double>((nkv+_model->beta[_model->bagOfWordsNum[d][i]])*(nak+_model->alpha[k]))/static_cast<double>((nk+_model->betaSum)*(na+_model->alphaSum)));
            }
            double sumProb=accumulate(topicProb.begin(),topicProb.end(),0.0);
            double randomValue=randN(mt);
            for(int k=0;k<_model->K;k++){
                cumProb+=topicProb[k]/sumProb;
                if(cumProb>randomValue || k==_model->K-1){
                    _model->z[d][i]=k;
                    break;
                }

            }
        }
    }
}//}}}

void CollapsedGibbsSamplerFromATM::runIteraion(){//{{{
    yUpdate();
    _model->nUpdate();
    zUpdate();
    _iteration++;
    // cout<<_iteration<<endl;
    for(int i=0;i<_nThIterationSample.size();i++){
        if(_iteration==_nThIterationSample[i]){
            samplingParam();
            break;
        }
    }
    hyperParamUpdate();
}//}}}
