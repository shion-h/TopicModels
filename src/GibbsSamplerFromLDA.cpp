//
// GibbsSamplerFromLDA.cpp
//
// Copyright (c) 2017 Shion Hosoda
//
// This software is released under the MIT License.
// http://opensource.org/licenses/mit-license.php
//

#include"include/GibbsSamplerFromLDA.hpp"

using namespace std;

LDA::LDA(const vector<vector<unsigned int> > &bagOfWordsNum,const unsigned int K,const unsigned int V)//{{{
    :bagOfWordsNum(bagOfWordsNum),
     K(K),
     V(V),
     z(bagOfWordsNum.size()),
     theta(bagOfWordsNum.size()),
     phi(K),
     thetaEx(bagOfWordsNum.size()),
     phiEx(K),
     alphaSum(0),
     betaSum(0),
     topicCounts(K),
     topicCountsInDoc(bagOfWordsNum.size()),
     vocaCountsFromTopic(K){
    initializeParam();
    initializeHyperParam();
}//}}}

void LDA::initializeParam(){//{{{
    random_device rnd;
    mt19937 mt(rnd());
    uniform_int_distribution<int> randN(0,K-1);

    // for(int d=0;d<z.size();d++)z[d].assign(bagOfWordsNum[d].size(),randN(mt));
    for(int d=0;d<bagOfWordsNum.size();d++){
        vector<unsigned int> ibuf;
        for(int i=0;i<bagOfWordsNum[d].size();i++){
            ibuf.push_back(randN(mt));
        }
        z[d]=ibuf;
        ibuf.clear();
    }
    for(int d=0;d<theta.size();d++)theta[d].assign(K,0);
    for(int k=0;k<phi.size();k++)phi[k].assign(V,0);
    for(int d=0;d<thetaEx.size();d++)thetaEx[d].assign(K,0);
    for(int k=0;k<phiEx.size();k++)phiEx[k].assign(V,0);
    topicCounts.assign(K,0);
    for(int d=0;d<topicCountsInDoc.size();d++)topicCountsInDoc[d].assign(K,0);
    for(int k=0;k<vocaCountsFromTopic.size();k++)vocaCountsFromTopic[k].assign(V,0);
}//}}}

void LDA::initializeHyperParam(){//{{{
    alpha.assign(K,0.1);
    beta.assign(V,0.05);
    calculateHyperParamSum();
}//}}}

void LDA::calculateEx(){//{{{
    for(int d=0;d<z.size();d++){
        for(int k=0;k<K;k++){
            thetaEx[d][k]=(topicCountsInDoc[d][k]+alpha[k])/(z[d].size()+alphaSum);
        }
    }
    for(int k=0;k<K;k++){
        for(int v=0;v<V;v++){
            phiEx[k][v]=(vocaCountsFromTopic[k][v]+beta[v])/(topicCounts[k]+betaSum);
        }
    }
}//}}}

void LDA::calculateHyperParamSum(){//{{{
    alphaSum=0;
    betaSum=0;
    alphaSum=accumulate(alpha.begin(),alpha.end(),0.0);
    betaSum=accumulate(beta.begin(),beta.end(),0.0);
}//}}}

void LDA::nUpdate(){//{{{
    topicCounts.assign(K,0);
    for(int d=0;d<topicCountsInDoc.size();d++)topicCountsInDoc[d].assign(K,0);
    for(int k=0;k<vocaCountsFromTopic.size();k++)vocaCountsFromTopic[k].assign(V,0);
    for(int d=0;d<z.size();d++){
        for(int i=0;i<z[d].size();i++){
            topicCountsInDoc[d][z[d][i]]++;
        }
    }
    for(int d=0;d<z.size();d++){
        for(int i=0;i<z[d].size();i++){
            vocaCountsFromTopic[z[d][i]][bagOfWordsNum[d][i]]++;
        }
    }
    for(int k=0;k<K;k++){
        for(int v=0;v<V;v++){
            topicCounts[k]+=vocaCountsFromTopic[k][v];
        }
    }

}//}}}
