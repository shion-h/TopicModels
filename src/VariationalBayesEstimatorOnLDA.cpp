//
// VariationalBayesEstimatorOnLDA.cpp
//
// Copyright (c) 2017 Shion Hosoda
//
// This software is released under the MIT License.
// http://opensource.org/licenses/mit-license.php
//

#include"include/VariationalBayesEstimatorOnLDA.hpp"

using namespace std;

VariationalBayesEstimatorOnLDA::VariationalBayesEstimatorOnLDA(const vector<vector<unsigned int> > &bagOfWordsNum,const vector<string> &wordList,const unsigned int K,const unsigned int V)//{{{
    :_wordList(wordList),
     bagOfWordsNum(bagOfWordsNum),
     K(K),
     V(V),
     thetaEx(bagOfWordsNum.size()),
     phiEx(K),
     alphaSum(0),
     betaSum(0),
     nk(K),
     nd(bagOfWordsNum.size()),
     nkv(K),
     ndk(bagOfWordsNum.size()){
    initializeParam();
    initializeHyperParam();
}//}}}

VariationalBayesEstimatorOnLDA::~VariationalBayesEstimatorOnLDA(){//{{{

}//}}}

void VariationalBayesEstimatorOnLDA::initializeParam(){//{{{
    for(int d=0;d<thetaEx.size();d++)thetaEx[d].assign(K,0);
    for(int k=0;k<phiEx.size();k++)phiEx[k].assign(V,0);
    nk.assign(K,0);
    nd.assign(bagOfWordsNum.size(),0);
    for(int k=0;k<nkv.size();k++)nkv[k].assign(V,0);
    for(int d=0;d<ndk.size();d++)ndk[d].assign(K,0);

    random_device rnd;
    mt19937 mt(rnd());
    uniform_real_distribution<double> randN(0,1);
    double wordCount = 0;
    for(int d=0;d<bagOfWordsNum.size();d++){
        wordCount += bagOfWordsNum[d].size();
    }
    for(int k=0;k<nkv.size();k++){
        for(int v=0;v<nkv[k].size();v++){
            double randomValue=randN(mt);
            nkv[k][v] = randomValue;
        }
    }
    double Z = 0;
    for(int k=0;k<nkv.size();k++){
        for(int v=0;v<nkv[k].size();v++){
            Z += nkv[k][v];
        }
    }
    for(int k=0;k<nkv.size();k++){
        for(int v=0;v<nkv[k].size();v++){
            nkv[k][v] /= Z/wordCount;
        }
    }
    for(int d=0;d<ndk.size();d++){
        for(int k=0;k<ndk[d].size();k++){
            double randomValue=randN(mt);
            ndk[d][k] = randomValue;
        }
    }
    Z = 0;
    for(int d=0;d<ndk.size();d++){
        for(int k=0;k<ndk[d].size();k++){
            Z += ndk[d][k];
        }
    }
    for(int d=0;d<ndk.size();d++){
        for(int k=0;k<ndk[d].size();k++){
            ndk[d][k] /= Z/wordCount;
        }
    }
    for(int d=0;d<bagOfWordsNum.size();d++){
        for(int k=0;k<K;k++){
            nd[d] += ndk[d][k];
        }
    }
    for(int k=0;k<K;k++){
        for(int v=0;v<V;v++){
            nk[k] += nkv[k][v];
        }
    }
}//}}}

void VariationalBayesEstimatorOnLDA::initializeHyperParam(){//{{{
    alpha.assign(K,0.011);
    beta.assign(V,0.003);
    calculateHyperParamSum();
}//}}}

void VariationalBayesEstimatorOnLDA::calculateEx(){//{{{
    for(int d=0;d<bagOfWordsNum.size();d++){
        for(int k=0;k<K;k++){
            thetaEx[d][k]=(ndk[d][k]+alpha[k])/(nd[d]+alphaSum);
        }
    }
    for(int k=0;k<K;k++){
        for(int v=0;v<V;v++){
            phiEx[k][v]=(nkv[k][v]+beta[v])/(nk[k]+betaSum);
        }
    }

}//}}}

void VariationalBayesEstimatorOnLDA::calculateHyperParamSum(){//{{{
    alphaSum=0;
    betaSum=0;
    alphaSum=accumulate(alpha.begin(),alpha.end(),0.0);
    betaSum=accumulate(beta.begin(),beta.end(),0.0);
}//}}}

vector<double> VariationalBayesEstimatorOnLDA::calculateQz(unsigned int d, unsigned int i){//{{{
    vector<double> qzdi;
    double Zq = 0;
    for(int k=0; k<K; k++){
        double numerator=0,denominator=0;
        // cout<<boost::math::digamma(nkv[k][bagOfWordsNum[d][i]]+beta[bagOfWordsNum[d][i]])<<endl;
        // cout<<boost::math::digamma(ndk[d][k]+alpha[k])<<endl;
        //
        // cout<<boost::math::digamma(nk[k])+betaSum<<endl;
        // cout<<boost::math::digamma(nd[d]+alphaSum)<<endl<<endl;
        numerator = exp(boost::math::digamma(nkv[k][bagOfWordsNum[d][i]]+beta[bagOfWordsNum[d][i]])) * exp(boost::math::digamma(ndk[d][k]+alpha[k]));
        denominator = exp(boost::math::digamma(nk[k])+betaSum) * exp(boost::math::digamma(nd[d]+alphaSum));
        double constProbability = numerator / denominator;
        Zq += constProbability;
        qzdi.push_back(constProbability);
    }
    for(int k=0; k<K; k++){
        qzdi[k] /= Zq;
        // cout<<qzdi[k]<<' ';
    }
        // cout<<endl;
    return(qzdi);
}//}}}

void VariationalBayesEstimatorOnLDA::nExUpdate(){//{{{
    vector<vector<double> > nkvBuf(nkv.size()), ndkBuf(ndk.size());
    for(int k=0;k<nkv.size();k++)nkvBuf[k].assign(V,0);
    for(int d=0;d<ndk.size();d++)ndkBuf[d].assign(K,0);
    for(int d=0;d<bagOfWordsNum.size();d++){
        for(int i=0;i<bagOfWordsNum[d].size();i++){
            vector<double> qzdi;
// TODO: 同じvなら呼ばない?
            qzdi = calculateQz(d,i);
            double temp = 0;
            for(int k=0; k<K; k++){
                ndkBuf[d][k] += qzdi[k];
                nkvBuf[k][bagOfWordsNum[d][i]] += qzdi[k];
                // cout<<qzdi[k]<<endl;
            }
        }
    }
    for(int d=0;d<bagOfWordsNum.size();d++){
        nd[d] = 0;
        for(int k=0;k<K;k++){
            // cout<<ndkBuf[d][k]<<endl;
            ndk[d][k] = ndkBuf[d][k];
            nd[d] += ndkBuf[d][k];
        }
    }
    for(int k=0;k<K;k++){
        nk[k] = 0;
        for(int v=0;v<V;v++){
            nkv[k][v] = nkvBuf[k][v];
            nk[k] += nkvBuf[k][v];
        }
    }

}//}}}

void VariationalBayesEstimatorOnLDA::hyperParamUpdate(){//{{{
    double numerator=0,denominator=0;
    for(int k=0;k<K;k++){
        for(int d=0;d<bagOfWordsNum.size();d++){
            // cout<<d<<endl;
            // cout<<k<<endl;
            // cout<<ndk[d][k]<<endl;
            // cout<<nd[d]<<endl;
            // cout<<alpha[k]<<endl;
            // cout<<endl;
            numerator+=(boost::math::digamma(ndk[d][k]+alpha[k])-boost::math::digamma(alpha[k]))*alpha[k];
            denominator+=boost::math::digamma(nd[d]+alphaSum)-boost::math::digamma(alphaSum);
        }
        // cout<<numerator<<endl;
        // cout<<denominator<<endl<<endl;
        alpha[k]=numerator/denominator;
    }
    double betaUpdated;
    numerator=0;
    denominator=0;
    for(int k=0;k<K;k++){
        for(int v=0;v<V;v++){
            // cout<<k<<endl;
            // cout<<v<<endl;
            // cout<<nkv[k][v]<<endl;
            // cout<<beta[v]<<endl;
            // cout<<endl;
            numerator+=(boost::math::digamma(nkv[k][v]+beta[v])-boost::math::digamma(beta[v]))*beta[v];
        }
        denominator+=boost::math::digamma(nk[k]+betaSum)-boost::math::digamma(betaSum);
    }
            // cout<<numerator<<endl;
            // cout<<denominator<<endl;
    betaUpdated=numerator/denominator/V;
    beta.assign(V,betaUpdated);
    calculateHyperParamSum();
}//}}}

double VariationalBayesEstimatorOnLDA::calculateVariationalLowerBound(){//{{{
    double term1=0, term2=0;
    for(int k=0; k<K; k++){
        // cout<<betaSum<<endl;
        // cout<<nk[k]+betaSum<<endl;
        // cout<<boost::math::lgamma(betaSum)<<endl;
        // cout<<boost::math::lgamma(nk[k]+betaSum)<<endl;
        term1 += boost::math::lgamma(betaSum) - boost::math::lgamma(nk[k]+betaSum);
        for(int v=0; v<V; v++){
            term1 += boost::math::lgamma(nkv[k][v]+beta[v]) - boost::math::lgamma(beta[v]);
        }
    }
    for(int d=0; d<bagOfWordsNum.size(); d++){
        // cout<<boost::math::tgamma(alphaSum)/boost::math::tgamma(nd[d]+alphaSum)<<endl;
        term2 += boost::math::lgamma(alphaSum) - boost::math::lgamma(nd[d]+alphaSum);
        for(int k=0; k<K; k++){
            term2 += boost::math::lgamma(ndk[d][k]+alpha[k]) - boost::math::lgamma(alpha[k]);
        }
    }
    // cout<<term1<<' '<<term2<<endl;
    double variationalLowerBound = term1 + term2;
    return(variationalLowerBound);
}//}}}

void VariationalBayesEstimatorOnLDA::printHyperParameter()const{//{{{
    cout<<fixed<<setprecision(5);
    cout<<"alpha:"<<endl;
    for(int i=0;i<alpha.size();i++){
        cout<<alpha[i]<<' ';
    }
    cout<<endl;
    cout<<"beta:"<<endl;
    for(int i=0;i<beta.size();i++){
        cout<<beta[i]<<' ';
    }
    cout<<endl;
}//}}}

void VariationalBayesEstimatorOnLDA::printTopFactor(int numOfTopFactor)const{//{{{
    vector<vector<double> > topFactorOfTheta;
    vector<vector<int> >  topFactorOfThetaIndex;
    vector<vector<double> > topFactorOfPhi;
    vector<vector<int> >  topFactorOfPhiIndex;
    vector<double> maxValue;
    vector<int> maxValueIndex;

    for(int d=0;d<thetaEx.size();d++){
        for(int k=0;k<thetaEx[d].size();k++){
            if(k<numOfTopFactor){
                maxValue.push_back(thetaEx[d][k]);
                maxValueIndex.push_back(k);
            }else{
                double minInMaxValue=maxValue[0];
                int minInMaxValueIndex=0;
                for(int i=1;i<maxValue.size();i++){//find min in maxValue
                    if(minInMaxValue>maxValue[i]){
                        minInMaxValue=maxValue[i];
                        minInMaxValueIndex=i;
                    }
                }
                if(minInMaxValue<thetaEx[d][k]){
                    // cout<<"min:"<<minInMaxValue<<"thetaEx:"<<thetaEx[d][k]<<endl;
                    maxValue[minInMaxValueIndex]=thetaEx[d][k];
                    maxValueIndex[minInMaxValueIndex]=k;
                }
            }
        }
        topFactorOfTheta.push_back(maxValue);
        topFactorOfThetaIndex.push_back(maxValueIndex);
        maxValue.clear();
        maxValueIndex.clear();
    }
    for(int k=0;k<phiEx.size();k++){
        for(int v=0;v<phiEx[k].size();v++){
            if(v<numOfTopFactor){
                maxValue.push_back(phiEx[k][v]);
                maxValueIndex.push_back(v);
            }else{
                double minInMaxValue=maxValue[0];
                int minInMaxValueIndex=0;
                for(int i=1;i<maxValue.size();i++){//find min in maxValue
                    if(minInMaxValue>maxValue[i]){
                        minInMaxValue=maxValue[i];
                        minInMaxValueIndex=i;
                    }
                }
                if(minInMaxValue<phiEx[k][v]){
                    // cout<<"min:"<<minInMaxValue<<"phiEx:"<<phiEx[k][v]<<endl;
                    maxValue[minInMaxValueIndex]=phiEx[k][v];
                    maxValueIndex[minInMaxValueIndex]=v;
                }
            }
        }
        topFactorOfPhi.push_back(maxValue);
        topFactorOfPhiIndex.push_back(maxValueIndex);
        maxValue.clear();
        maxValueIndex.clear();
    }
    for(int d=0;d<topFactorOfTheta.size();d++){
        cout<<"doc"<<d<<':'<<endl;
        for(int n=0;n<topFactorOfTheta[d].size();n++){
            cout<<"topic"<<topFactorOfThetaIndex[d][n]<<": "<<topFactorOfTheta[d][n]<<endl;
        }
    }
    for(int k=0;k<topFactorOfPhi.size();k++){
        cout<<"topic"<<k<<':'<<endl;
        for(int n=0;n<topFactorOfPhi[k].size();n++){
            cout<<_wordList[topFactorOfPhiIndex[k][n]]<<": "<<topFactorOfPhi[k][n]<<endl;
        }
    }
    cout<<"alpha"<<endl;
    for(int i=0;i<alpha.size();i++){
        cout<<alpha[i]<<',';
    }
    cout<<endl;

    cout<<"beta"<<endl;
    cout<<beta[0];
    cout<<endl;
}//}}}

void VariationalBayesEstimatorOnLDA::printThetaEx()const{//{{{
    cout<<fixed<<setprecision(5);
    cout<<"thetaEx:"<<endl;
    for(int i=0;i<thetaEx.size();i++){
        for(int j=0;j<thetaEx[i].size();j++){
            cout<<thetaEx[i][j]<<' ';
        }
        cout<<endl;
    }
}//}}}

void VariationalBayesEstimatorOnLDA::printPhiEx()const{//{{{
    cout<<fixed<<setprecision(5);
    cout<<"phiEx:"<<endl;
    for(int i=0;i<phiEx.size();i++){
        for(int j=0;j<phiEx[i].size();j++){
            cout<<phiEx[i][j]<<' ';
        }
        cout<<endl;
    }
}//}}}

void VariationalBayesEstimatorOnLDA::printNum()const{//{{{
    cout<<fixed<<setprecision(5);
    cout<<"numOfOccurencesOfTopic:"<<endl;
    for(int i=0;i<nk.size();i++){
        cout<<nk[i]<<' ';
    }
    cout<<endl;

    cout<<"numOfOccurencesOfTopicInDoc:"<<endl;
    for(int i=0;i<ndk.size();i++){
        for(int j=0;j<ndk[i].size();j++){
            cout<<ndk[i][j]<<' ';
        }
        cout<<endl;
    }

    cout<<"numOfOccurencesOfVocaFromTopic:"<<endl;
    for(int i=0;i<nkv.size();i++){
        for(int j=0;j<nkv[i].size();j++){
            cout<<nkv[i][j]<<' ';
        }
        cout<<endl;
    }
}//}}}

void VariationalBayesEstimatorOnLDA::writeParameter(string thetaFilename,string phiFilename,string alphaFilename,string betaFilename)const{//{{{
// TODO:VLB書き出し
    ofstream thetaOutput;
    ofstream phiOutput;
    ofstream alphaOutput;
    ofstream betaOutput;
    thetaOutput.open(thetaFilename,ios::out);
    phiOutput.open(phiFilename,ios::out);
    alphaOutput.open(alphaFilename,ios::out);
    betaOutput.open(betaFilename,ios::out);

    for(int i=0;i<thetaEx.size();i++){
        for(int j=0;j<thetaEx[i].size();j++){
            thetaOutput<<thetaEx[i][j];
            if(j!=(thetaEx[i].size()-1)){
                thetaOutput<<',';
            }
        }
        thetaOutput<<endl;
    }

    for(int i=0;i<phiEx.size();i++){
        for(int j=0;j<phiEx[i].size();j++){
            phiOutput<<phiEx[i][j];
            if(j!=(phiEx[i].size()-1)){
                phiOutput<<',';
            }
        }
        phiOutput<<endl;
    }
    for(int i=0;i<alpha.size();i++){
        alphaOutput<<alpha[i];
        alphaOutput<<endl;
    }
    for(int i=0;i<beta.size();i++){
        betaOutput<<beta[i];
        betaOutput<<endl;
    }
    thetaOutput.close();
    phiOutput.close();
    alphaOutput.close();
    betaOutput.close();
}//}}}

void VariationalBayesEstimatorOnLDA::writeThetaEx(string thetaFilename)const{//{{{
    ofstream thetaOutput;
    thetaOutput.open(thetaFilename,ios::out);

    for(int i=0;i<thetaEx.size();i++){
        for(int j=0;j<thetaEx[i].size();j++){
            thetaOutput<<thetaEx[i][j];
            if(j!=(thetaEx[i].size()-1)){
                thetaOutput<<',';
            }
        }
        thetaOutput<<endl;
    }
}//}}}

void VariationalBayesEstimatorOnLDA::writePhiEx(string phiFilename)const{//{{{
    ofstream phiOutput;
    phiOutput.open(phiFilename,ios::out);

    for(int i=0;i<phiEx.size();i++){
        for(int j=0;j<phiEx[i].size();j++){
            phiOutput<<phiEx[i][j];
            if(j!=(phiEx[i].size()-1)){
                phiOutput<<',';
            }
        }
        phiOutput<<endl;
    }
}//}}}

void VariationalBayesEstimatorOnLDA::writeVariationalLowerBound(string VLBFilename)const{
    ofstream VLBOutput;
    VLBOutput.open(VLBFilename,ios::out);
    VLBOutput<<_variationalLowerBound<<endl;
    VLBOutput.close();
}

void VariationalBayesEstimatorOnLDA::runIteraions(){//{{{
    double prevVariationalLowerBound = 1;
    double thisVariationalLowerBound = 2;
    unsigned int count = 0;
    while(1){
        if(count<2){
        }else if((thisVariationalLowerBound - prevVariationalLowerBound) / abs(thisVariationalLowerBound) < 0.001){
            _variationalLowerBound = thisVariationalLowerBound;
            break;
        }
        prevVariationalLowerBound = thisVariationalLowerBound;
        calculateEx();
        nExUpdate();
        hyperParamUpdate();
        thisVariationalLowerBound = calculateVariationalLowerBound();
        cout<<"VLB"<<thisVariationalLowerBound<<endl;
        count++;
    }
    cout<<"iter:"<<count<<endl;
}//}}}
