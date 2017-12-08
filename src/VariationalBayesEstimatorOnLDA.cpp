//
// VariationalBayesEstimatorOnLDA.cpp
//
// Copyright (c) 2017 Shion Hosoda
//
// This software is released u_nder the MIT License.
// http://opensource.org/licenses/mit-license.php
//

#include"include/VariationalBayesEstimatorOnLDA.hpp"

using namespace std;

VariationalBayesEstimatorOnLDA::VariationalBayesEstimatorOnLDA(const vector<vector<unsigned int> > &bagOfWordsNum, const vector<string> &wordList, const unsigned int K, const unsigned int V, const double convergenceDiterminationRate)//{{{
    :_wordList(wordList),
     _bagOfWordsNum(bagOfWordsNum),
     _K(K),
     _V(V),
     _convergenceDiterminationRate(convergenceDiterminationRate),
     _thetaEx(bagOfWordsNum.size()),
     _phiEx(K),
     _alphaSum(0),
     _betaSum(0),
     _nk(K),
     _nd(bagOfWordsNum.size()),
     _nkv(K),
     _ndk(bagOfWordsNum.size()){
    initializeParam();
    initializeHyperParam();
}//}}}

VariationalBayesEstimatorOnLDA::~VariationalBayesEstimatorOnLDA(){//{{{

}//}}}

void VariationalBayesEstimatorOnLDA::initializeParam(){//{{{
    for(int d=0;d<_thetaEx.size();d++)_thetaEx[d].assign(_K, 0);
    for(int k=0;k<_phiEx.size();k++)_phiEx[k].assign(_V, 0);
    _nk.assign(_K, 0);
    _nd.assign(_bagOfWordsNum.size(), 0);
    for(int k=0;k<_nkv.size();k++)_nkv[k].assign(_V, 0);
    for(int d=0;d<_ndk.size();d++)_ndk[d].assign(_K, 0);

    random_device rnd;
    mt19937 mt(rnd());
    uniform_real_distribution<double> randN(0, 1);
    double wordCount = 0;
    for(int d=0;d<_bagOfWordsNum.size();d++){
        wordCount += _bagOfWordsNum[d].size();
    }
    for(int k=0;k<_nkv.size();k++){
        for(int v=0;v<_nkv[k].size();v++){
            double randomValue=randN(mt);
            _nkv[k][v] = randomValue;
        }
    }
    double Z = 0;
    for(int k=0;k<_nkv.size();k++){
        for(int v=0;v<_nkv[k].size();v++){
            Z += _nkv[k][v];
        }
    }
    for(int k=0;k<_nkv.size();k++){
        for(int v=0;v<_nkv[k].size();v++){
            _nkv[k][v] /= Z/wordCount;
        }
    }
    for(int d=0;d<_ndk.size();d++){
        for(int k=0;k<_ndk[d].size();k++){
            double randomValue=randN(mt);
            _ndk[d][k] = randomValue;
        }
    }
    Z = 0;
    for(int d=0;d<_ndk.size();d++){
        for(int k=0;k<_ndk[d].size();k++){
            Z += _ndk[d][k];
        }
    }
    for(int d=0;d<_ndk.size();d++){
        for(int k=0;k<_ndk[d].size();k++){
            _ndk[d][k] /= Z/wordCount;
        }
    }
    for(int d=0;d<_bagOfWordsNum.size();d++){
        for(int k=0;k<_K;k++){
            _nd[d] += _ndk[d][k];
        }
    }
    for(int k=0;k<_K;k++){
        for(int v=0;v<_V;v++){
            _nk[k] += _nkv[k][v];
        }
    }
}//}}}

void VariationalBayesEstimatorOnLDA::initializeHyperParam(){//{{{
    _alpha.assign(_K, 0.1);
    _beta.assign(_V, 0.05);
    calculateHyperParamSum();
}//}}}

void VariationalBayesEstimatorOnLDA::calculateEx(){//{{{
    for(int d=0;d<_bagOfWordsNum.size();d++){
        for(int k=0;k<_K;k++){
            _thetaEx[d][k]=(_ndk[d][k]+_alpha[k])/(_nd[d]+_alphaSum);
        }
    }
    for(int k=0;k<_K;k++){
        for(int v=0;v<_V;v++){
            _phiEx[k][v]=(_nkv[k][v]+_beta[v])/(_nk[k]+_betaSum);
        }
    }

}//}}}

void VariationalBayesEstimatorOnLDA::calculateHyperParamSum(){//{{{
    _alphaSum=0;
    _betaSum=0;
    _alphaSum=accumulate(_alpha.begin(), _alpha.end(), 0.0);
    _betaSum=accumulate(_beta.begin(), _beta.end(), 0.0);
}//}}}

vector<double> VariationalBayesEstimatorOnLDA::calculateQz(unsigned int d, unsigned int i){//{{{
    vector<double> qzdi;
    double Zq = 0;
    for(int k=0; k<_K; k++){
        double numerator=0, denominator=0;
        numerator = exp(boost::math::digamma(_nkv[k][_bagOfWordsNum[d][i]]+_beta[_bagOfWordsNum[d][i]])) * exp(boost::math::digamma(_ndk[d][k]+_alpha[k]));
        denominator = exp(boost::math::digamma(_nk[k])+_betaSum) * exp(boost::math::digamma(_nd[d]+_alphaSum));
        double constProbability = numerator / denominator;
        Zq += constProbability;
        qzdi.push_back(constProbability);
    }
    for(int k=0; k<_K; k++){
        qzdi[k] /= Zq;
    }
    return(qzdi);
}//}}}

void VariationalBayesEstimatorOnLDA::nExUpdate(){//{{{
    vector<vector<double> > _nkvBuf(_nkv.size()), _ndkBuf(_ndk.size());
    for(int k=0;k<_nkv.size();k++)_nkvBuf[k].assign(_V, 0);
    for(int d=0;d<_ndk.size();d++)_ndkBuf[d].assign(_K, 0);
    _variationalLowerBoundOfQz = 0;
    for(int d=0;d<_bagOfWordsNum.size();d++){
        for(int i=0;i<_bagOfWordsNum[d].size();i++){
            vector<double> qzdi;
// TODO: 同じvなら呼ばない?
            qzdi = calculateQz(d, i);
            double temp = 0;
            for(int k=0; k<_K; k++){
                _ndkBuf[d][k] += qzdi[k];
                _nkvBuf[k][_bagOfWordsNum[d][i]] += qzdi[k];
                _variationalLowerBoundOfQz += qzdi[k] * log(qzdi[k]);
                // cout<<qzdi[k]<<endl;
            }
        }
    }
    for(int d=0;d<_bagOfWordsNum.size();d++){
        _nd[d] = 0;
        for(int k=0;k<_K;k++){
            // cout<<_ndkBuf[d][k]<<endl;
            _ndk[d][k] = _ndkBuf[d][k];
            _nd[d] += _ndkBuf[d][k];
        }
    }
    for(int k=0;k<_K;k++){
        _nk[k] = 0;
        for(int v=0;v<_V;v++){
            _nkv[k][v] = _nkvBuf[k][v];
            _nk[k] += _nkvBuf[k][v];
        }
    }

}//}}}

void VariationalBayesEstimatorOnLDA::hyperParamUpdate(){//{{{
    double numerator=0, denominator=0;
    for(int k=0;k<_K;k++){
        for(int d=0;d<_bagOfWordsNum.size();d++){
            // cout<<d<<endl;
            // cout<<k<<endl;
            // cout<<_ndk[d][k]<<endl;
            // cout<<_nd[d]<<endl;
            // cout<<_alpha[k]<<endl;
            // cout<<endl;
            numerator+=(boost::math::digamma(_ndk[d][k]+_alpha[k])-boost::math::digamma(_alpha[k]))*_alpha[k];
            denominator+=boost::math::digamma(_nd[d]+_alphaSum)-boost::math::digamma(_alphaSum);
        }
        // cout<<numerator<<endl;
        // cout<<denominator<<endl<<endl;
        _alpha[k]=numerator/denominator;
    }
    double _betaUpdated;
    numerator=0;
    denominator=0;
    for(int k=0;k<_K;k++){
        for(int v=0;v<_V;v++){
            // cout<<k<<endl;
            // cout<<v<<endl;
            // cout<<_nkv[k][v]<<endl;
            // cout<<_beta[v]<<endl;
            // cout<<endl;
            numerator+=(boost::math::digamma(_nkv[k][v]+_beta[v])-boost::math::digamma(_beta[v]))*_beta[v];
        }
        denominator+=boost::math::digamma(_nk[k]+_betaSum)-boost::math::digamma(_betaSum);
    }
            // cout<<numerator<<endl;
            // cout<<denominator<<endl;
    _betaUpdated=numerator/denominator/_V;
    _beta.assign(_V, _betaUpdated);
    calculateHyperParamSum();
}//}}}

double VariationalBayesEstimatorOnLDA::calculateVariationalLowerBound(){//{{{
    double term1=0, term2=0;
    for(int k=0; k<_K; k++){
        term1 += boost::math::lgamma(_betaSum) - boost::math::lgamma(_nk[k]+_betaSum);
        for(int v=0; v<_V; v++){
            term1 += boost::math::lgamma(_nkv[k][v]+_beta[v]) - boost::math::lgamma(_beta[v]);
        }
    }
    for(int d=0; d<_bagOfWordsNum.size(); d++){
        term2 += boost::math::lgamma(_alphaSum) - boost::math::lgamma(_nd[d]+_alphaSum);
        for(int k=0; k<_K; k++){
            term2 += boost::math::lgamma(_ndk[d][k]+_alpha[k]) - boost::math::lgamma(_alpha[k]);
        }
    }
    double variationalLowerBound = term1 + term2 - _variationalLowerBoundOfQz;
    return(variationalLowerBound);
}//}}}

void VariationalBayesEstimatorOnLDA::printHyperParameter()const{//{{{
    cout<<"_alpha:"<<endl;
    for(int i=0;i<_alpha.size();i++){
        cout<<_alpha[i]<<' ';
    }
    cout<<endl;
    cout<<"_beta:"<<endl;
    for(int i=0;i<_beta.size();i++){
        cout<<_beta[i]<<' ';
    }
    cout<<endl;
}//}}}

void VariationalBayesEstimatorOnLDA::printTopFactor(int numOfTopFactor)const{//{{{
    vector<vector<double> > topFactorOfTheta;
    vector<vector<int> >  topFactorOfThetaI_ndex;
    vector<vector<double> > topFactorOfPhi;
    vector<vector<int> >  topFactorOfPhiI_ndex;
    vector<double> maxValue;
    vector<int> maxValueI_ndex;

    for(int d=0;d<_thetaEx.size();d++){
        for(int k=0;k<_thetaEx[d].size();k++){
            if(k<numOfTopFactor){
                maxValue.push_back(_thetaEx[d][k]);
                maxValueI_ndex.push_back(k);
            }else{
                double minInMaxValue=maxValue[0];
                int minInMaxValueI_ndex=0;
                for(int i=1;i<maxValue.size();i++){//fi_nd min in maxValue
                    if(minInMaxValue>maxValue[i]){
                        minInMaxValue=maxValue[i];
                        minInMaxValueI_ndex=i;
                    }
                }
                if(minInMaxValue<_thetaEx[d][k]){
                    // cout<<"min:"<<minInMaxValue<<"_thetaEx:"<<_thetaEx[d][k]<<endl;
                    maxValue[minInMaxValueI_ndex]=_thetaEx[d][k];
                    maxValueI_ndex[minInMaxValueI_ndex]=k;
                }
            }
        }
        topFactorOfTheta.push_back(maxValue);
        topFactorOfThetaI_ndex.push_back(maxValueI_ndex);
        maxValue.clear();
        maxValueI_ndex.clear();
    }
    for(int k=0;k<_phiEx.size();k++){
        for(int v=0;v<_phiEx[k].size();v++){
            if(v<numOfTopFactor){
                maxValue.push_back(_phiEx[k][v]);
                maxValueI_ndex.push_back(v);
            }else{
                double minInMaxValue=maxValue[0];
                int minInMaxValueI_ndex=0;
                for(int i=1;i<maxValue.size();i++){//fi_nd min in maxValue
                    if(minInMaxValue>maxValue[i]){
                        minInMaxValue=maxValue[i];
                        minInMaxValueI_ndex=i;
                    }
                }
                if(minInMaxValue<_phiEx[k][v]){
                    // cout<<"min:"<<minInMaxValue<<"_phiEx:"<<_phiEx[k][v]<<endl;
                    maxValue[minInMaxValueI_ndex]=_phiEx[k][v];
                    maxValueI_ndex[minInMaxValueI_ndex]=v;
                }
            }
        }
        topFactorOfPhi.push_back(maxValue);
        topFactorOfPhiI_ndex.push_back(maxValueI_ndex);
        maxValue.clear();
        maxValueI_ndex.clear();
    }
    for(int d=0;d<topFactorOfTheta.size();d++){
        cout<<"doc"<<d<<':'<<endl;
        for(int n=0;n<topFactorOfTheta[d].size();n++){
            cout<<"topic"<<topFactorOfThetaI_ndex[d][n]<<": "<<topFactorOfTheta[d][n]<<endl;
        }
    }
    for(int k=0;k<topFactorOfPhi.size();k++){
        cout<<"topic"<<k<<':'<<endl;
        for(int n=0;n<topFactorOfPhi[k].size();n++){
            cout<<_wordList[topFactorOfPhiI_ndex[k][n]]<<": "<<topFactorOfPhi[k][n]<<endl;
        }
    }
    cout<<"_alpha"<<endl;
    for(int i=0;i<_alpha.size();i++){
        cout<<_alpha[i]<<',';
    }
    cout<<endl;

    cout<<"_beta"<<endl;
    cout<<_beta[0];
    cout<<endl;
}//}}}

void VariationalBayesEstimatorOnLDA::printThetaEx()const{//{{{
    cout<<"_thetaEx:"<<endl;
    for(int i=0;i<_thetaEx.size();i++){
        for(int j=0;j<_thetaEx[i].size();j++){
            cout<<_thetaEx[i][j]<<' ';
        }
        cout<<endl;
    }
}//}}}

void VariationalBayesEstimatorOnLDA::printPhiEx()const{//{{{
    cout<<"_phiEx:"<<endl;
    for(int i=0;i<_phiEx.size();i++){
        for(int j=0;j<_phiEx[i].size();j++){
            cout<<_phiEx[i][j]<<' ';
        }
        cout<<endl;
    }
}//}}}

void VariationalBayesEstimatorOnLDA::printNum()const{//{{{
    cout<<"numOfOccurencesOfTopic:"<<endl;
    for(int i=0;i<_nk.size();i++){
        cout<<_nk[i]<<' ';
    }
    cout<<endl;

    cout<<"numOfOccurencesOfTopicInDoc:"<<endl;
    for(int i=0;i<_ndk.size();i++){
        for(int j=0;j<_ndk[i].size();j++){
            cout<<_ndk[i][j]<<' ';
        }
        cout<<endl;
    }

    cout<<"numOfOccurencesOfVocaFromTopic:"<<endl;
    for(int i=0;i<_nkv.size();i++){
        for(int j=0;j<_nkv[i].size();j++){
            cout<<_nkv[i][j]<<' ';
        }
        cout<<endl;
    }
}//}}}

void VariationalBayesEstimatorOnLDA::writeParameter(string thetaFilename, string phiFilename, string alphaFilename, string betaFilename)const{//{{{
    ofstream thetaOutput;
    ofstream phiOutput;
    ofstream alphaOutput;
    ofstream betaOutput;
    thetaOutput.open(thetaFilename, ios::out);
    phiOutput.open(phiFilename, ios::out);
    alphaOutput.open(alphaFilename, ios::out);
    betaOutput.open(betaFilename, ios::out);

    for(int i=0;i<_thetaEx.size();i++){
        for(int j=0;j<_thetaEx[i].size();j++){
            thetaOutput<<_thetaEx[i][j];
            if(j!=(_thetaEx[i].size()-1)){
                thetaOutput<<',';
            }
        }
        thetaOutput<<endl;
    }

    for(int i=0;i<_phiEx.size();i++){
        for(int j=0;j<_phiEx[i].size();j++){
            phiOutput<<_phiEx[i][j];
            if(j!=(_phiEx[i].size()-1)){
                phiOutput<<',';
            }
        }
        phiOutput<<endl;
    }
    for(int i=0;i<_alpha.size();i++){
        alphaOutput<<_alpha[i];
        alphaOutput<<endl;
    }
    for(int i=0;i<_beta.size();i++){
        betaOutput<<_beta[i];
        betaOutput<<endl;
    }
    thetaOutput.close();
    phiOutput.close();
    alphaOutput.close();
    betaOutput.close();
}//}}}

void VariationalBayesEstimatorOnLDA::writeVariationalLowerBound(string VLBFilename, string VLBTimeSeriesFilename)const{//{{{
    ofstream VLBOutput;
    ofstream VLBTimeSeriesOutput;
    VLBOutput.open(VLBFilename, ios::out);
    VLBTimeSeriesOutput.open(VLBTimeSeriesFilename, ios::out);
    VLBOutput<<_variationalLowerBound<<endl;
    for(int i=0;i<_VLBTimeSeries.size();i++){
        VLBTimeSeriesOutput<<_VLBTimeSeries[i];
        VLBTimeSeriesOutput<<endl;
    }
    VLBOutput.close();
    VLBTimeSeriesOutput.close();
}//}}}

void VariationalBayesEstimatorOnLDA::runIteraions(){//{{{
    double prevVariationalLowerBound = 1;
    double thisVariationalLowerBound = 2;
    unsigned int count = 0;
    while(1){
        if(count<2){
        }else if((thisVariationalLowerBound - prevVariationalLowerBound) / abs(thisVariationalLowerBound) < _convergenceDiterminationRate){
            _variationalLowerBound = thisVariationalLowerBound;
            break;
        }
        prevVariationalLowerBound = thisVariationalLowerBound;
        calculateEx();
        nExUpdate();
        hyperParamUpdate();
        thisVariationalLowerBound = calculateVariationalLowerBound();
        cout<<"VLB"<<thisVariationalLowerBound<<endl;
        _VLBTimeSeries.push_back(thisVariationalLowerBound);
        count++;
    }
    cout<<"iter:"<<count<<endl;
}//}}}
