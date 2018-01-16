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
    this->initializeParam();
    this->initializeHyperParam();
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
    for(int d=0;d<_bagOfWordsNum.size();d++){
        _nd[d] = _bagOfWordsNum[d].size();
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
    vector<double> ndkSum(_bagOfWordsNum.size());
    for(int d=0;d<_ndk.size();d++){
        ndkSum[d] = 0;
        for(int k=0;k<_ndk[d].size();k++){
            ndkSum[d] += _ndk[d][k];
        }
    }
    for(int d=0;d<_ndk.size();d++){
        for(int k=0;k<_ndk[d].size();k++){
            _ndk[d][k] /= ndkSum[d]/_nd[d];
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
    _alphaTimeSeries.push_back(_alpha);
    _betaTimeSeries.push_back(_beta);
    this->calculateHyperParamSum();
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

vector<double> VariationalBayesEstimatorOnLDA::calculateQz(unsigned int d, unsigned int i)const{//{{{
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

void VariationalBayesEstimatorOnLDA::updateNEx(){//{{{
    vector<vector<double> > nkvBuf(_nkv.size()), ndkBuf(_ndk.size());
    for(int k=0;k<_nkv.size();k++)nkvBuf[k].assign(_V, 0);
    for(int d=0;d<_ndk.size();d++)ndkBuf[d].assign(_K, 0);
    _variationalLowerBoundOfQz = 0;
    for(int d=0;d<_bagOfWordsNum.size();d++){
        for(int i=0;i<_bagOfWordsNum[d].size();i++){
            vector<double> qzdi;
            qzdi = this->calculateQz(d, i);
            double temp = 0;
            for(int k=0; k<_K; k++){
                ndkBuf[d][k] += qzdi[k];
                nkvBuf[k][_bagOfWordsNum[d][i]] += qzdi[k];
                _variationalLowerBoundOfQz += qzdi[k] * log(qzdi[k]);
            }
        }
    }
    for(int d=0;d<_bagOfWordsNum.size();d++){
        for(int k=0;k<_K;k++){
            _ndk[d][k] = ndkBuf[d][k];
        }
    }
    for(int k=0;k<_K;k++){
        _nk[k] = 0;
        for(int v=0;v<_V;v++){
            _nkv[k][v] = nkvBuf[k][v];
            _nk[k] += nkvBuf[k][v];
        }
    }

}//}}}

void VariationalBayesEstimatorOnLDA::updateBeta(unsigned int isAsymmetry){
    double numerator=0, denominator=0;
    if(isAsymmetry == 0){
        double commonBeta;
        numerator=0;
        denominator=0;
        for(int k=0;k<_K;k++){
            for(int v=0;v<_V;v++){
                numerator+=(boost::math::digamma(_nkv[k][v]+_beta[v])-boost::math::digamma(_beta[v]))*_beta[v];
            }
            denominator+=boost::math::digamma(_nk[k]+_betaSum)-boost::math::digamma(_betaSum);
        }
        commonBeta = numerator/denominator/_V;
        _beta.assign(_V, commonBeta);
    }else{
        double _bupdateEtad;
        numerator=0;
        denominator=0;
        for(int v=0;v<_V;v++){
            for(int k=0;k<_K;k++){
                numerator+=(boost::math::digamma(_nkv[k][v]+_beta[v])-boost::math::digamma(_beta[v]))*_beta[v];
                denominator+=boost::math::digamma(_nk[k]+_betaSum)-boost::math::digamma(_betaSum);
            }
            _beta[v] = numerator/denominator;
        }
    }
    _betaTimeSeries.push_back(_beta);
}

void VariationalBayesEstimatorOnLDA::updateHyperParameters(){//{{{
    double numerator=0, denominator=0;
    for(int k=0;k<_K;k++){
        for(int d=0;d<_bagOfWordsNum.size();d++){
            numerator+=(boost::math::digamma(_ndk[d][k]+_alpha[k])-boost::math::digamma(_alpha[k]))*_alpha[k];
            denominator+=boost::math::digamma(_nd[d]+_alphaSum)-boost::math::digamma(_alphaSum);
        }
        _alpha[k]=numerator/denominator;
    }
    _alphaTimeSeries.push_back(_alpha);
    this->updateBeta(1);
    this->calculateHyperParamSum();
}//}}}

double VariationalBayesEstimatorOnLDA::calculateVariationalLowerBound()const{//{{{
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
    vector<vector<int> >  topFactorOfThetaIndex;
    vector<vector<double> > topFactorOfPhi;
    vector<vector<int> >  topFactorOfPhiIndex;
    vector<double> maxValue;
    vector<int> maxValueIndex;

    for(int d=0;d<_thetaEx.size();d++){
        for(int k=0;k<_thetaEx[d].size();k++){
            if(k<numOfTopFactor){
                maxValue.push_back(_thetaEx[d][k]);
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
                if(minInMaxValue<_thetaEx[d][k]){
                    // cout<<"min:"<<minInMaxValue<<"_thetaEx:"<<_thetaEx[d][k]<<endl;
                    maxValue[minInMaxValueIndex]=_thetaEx[d][k];
                    maxValueIndex[minInMaxValueIndex]=k;
                }
            }
        }
        topFactorOfTheta.push_back(maxValue);
        topFactorOfThetaIndex.push_back(maxValueIndex);
        maxValue.clear();
        maxValueIndex.clear();
    }
    for(int k=0;k<_phiEx.size();k++){
        for(int v=0;v<_phiEx[k].size();v++){
            if(v<numOfTopFactor){
                maxValue.push_back(_phiEx[k][v]);
                maxValueIndex.push_back(v);
            }else{
                double minInMaxValue=maxValue[0];
                int minInMaxValueIndex=0;
                for(int i=1;i<maxValue.size();i++){//fi_nd min in maxValue
                    if(minInMaxValue>maxValue[i]){
                        minInMaxValue=maxValue[i];
                        minInMaxValueIndex=i;
                    }
                }
                if(minInMaxValue<_phiEx[k][v]){
                    // cout<<"min:"<<minInMaxValue<<"_phiEx:"<<_phiEx[k][v]<<endl;
                    maxValue[minInMaxValueIndex]=_phiEx[k][v];
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
    outputVector(_thetaEx, thetaFilename);
    outputVector(_phiEx, phiFilename);
    outputVector(_alphaTimeSeries, alphaFilename);
    outputVector(_betaTimeSeries, betaFilename);
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
    double prevVariationalLowerBound = 0;
    double thisVariationalLowerBound = 0;
    unsigned int count = 0;
    while(1){
        prevVariationalLowerBound = thisVariationalLowerBound;
        this->updateNEx();
        thisVariationalLowerBound = this->calculateVariationalLowerBound();
        cout<<"VLB"<<thisVariationalLowerBound<<endl;
        _VLBTimeSeries.push_back(thisVariationalLowerBound);
        count++;
        if(count<2){
        }else if((thisVariationalLowerBound - prevVariationalLowerBound) / abs(thisVariationalLowerBound) < _convergenceDiterminationRate){
            _variationalLowerBound = thisVariationalLowerBound;
            break;
        }
        this->updateHyperParameters();
    }
    this->calculateEx();
    cout<<"iter:"<<count<<endl;
}//}}}
