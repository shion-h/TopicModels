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

VariationalBayesEstimatorOnLDA::VariationalBayesEstimatorOnLDA(const BOWFileParser &parser, const unsigned int K, const double convergenceDiterminationRate)//{{{
    :_frequencyMatrix(parser.getFrequencyMatrix()),
     _docVoca(parser.getDocVoca()),
     _D(_frequencyMatrix.size()),
     _K(K),
     _V(parser.getV()),
     _qz(_D),
     _convergenceDiterminationRate(convergenceDiterminationRate),
     _thetaEx(_D),
     _phiEx(K),
     _alphaSum(0),
     _betaSum(0),
     _nk(K),
     _nd(_D),
     _nkv(K),
     _ndk(_D){
    this->initializeParam();
    this->initializeHyperParam();
}//}}}

VariationalBayesEstimatorOnLDA::~VariationalBayesEstimatorOnLDA(){//{{{

}//}}}

void VariationalBayesEstimatorOnLDA::initializeParam(){//{{{
    for(int d=0; d<_thetaEx.size(); d++)_thetaEx[d].assign(_K, 0);
    for(int k=0; k<_phiEx.size(); k++)_phiEx[k].assign(_V, 0);
    for(int d=0; d<_ndk.size(); d++)_ndk[d].assign(_K, 0);
    for(int k=0; k<_nkv.size(); k++)_nkv[k].assign(_V, 0);
    _nk.assign(_K, 0);

    random_device rnd;
    mt19937 mt(rnd());
    uniform_real_distribution<double> randN(0, 1);
    for(int d=0; d<_docVoca.size(); d++){
        _qz[d] = vector<vector<double> >(_docVoca[d].size());
        for(int l=0; l<_docVoca[d].size(); l++){
            _qz[d][l] = vector<double>(_K);
            for(int k=0; k<_K; k++){
                double randomValue = randN(mt);
                _qz[d][l][k] = randomValue;
            }
            double qzSum = accumulate(_qz[d][l].begin(), _qz[d][l].end(), 0.0);
            for(int k=0; k<_K; k++){
                _qz[d][l][k] /= qzSum;
            }
        }
    }

    _nd.assign(_D, 0);
    for(int d=0;d<_D;d++){
        _nd[d] = accumulate(_frequencyMatrix[d].begin(), _frequencyMatrix[d].end(), 0.0);
    }

    for(int d=0;d<_ndk.size();d++)_ndk[d].assign(_K, 0);
    for(int d=0;d<_ndk.size();d++){
        for(int k=0;k<_ndk[d].size();k++){
            double randomValue=randN(mt);
            _ndk[d][k] = randomValue;
        }
    }
    vector<double> ndkSum(_D, 0);
    for(int d=0;d<_ndk.size();d++){
        ndkSum[d] = accumulate(_ndk[d].begin(), _ndk[d].end(), 0.0);
    }
    for(int d=0;d<_ndk.size();d++){
        for(int k=0;k<_ndk[d].size();k++){
            _ndk[d][k] /= ndkSum[d]/_nd[d];
        }
    }

    for(int k=0;k<_nkv.size();k++)_nkv[k].assign(_V, 0);
    for(int k=0;k<_nkv.size();k++){
        for(int v=0;v<_nkv[k].size();v++){
            double randomValue=randN(mt);
            _nkv[k][v] = randomValue;
        }
    }
    double Z = 0;
    for(int k=0;k<_nkv.size();k++){
        Z += accumulate(_nkv[k].begin(), _nkv[k].end(), 0.0);
    }
    double wordCount = accumulate(_nd.begin(), _nd.end(), 0.0);
    for(int k=0;k<_nkv.size();k++){
        for(int v=0;v<_nkv[k].size();v++){
            _nkv[k][v] /= Z/wordCount;
        }
    }

    _nk.assign(_K, 0);
    for(int k=0;k<_K;k++){
        _nk[k] = accumulate(_nkv[k].begin(), _nkv[k].end(), 0.0);
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
    for(int d=0;d<_D;d++){
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
    _alphaSum = accumulate(_alpha.begin(), _alpha.end(), 0.0);
    _betaSum = accumulate(_beta.begin(), _beta.end(), 0.0);
}//}}}

void VariationalBayesEstimatorOnLDA::updateQz(){//{{{
    for(int d=0;d<_docVoca.size();d++){
        for(int l=0; l<_docVoca[d].size(); l++){
            double Zq = 0;
            unsigned int v = _docVoca[d][l];
            for(int k=0; k<_K; k++){
                double term1 = 0, term2 = 0, term3 = 0, term4 = 0;
                double xikv = _nkv[k][v] + _beta[v];
                double xidk = _ndk[d][k] + _alpha[k];
                double xik = _nk[k] + _betaSum;
                double xid = _nd[d] + _alphaSum;
                try{
                    term1 = exp(boost::math::digamma(xikv));
                }catch(...){
                    term1 = 0;
                }
                try{
                    term2 = exp(boost::math::digamma(xidk));
                }catch(...){
                    term2 = 0;
                }
                try{
                    term3 = exp(boost::math::digamma(xik));
                }catch(...){
                    term3 = 0;
                }
                try{
                    term4 = exp(boost::math::digamma(xid));
                }catch(...){
                    term4 = 0;
                }
                double constProbability = (term1 * term2) / (term3 * term4);
                Zq += constProbability;
                _qz[d][l][k] = constProbability;
            }
            for(int k=0; k<_K; k++){
                _qz[d][l][k] /= Zq;
            }
        }
    }
}//}}}

void VariationalBayesEstimatorOnLDA::updateNEx(){//{{{
    vector<vector<double> > nkvBuf(_nkv.size()), ndkBuf(_ndk.size());
    for(int k=0;k<_nkv.size();k++)nkvBuf[k].assign(_V, 0);
    for(int d=0;d<_ndk.size();d++)ndkBuf[d].assign(_K, 0);
    for(int d=0;d<_docVoca.size();d++){
        for(int l=0; l<_docVoca[d].size(); l++){
            unsigned int v = _docVoca[d][l];
            for(int k=0; k<_K; k++){
                ndkBuf[d][k] += _qz[d][l][k] * _frequencyMatrix[d][v];
                nkvBuf[k][v] += _qz[d][l][k] * _frequencyMatrix[d][v];
            }
        }
    }
    for(int d=0;d<_D;d++){
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

void VariationalBayesEstimatorOnLDA::updateBeta(BetaUpdateManner manner){//{{{
    if(manner == SYMMETRY){
        double numerator=0;
        double denominator=0;
        for(int k=0;k<_K;k++){
            for(int v=0;v<_V;v++){
                numerator += (boost::math::digamma(_nkv[k][v]+_beta[v]) - boost::math::digamma(_beta[v]))*_beta[v];
            }
            denominator += boost::math::digamma(_nk[k]+_betaSum) - boost::math::digamma(_betaSum);
        }
        double commonBeta = numerator/denominator/_V;
        _beta.assign(_V, commonBeta);
    }else{
        double denominator=0;
        for(int k=0;k<_K;k++){
            denominator += boost::math::digamma(_nk[k]+_betaSum) - boost::math::digamma(_betaSum);
        }
        for(int v=0;v<_V;v++){
            double numerator=0;
            for(int k=0;k<_K;k++){
                numerator += (boost::math::digamma(_nkv[k][v]+_beta[v])-boost::math::digamma(_beta[v]))*_beta[v];
            }
            _beta[v] = numerator/denominator;
        }
    }
    _betaTimeSeries.push_back(_beta);
}//}}}

void VariationalBayesEstimatorOnLDA::updateHyperParameters(){//{{{
    double denominator = 0;
    for(int d=0;d<_D;d++){
        denominator += boost::math::digamma(_nd[d]+_alphaSum) - boost::math::digamma(_alphaSum);
    }
    for(int k=0;k<_K;k++){
        double numerator = 0;
        for(int d=0;d<_D;d++){
            //ndk == 0 -> alphak == 0  digamma(0) - digamma(0) = 0
            try{
                numerator += (boost::math::digamma(_ndk[d][k]+_alpha[k]) - boost::math::digamma(_alpha[k])) * _alpha[k];
            }catch(...){
                numerator += 0;
            }
        }
        _alpha[k] = numerator / denominator;
    }
    _alphaTimeSeries.push_back(_alpha);
    this->updateBeta(ASYMMETRY);
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
    for(int d=0; d<_D; d++){
        term2 += boost::math::lgamma(_alphaSum) - boost::math::lgamma(_nd[d]+_alphaSum);
        for(int k=0; k<_K; k++){
        //ndk == 0 -> alphak == 0  lgamma(0) - lgamma(0) = 0
            try{
                term2 += boost::math::lgamma(_ndk[d][k]+_alpha[k]) - boost::math::lgamma(_alpha[k]);
            }catch(...){
                term2 += 0;
            }
        }
    }
    double term3 = 0;
    int check = 0;
    for(int d=0;d<_docVoca.size();d++){
        for(int l=0; l<_docVoca[d].size(); l++){
            unsigned int v = _docVoca[d][l];
            for(int k=0; k<_K; k++){
                if(_qz[d][l][k] != 0){
                    term3 += _qz[d][l][k] * log(_qz[d][l][k]) * _frequencyMatrix[d][v];
                }
            }
        }
    }
    double variationalLowerBound = term1 + term2 - term3;
    return(variationalLowerBound);
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
        this->updateQz();
        this->updateNEx();
        thisVariationalLowerBound = this->calculateVariationalLowerBound();
        cout<<"VLB"<<thisVariationalLowerBound<<endl;
        _VLBTimeSeries.push_back(thisVariationalLowerBound);
        count++;
        if(std::isnan(thisVariationalLowerBound))break;
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
