//
// VariationalBayesEstimatorOnATM.cpp
//
// Copyright (c) 2017 Shion Hosoda
//
// This software is released u_nder the MIT License.
// http://opensource.org/licenses/mit-license.php
//

#include"include/VariationalBayesEstimatorOnATM.hpp"

using namespace std;

VariationalBayesEstimatorOnATM::VariationalBayesEstimatorOnATM(const BOWFileParser &parser, const AuthorFileParser &aParser, const unsigned int K, const double convergenceDiterminationRate)//{{{
    :_frequencyMatrix(parser.getFrequencyMatrix()),
     _docVoca(parser.getDocVoca()),
     _K(K),
     _D(_frequencyMatrix.size()),
     _V(parser.getV()),
     _M(aParser.getM()),
     _docAuth(aParser.getDocAuth()),
     _convergenceDiterminationRate(convergenceDiterminationRate),
     _qzy(_D),
     _thetaEx(_M),
     _phiEx(K),
     _alphaSum(0),
     _betaSum(0),
     _nk(K),
     _nm(_M),
     _nkv(K),
     _nmk(_M){
    this->initializeParam();
    this->initializeHyperParam();
}//}}}

VariationalBayesEstimatorOnATM::~VariationalBayesEstimatorOnATM(){//{{{

}//}}}

void VariationalBayesEstimatorOnATM::initializeParam(){//{{{
    for(int m=0; m<_thetaEx.size(); m++)_thetaEx[m].assign(_K, 0);
    for(int k=0; k<_phiEx.size(); k++)_phiEx[k].assign(_V, 0);
    for(int m=0; m<_nmk.size(); m++)_nmk[m].assign(_K, 0);
    for(int k=0; k<_nkv.size(); k++)_nkv[k].assign(_V, 0);
    _nk.assign(_K, 0);
    _nm.assign(_M, 0);

    random_device rnd;
    mt19937 mt(rnd());
    uniform_real_distribution<double> randN(0, 1);
    for(int d=0; d<_docVoca.size(); d++){
        _qzy[d] = vector<vector<vector<double> > >(_docVoca[d].size());
        for(int l=0; l<_docVoca[d].size(); l++){
            _qzy[d][l] = vector<vector<double> >(_K);
            double qzySum = 0.0;
            for(int k=0; k<_K; k++){
                _qzy[d][l][k] = vector<double>(_docAuth.size());
                for(int r=0; r<_docAuth[d].size(); r++){
                    double randomValue = randN(mt);
                    _qzy[d][l][k][r] = randomValue;
                }
                qzySum = accumulate(_qzy[d][l][k].begin(), _qzy[d][l][k].end(), 0.0);
            }
            for(int k=0; k<_K; k++){
                for(int r=0; r<_docAuth[d].size(); r++){
                    _qzy[d][l][k][r] /= qzySum;
                }
            }
        }
    }

}//}}}

void VariationalBayesEstimatorOnATM::initializeHyperParam(){//{{{
    _alpha.assign(_K, 0.1);
    _beta.assign(_V, 0.05);
    _alphaTimeSeries.push_back(_alpha);
    _betaTimeSeries.push_back(_beta);
    this->calculateHyperParamSum();
}//}}}

void VariationalBayesEstimatorOnATM::calculateEx(){//{{{
    for(int m=0;m<_M;m++){
        for(int k=0;k<_K;k++){
            _thetaEx[m][k]=(_nmk[m][k]+_alpha[k])/(_nm[m]+_alphaSum);
        }
    }
    for(int k=0;k<_K;k++){
        for(int v=0;v<_V;v++){
            _phiEx[k][v]=(_nkv[k][v]+_beta[v])/(_nk[k]+_betaSum);
        }
    }

}//}}}

void VariationalBayesEstimatorOnATM::calculateHyperParamSum(){//{{{
    _alphaSum = accumulate(_alpha.begin(), _alpha.end(), 0.0);
    _betaSum = accumulate(_beta.begin(), _beta.end(), 0.0);
}//}}}

void VariationalBayesEstimatorOnATM::updateQzy(){//{{{
    for(int d=0;d<_docVoca.size();d++){
        double pidm = 1.0 / _docAuth[d].size();
        for(int l=0; l<_docVoca[d].size(); l++){
            double Zq = 0;
            unsigned int v = _docVoca[d][l];
            for(int k=0; k<_K; k++){
                for(int r=0; r<_docAuth[d].size(); r++){
                    unsigned int m = _docAuth[d][r];
                    double xikv = _nkv[k][v] + _beta[v];
                    double ximk = _nmk[m][k] + _alpha[k];
                    double xik = _nk[k] + _betaSum;
                    double xim = _nm[m] + _alphaSum;
                    double numerator = 0, denominator = 0;
                    numerator = exp(boost::math::digamma(xikv)) * exp(boost::math::digamma(ximk));
                    denominator = exp(boost::math::digamma(xik)) * exp(boost::math::digamma(xim));
                    double constProbability = pidm * numerator / denominator;
                    Zq += constProbability;
                    _qzy[d][l][k][r] = constProbability;
                }
            }
            for(int k=0; k<_K; k++){
                for(int r=0; r<_docAuth[d].size(); r++){
                    _qzy[d][l][k][r] /= Zq;
                }
            }
        }
    }
}//}}}

void VariationalBayesEstimatorOnATM::updateNEx(){//{{{
    vector<vector<double> > nkvBuf(_nkv.size()), nmkBuf(_nmk.size());
    for(int k=0;k<_nkv.size();k++)nkvBuf[k].assign(_V, 0);
    for(int m=0;m<_nmk.size();m++)nmkBuf[m].assign(_K, 0);
    for(int d=0;d<_D;d++){
        for(int l=0; l<_docVoca[d].size(); l++){
            unsigned int v = _docVoca[d][l];
            for(int k=0; k<_K; k++){
                for(int r=0; r<_docAuth[d].size(); r++){
                    unsigned int m = _docAuth[d][r];
                    nmkBuf[m][k] += _qzy[d][l][k][r] * _frequencyMatrix[d][v];
                    nkvBuf[k][v] += _qzy[d][l][k][r] * _frequencyMatrix[d][v];
                }
            }
        }
    }
    for(int m=0; m<nmkBuf.size(); m++){
        _nm[m] = 0;
        for(int k=0;k<_K;k++){
            _nmk[m][k] = nmkBuf[m][k];
            _nm[m] += nmkBuf[m][k];
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

void VariationalBayesEstimatorOnATM::updateBeta(BetaUpdateManner manner){//{{{
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

void VariationalBayesEstimatorOnATM::updateHyperParameters(){//{{{
    double denominator = 0;
    for(int m=0; m<_nm.size(); m++){
        denominator += boost::math::digamma(_nm[m]+_alphaSum) - boost::math::digamma(_alphaSum);
    }
    for(int k=0;k<_K;k++){
        double numerator = 0;
        for(int m=0; m<_nm.size(); m++){
            numerator += (boost::math::digamma(_nmk[m][k]+_alpha[k]) - boost::math::digamma(_alpha[k])) * _alpha[k];
        }
        _alpha[k] = numerator / denominator;
    }
    _alphaTimeSeries.push_back(_alpha);
    this->updateBeta(ASYMMETRY);
    this->calculateHyperParamSum();
}//}}}

double VariationalBayesEstimatorOnATM::calculateVariationalLowerBound()const{//{{{
    double term1=0;
    for(int k=0; k<_K; k++){
        term1 += boost::math::lgamma(_betaSum) - boost::math::lgamma(_nk[k]+_betaSum);
        for(int v=0; v<_V; v++){
            term1 += boost::math::lgamma(_nkv[k][v]+_beta[v]) - boost::math::lgamma(_beta[v]);
        }
    }
    double term2=0;
    for(int m=0; m<_M; m++){
        term2 += boost::math::lgamma(_alphaSum) - boost::math::lgamma(_nm[m]+_alphaSum);
        for(int k=0; k<_K; k++){
            term2 += boost::math::lgamma(_nmk[m][k]+_alpha[k]) - boost::math::lgamma(_alpha[k]);
        }
    }
    double term3 = 0;
    for(int d=0;d<_docVoca.size();d++){
        for(int l=0; l<_docVoca[d].size(); l++){
            unsigned int v = _docVoca[d][l];
            for(int k=0; k<_K; k++){
                for(int r=0; r<_docAuth[d].size(); r++){
                    if(_qzy[d][l][k][r] != 0){
                        term3 += _qzy[d][l][k][r] * log(_qzy[d][l][k][r]) * _frequencyMatrix[d][v];
                    }
                }
            }
        }
    }
    double term4 = 0;
    for(int d=0;d<_docVoca.size();d++){
        for(int l=0; l<_docVoca[d].size(); l++){
            unsigned int v = _docVoca[d][l];
            for(int r=0; r<_docAuth[d].size(); r++){
                double qydim = 0.0;
                for(int k=0; k<_K; k++){
                    qydim += _qzy[d][l][k][r];
                }
                term4 += qydim * log(1.0 / _docAuth[d].size()) * _frequencyMatrix[d][v];
            }
        }
    }

    double variationalLowerBound = term1 + term2 - term3 + term4;
    return(variationalLowerBound);
}//}}}

void VariationalBayesEstimatorOnATM::writeParameter(string thetaFilename, string phiFilename, string alphaFilename, string betaFilename)const{//{{{
    outputVector(_thetaEx, thetaFilename);
    outputVector(_phiEx, phiFilename);
    outputVector(_alphaTimeSeries, alphaFilename);
    outputVector(_betaTimeSeries, betaFilename);
}//}}}

void VariationalBayesEstimatorOnATM::writeVariationalLowerBound(string VLBFilename, string VLBTimeSeriesFilename)const{//{{{
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

void VariationalBayesEstimatorOnATM::runIteraions(){//{{{
    double prevVariationalLowerBound = 0;
    double thisVariationalLowerBound = 0;
    unsigned int count = 0;
    while(1){
        prevVariationalLowerBound = thisVariationalLowerBound;
        this->updateQzy();
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
