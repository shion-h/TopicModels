//
// GibbsSamplerFromHDP.cpp
//
// Copyright (c) 2018 Shion Hosoda
//
// This software is released u_nder the MIT License.
// http://opensource.org/licenses/mit-license.php
//

#include"include/GibbsSamplerFromHDP.hpp"

using namespace std;

GibbsSamplerFromHDP::GibbsSamplerFromHDP(const BOWFileParser &parser, double alpha, double beta, double gamma, unsigned int iterationNumber, unsigned int burnIn, unsigned int samplingInterval)//{{{
    :_frequencyMatrix(parser.getFrequencyMatrix()),
     _docVoca(parser.getDocVoca()),
     _iterationNumber(iterationNumber),
     _burnIn(burnIn),
     _samplingInterval(samplingInterval),
     _samplingCount(0),
     _alpha(alpha),
     _beta(beta),
     _gamma(gamma),
     _J(_frequencyMatrix.size()),
     _V(parser.getV()),
     _K(1),
     _Tj(_J, 1),
     _nj(_J),
     _njtv(_J),
     _nk(_K),
     _mk(_K)
{
    this->initializeParameters();
    this->countN();
}//}}}

GibbsSamplerFromHDP::~GibbsSamplerFromHDP(){//{{{

}//}}}

void GibbsSamplerFromHDP::initializeParameters(){//{{{
    _Tj.assign(_J, 1);
    for(int j=0; j<_J; j++){
        _nj[j] = accumulate(_frequencyMatrix[j].begin(), _frequencyMatrix[j].end(), 0u);
    }
    // initialize xji with frequencyMatrix
    for(int j=0; j<_J; j++){
        vector<unsigned int> bufVector;
        for(int v=0; v<_V; v++){
            for(int i=0; i<_frequencyMatrix[j][v]; i++){
                bufVector.push_back(v);
            }
        }
        _xji.push_back(bufVector);
    }
    for(int j=0; j<_J; j++){
        vector<unsigned int> bufVector;
        for(int i=0; i<_nj[j]; i++){
            bufVector.push_back(0);
        }
        _tji.push_back(bufVector);
    }
    unsigned int initialTNumber = 1;
    for(int j=0; j<_J; j++){
        vector<unsigned int> bufVector(initialTNumber, 0);
        _kjt.push_back(bufVector);
        _njt.push_back(bufVector);
    }
    _nkv = vector<vector<unsigned int> >(_K, vector<unsigned int>(_V, 0));
    _nkv = vector<vector<unsigned int> >(_K, vector<unsigned int>(_V, 0));
    _njtv = vector<vector<vector<unsigned int> > >(_J, vector<vector<unsigned int> >(initialTNumber, vector<unsigned int>(_V, 0)));
    unsigned int nSum = 0;
    for(int j=0; j<_J; j++){
        nSum += _nj[j];
    }
    _nk = vector<unsigned int>(_K, nSum/_K);
    _mk = vector<unsigned int>(_K, _J * 1 / _K);
    _sumOfNkSampled = vector<double>(_K, 0);
    _sumOfNjkSampled = vector<vector<double> >(_J, vector<double>(_K, 0));
    _sumOfNkvSampled = vector<vector<double> >(_K, vector<double>(_V, 0));
}//}}}

void GibbsSamplerFromHDP::sampleTji(){//{{{
    unsigned int mkSum = accumulate(_mk.begin(), _mk.end(), 0u);
    for(int j=0; j<_J; j++){
        for(int i=0; i<_nj[j]; i++){
            unsigned int v = _xji[j][i];
            unsigned int tOld = _tji[j][i];
            unsigned int kOld = _kjt[j][tOld];
            _njt[j][tOld]--;
            _njtv[j][tOld][v]--;
            _nkv[kOld][v]--;
            _nk[kOld]--;
            if(_njt[j][tOld] == 0) this->deleteT(j, tOld);
            // calculete pOfXjiGivenTnew and substitute pOfXjiGivenKused
            vector<double> pOfXjiGivenKused(_K);
            double pOfXjiGivenTnew = 0;
            for(int k=0; k<_K; k++){
                pOfXjiGivenKused[k] = (_beta + _nkv[k][v]) / (_V * _beta + _nk[k]);
                pOfXjiGivenTnew += _mk[k] * pOfXjiGivenKused[k];
            }
            pOfXjiGivenTnew += _gamma / _V;
            pOfXjiGivenTnew /= mkSum + _gamma;
            // t previously used
            vector<double> distributionOfT;
            for(int t=0; t<_Tj[j]; t++){
                unsigned int k = _kjt[j][t];
                distributionOfT.push_back(_njt[j][t] * pOfXjiGivenKused[k]);
            }
            // t new
            distributionOfT.push_back(_alpha * pOfXjiGivenTnew);
            _tji[j][i] = sampleDiscreteValues(distributionOfT);
            if(_tji[j][i] > _Tj[j]) exit(0);
            bool isNewT = false;
            if(_tji[j][i] == _Tj[j]){
                isNewT = true;
                this->addT(j, i);
            }
            //if t is new, already added so far
            if(not isNewT){
                _njt[j][_tji[j][i]]++;
                _njtv[j][_tji[j][i]][v]++;
            }
            unsigned int kNew = _kjt[j][_tji[j][i]];
            _nkv[kNew][v]++;
            _nk[kNew]++;
        }
    }
}//}}}

void GibbsSamplerFromHDP::sampleKjt(unsigned int j, unsigned int t, bool isNewT/*=false*/){//{{{
    if(not isNewT){
        unsigned int kOld = _kjt[j][t];
        for(int v=0; v<_V; v++){
            _nkv[kOld][v] -= _njtv[j][t][v];
        }
        _nk[kOld] -= _njt[j][t];
        _mk[kOld]--;
        if(_mk[kOld] == 0){
            if(_nk[kOld] != 0){
                cout<<"nk error"<<endl;
                cout<<_mk[kOld]<<' '<<_nk[kOld]<<endl;
                exit(0);
            }
            this->deleteK(kOld);
        }
    }
    // k used
    vector<double> logDistributionOfK;
    for(int k=0; k<_K; k++){
        double term1 = 0;
        for(int v=0; v<_V; v++){
            term1 += boost::math::lgamma(_beta + _nkv[k][v] + _njtv[j][t][v]);
        }
        double term2 = boost::math::lgamma(_V * _beta + _nk[k] + _njt[j][t]);
        double term3 = 0;
        for(int v=0; v<_V; v++){
            term3 += boost::math::lgamma(_beta + _nkv[k][v]);
        }
        double term4 = boost::math::lgamma(_V * _beta + _nk[k]);
        double logPOfXjtGivenKused = term1 - term2 - term3 + term4;

        logDistributionOfK.push_back(log(_mk[k]) +  logPOfXjtGivenKused);
    }
    // k new
    double term1 = boost::math::lgamma(_V * _beta);
    double term2 = 0;
    for(int v=0; v<_V; v++){
        term2 += boost::math::lgamma(_beta + _njtv[j][t][v]);
    }
    double term3 = boost::math::lgamma(_V * _beta + _njt[j][t]);
    double term4 = _V * boost::math::lgamma(_beta);
    double logPOfXjtGivenKnew = term1 + term2 - term3 - term4;

    logDistributionOfK.push_back(log(_gamma) + logPOfXjtGivenKnew);
    _kjt[j][t] = sampleDiscreteValues(logDistributionOfK, true);
    if(_kjt[j][t] > _K){
        cout<<"K discrete distribution error"<<endl;
        cout<<_kjt[j][t]<<' '<<_K<<endl;
        printVector(logDistributionOfK);
        exit(0);
    }
    if(_kjt[j][t] == _K) this->addK();
    if(not isNewT){
        for(int v=0; v<_V; v++){
            _nkv[_kjt[j][t]][v] += _njtv[j][t][v];
        }
        _nk[_kjt[j][t]] += _njt[j][t];
    }
    _mk[_kjt[j][t]]++;
}//}}}

void GibbsSamplerFromHDP::countN(){//{{{
    // count njt from tji
    fillMatrix(_njt, 0u);
    fillTensor(_njtv, 0u);
    for(int j=0; j<_J; j++){
        for(int i=0; i<_nj[j]; i++){
            unsigned int t = _tji[j][i];
            unsigned int v = _xji[j][i];
            _njtv[j][t][v]++;
        }
    }
    for(int j=0; j<_J; j++){
        for(int t=0; t<_Tj[j]; t++){
            _njt[j][t] = accumulate(_njtv[j][t].begin(), _njtv[j][t].end(), 0u);
        }
    }
    // delete t & update tji
    // count njt from tji
    fill(_nk.begin(), _nk.end(), 0);
    // count nkv from kjtji & xji
    fillMatrix(_nkv, 0u);
    for(int j=0; j<_J; j++){
        for(int i=0; i<_nj[j]; i++){
            unsigned int t = _tji[j][i];
            unsigned int k = _kjt[j][t];
            unsigned int v = _xji[j][i];
            _nkv[k][v]++;
        }
    }
    for(int k=0; k<_K; k++){
        _nk[k] = accumulate(_nkv[k].begin(), _nkv[k].end(), 0u);
    }
    fill(_mk.begin(), _mk.end(), 0);
    // mk
    for(int j=0; j<_J; j++){
        for(int t=0; t<_Tj[j]; t++){
            _mk[_kjt[j][t]]++;
        }
    }
}//}}}

void GibbsSamplerFromHDP::addT(unsigned int j, unsigned int i){//{{{
    _Tj[j]++;
    _kjt[j].push_back(-1);
    _njt[j].push_back(1);
    _njtv[j].push_back(vector<unsigned int>(_V, 0));
    _njtv[j][_Tj[j]-1][_xji[j][i]]++;
    this->sampleKjt(j, _Tj[j]-1, true);
}//}}}

void GibbsSamplerFromHDP::addK(){//{{{
    _K++;
    _nk.push_back(0);
    _mk.push_back(0);
    _nkv.push_back(vector<unsigned int>(_V, 0));
}//}}}

void GibbsSamplerFromHDP::deleteT(unsigned int j, unsigned int t){//{{{ // revised tji
    _mk[_kjt[j][t]]--;
    if(_mk[_kjt[j][t]] == 0)this->deleteK(_kjt[j][t]);
    for(int i=0; i<_nj[j]; i++){
        if(_tji[j][i] > t) _tji[j][i]--;
    }
    _kjt[j].erase(_kjt[j].begin() + t);
    _njt[j].erase(_njt[j].begin() + t);
    _njtv[j].erase(_njtv[j].begin() + t);
    _Tj[j]--;
}//}}}

void GibbsSamplerFromHDP::deleteK(unsigned int k){//{{{
    for(int j=0; j<_J; j++){
        for(int t=0; t<_Tj[j]; t++){
            if(_kjt[j][t] > k) _kjt[j][t]--;
        }
    }
    _nkv.erase(_nkv.begin() + k);
    _nk.erase(_nk.begin() + k);
    _mk.erase(_mk.begin() + k);
    _K--;
}//}}}

void GibbsSamplerFromHDP::calculateExpectations(bool isUseOnlyLastSample/*=false*/){//{{{
    vector<double> nk(_K, 0);
    for(int k=0; k<_K; k++){
        // nk[k] = _sumOfNkSampled[k] / _samplingCount;
        if(isUseOnlyLastSample) nk[k] = _nk[k];
    }
    double nkSum = accumulate(nk.begin(), nk.end(), 0.0);
    _thetaEx = vector<vector<double> >(_J, vector<double>(_K, 0));
    _phiEx = vector<vector<double> >(_K, vector<double>(_V, 0));
    for(int j=0; j<_J; j++){
        for(int k=0; k<_K; k++){
            // double njk = _sumOfNjkSampled[j][k] / _samplingCount;
            double njk = 0;
            if(isUseOnlyLastSample){
                for(int t=0; t<_Tj[j]; t++){
                    if(k == _kjt[j][t]) njk += _njt[j][t];
                }
            }
            _thetaEx[j][k] = ((nkSum + _gamma) * njk + _alpha * nk[k]) / ((nkSum + _gamma) * _nj[j] + _alpha * nkSum);
        }
    }
    for(int k=0; k<_K; k++){
        for(int v=0; v<_V; v++){
            // double nkv = _sumOfNkvSampled[k][v] / _samplingCount;
            double nkv;
            if(isUseOnlyLastSample) nkv = _nkv[k][v];
            _phiEx[k][v] = (nkv + _beta) / (nk[k] + _V * _beta);
        }
    }
}//}}}

void GibbsSamplerFromHDP::calculateLogLikelihood(){//{{{
    this->calculateExpectations(true);
    double logLikelihood = 0;
    vector<vector<double> > probabilityMatrix = dot(_thetaEx, _phiEx);
    for(int j=0; j<_J; j++){
        for(int i=0; i<_nj[j]; i++){
            for(int k=0; k<_K; k++){
                logLikelihood += log(probabilityMatrix[j][_xji[j][i]]);
            }
        }
    }
    cout<<logLikelihood<<endl;
    _logLikelihood.push_back(logLikelihood);
}//}}}

void GibbsSamplerFromHDP::checkNumbers(){//{{{
    if(sum(_nkv) != sum(_nk) || sum(_nk) != sum(_nj) || sum(_nj) != sum(_njt) || sum(_njt) != sum(_njtv))exit(0);
    if(_K != _nkv.size() || _nkv.size() != _mk.size() || _mk.size() != _nk.size())exit(0);
}//}}}

void GibbsSamplerFromHDP::writeParameters(string thetaFilename, string phiFilename)const{//{{{
    outputVector(_thetaEx, thetaFilename);
    outputVector(_phiEx, phiFilename);
}//}}}

void GibbsSamplerFromHDP::writeLogLikelihood(string logLikelihoodFilename)const{//{{{
    outputVector(_logLikelihood, logLikelihoodFilename);
}//}}}

void GibbsSamplerFromHDP::storeSamples(){//{{{
    _samplingCount++;
    for(int k=0; k<_K; k++){
        _sumOfNkSampled[k] += _nk[k];
    }
    //nkv
    for(int k=0; k<_K; k++){
        for(int v=0; v<_V; v++){
            _sumOfNkvSampled[k][v] += _nkv[k][v];
        }
    }
    //njk
    for(int j=0; j<_J; j++){
        for(int t=0; t<_Tj[j]; t++){
            _sumOfNjkSampled[j][_kjt[j][t]] += _njt[j][t];
        }
    }
}//}}}

void GibbsSamplerFromHDP::runIteraions(){//{{{
    for(int i=0; i<_iterationNumber; i++){
        this->sampleTji();
        for(int j=0; j<_J; j++){
            for(int t=0; t<_Tj[j]; t++){
                this->sampleKjt(j, t);
            }
        }
        this->calculateExpectations(true);
        this->calculateLogLikelihood();
        // this->checkNumbers();
        cout<<_K<<endl;
        // if((i > _burnIn - 1) && ((i + 1) % _samplingInterval == 0)){
            // this->storeSamples();
        // }
    }
    // this->calculateExpectations();
}//}}}
