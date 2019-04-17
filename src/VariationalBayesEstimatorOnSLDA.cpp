//
// VariationalBayesEstimatorOnSLDA.cpp
//
// Copyright (c) 2017 Shion Hosoda
//
// This software is released u_nder the MIT License.
// http://opensource.org/licenses/mit-license.php
//

#include"include/VariationalBayesEstimatorOnSLDA.hpp"

using namespace std;

VariationalBayesEstimatorOnSLDA::VariationalBayesEstimatorOnSLDA(const BOWFileParser &parser, const LabelFileParser &labelFileParser, const unsigned int K, const double convergenceDiterminationRate)//{{{
    :VariationalBayesEstimatorOnLDA(parser, K, convergenceDiterminationRate),
     _y(labelFileParser.getY()),
     _C(labelFileParser.getC()),
     _xi(_D, 1),
     _eta(_C),
     _expEtaZ(_D){
    this->initializeParam();
    this->initializeHyperParam();
}//}}}

VariationalBayesEstimatorOnSLDA::~VariationalBayesEstimatorOnSLDA(){//{{{

}//}}}

void VariationalBayesEstimatorOnSLDA::initializeParam(){//{{{
    VariationalBayesEstimatorOnLDA::initializeParam();
    random_device rnd;
    mt19937 mt(rnd());
    uniform_real_distribution<double> rand0(-1, 1);
    for(int c=0; c<_eta.size();c++) _eta[c].assign(_K, 0);
    for(int c=0; c<_eta.size(); c++){
        for(int k=0; k<_eta[c].size(); k++){
            double randomValue = rand0(mt);
            _eta[c][k] = randomValue;
        }
    }
    for(int d=0; d<_expEtaZ.size();d++) _expEtaZ[d].assign(_C, 0);
}//}}}

void VariationalBayesEstimatorOnSLDA::updateExpEtaZ(int c){//{{{
    for(int d=0; d<_qz.size(); d++){
        _expEtaZ[d][c] = 1;
        for(int l=0; l<_qz[d].size(); l++){
            double expEtacZdi = 0;
            for(int k=0; k<_K; k++){
                expEtacZdi += _qz[d][l][k] * exp(_eta[c][k]/_nd[d]);
            }
            unsigned int v = _docVoca[d][l];
            _expEtaZ[d][c] *= pow(expEtacZdi, _frequencyMatrix[d][v]);
        }
    }
}//}}}

void VariationalBayesEstimatorOnSLDA::updateExpEtaZ(){//{{{
    for(int c=0; c<_C; c++){
        this->updateExpEtaZ(c);
    }
}//}}}

void VariationalBayesEstimatorOnSLDA::updateXi(){//{{{
    this->updateExpEtaZ();
    for(int d=0; d<_qz.size(); d++){
        _xi[d] = 0;
        for(int c=0; c<_C; c++){
            _xi[d] += _expEtaZ[d][c];
        }
    }
}//}}}

void VariationalBayesEstimatorOnSLDA::updateQz(){//{{{
    this->updateExpEtaZ();
    for(int d=0; d<_docVoca.size(); d++){
        for(int l=0; l<_docVoca[d].size(); l++){
            vector<double> expEtaZdi(_C, 0);
            for(int c=0; c<_C; c++){
                for(int k=0; k<_K; k++){
                    expEtaZdi[c] += _qz[d][l][k] * exp(_eta[c][k] / _nd[d]);
                }
            }
            double Zq = 0;
            unsigned int v = _docVoca[d][l];
            for(int k=0; k<_K; k++){
                double term1 = 0, term2 = 0, term3 = 0, term4 = 0;
                double xikv = _nkv[k][v]+_beta[v];
                double xidk = _ndk[d][k]+_alpha[k];
                double xik = _nk[k]+_betaSum;
                double xid = _nd[d]+_alphaSum;
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
                double dxidqzdik = 0;
                for(int c=0; c<_C; c++){
                    dxidqzdik -= exp(_eta[c][k] / _nd[d]) * (_expEtaZ[d][c] / expEtaZdi[c]);
                }
                double logisticTerm = _eta[_y[d]][k]/_nd[d] - dxidqzdik/_xi[d];
                double constProbability = (term1 * term2) / (term3 * term4) * logisticTerm;
                Zq += constProbability;
                _qz[d][l][k] = constProbability;
            }
            for(int k=0; k<_K; k++){
                _qz[d][l][k] /= Zq;
            }
        }
    }
}//}}}

vector<double> VariationalBayesEstimatorOnSLDA::calculateNablaL(unsigned int c)const{//{{{
    vector<double> nablaLc(_K, 0);
    for(int d=0; d<_qz.size(); d++){
        // yd == c
        if(_y[d] == c){
            vector<double> zd(_K);
            // q(zd)[k]
            for(int l=0; l<_qz[d].size(); l++){
                unsigned int v = _docVoca[d][l];
                for(int k=0; k<_K; k++){
                    zd[k] += _qz[d][l][k] * _frequencyMatrix[d][v];
                }
            }
            // q(zd_bar)[k]
            for(int k=0; k<_K; k++){
                nablaLc[k] += zd[k] / _nd[d];
            }
            //d E[exp(eta[c] . zd_bar)] / d eta
            vector<double> dExpEtacZddEta(_K);
            for(int l=0; l<_qz[d].size(); l++){
                double expEtacZdi = 0;
                // exp(eta[c] . q(zd_bar))
                for(int k=0; k<_K; k++){
                    expEtacZdi += _qz[d][l][k] * exp(_eta[c][k] / _nd[d]);
                }
                // q(z[d,i] = k)/nd * exp(eta[c,k]/nd) * E\i [exp(eta[c] . zd_bar)]
                unsigned int v = _docVoca[d][l];
                for(int k=0; k<_K; k++){
                    dExpEtacZddEta[k] += (_qz[d][l][k] / _nd[d]) * exp(_eta[c][k] / _nd[d]) * _expEtaZ[d][c] / expEtacZdi * _frequencyMatrix[d][v];
                }
            }
            for(int k=0; k<_K; k++){
                nablaLc[k] -= dExpEtacZddEta[k] / _xi[d];
            }
        }
    }
    return nablaLc;
}//}}}

void VariationalBayesEstimatorOnSLDA::updateEta(){//{{{

    class LambdaUpdater{//{{{
    private:
        unsigned int _c;
        vector<double> _d;
        VariationalBayesEstimatorOnSLDA *_estimator;
    public:
        LambdaUpdater(unsigned int c, vector<double> d, VariationalBayesEstimatorOnSLDA *estimator)//{{{
        :_c(c),
         _d(d),
         _estimator(estimator){}//}}}

        double calculateLc(double lambda){//{{{
            // calculate the terms in L having etac
            double Lc = 0;
            vector<double> etac(_estimator->_K);
            for(int k=0; k<etac.size(); k++){
                etac[k] = _estimator->_eta[_c][k] + lambda * _d[k];
            }
/* if(_c == 0){
outputVectorForDebug(etac ,"etac");
// for(int k=0; k<etac.size(); k++) cout<<etac[k]<<',';
// cout<<endl;
} */
            for(int d=0; d<_estimator->_qz.size(); d++){
                double etaydZd = 0;
                if(_estimator->_y[d] == _c){
                    for(int l=0; l<_estimator->_qz[d].size(); l++){
                        double etaydZdi = 0;
                        for(int k=0; k<_estimator->_K; k++){
                            etaydZdi += etac[k] * _estimator->_qz[d][l][k];
                        }
                        unsigned int v = _estimator->_docVoca[d][l];
                        etaydZd += etaydZdi * _estimator->_frequencyMatrix[d][v];
                    }
                    etaydZd /= _estimator->_nd[d];
                }
                double expEtacZd = 1.0;
                for(int l=0; l<_estimator->_qz[d].size(); l++){
                    double expEtacZdi = 0;
                    for(int k=0; k<_estimator->_K; k++){
                        expEtacZdi += _estimator->_qz[d][l][k] * exp(etac[k]/_estimator->_nd[d]);
                    }
                    unsigned int v = _estimator->_docVoca[d][l];
                    expEtacZd *= pow(expEtacZdi, _estimator->_frequencyMatrix[d][v]);
                }
                Lc += etaydZd - expEtacZd/_estimator->_xi[d];
            }
            return Lc;
        }//}}}

        double calculateLambdaArgmaxLc(){//{{{
            //maxとるようになってるか？
            double min = 0.0, max = 100, goldenRatio = (1 + sqrt(5)) / 2;
            double lambda1, lambda2;
            double L1, L2;
            lambda1 = (max - min)/(goldenRatio + 1.0) + min;
            L1 = this->calculateLc(lambda1);
            lambda2 = (max - min)/goldenRatio + min;
            L2 = this->calculateLc(lambda2);
            while((max - min) > 0.00001){
// cout<<lambda1<<' '<<lambda2<<endl;
// cout<<L1<<' '<<L2<<endl;
                if(L1 < L2){
                    min = lambda1;
                    lambda1 = lambda2;
                    L1 = L2;
                    lambda2 = (max - min)/goldenRatio + min;
                    L2 = this->calculateLc(lambda2);
                }else{
                    max = lambda2;
                    lambda2 = lambda1;
                    L2 = L1;
                    lambda1 = (max - min)/(goldenRatio + 1.0) + min;
                    L1 = this->calculateLc(lambda1);
                }
            }
            return (min + max)/2.0;
        }//}}}
    };//}}}

    for(int c=0; c<_C; c++){
// cout<<"c: "<<c<<endl;
        this->updateExpEtaZ(c);
        vector<double> thisNablaLc = this->calculateNablaL(c);
        vector<double> prevNablaLc;
        double zeta = 0, lambda = 0;
        vector<double> d(thisNablaLc);
        while(1){
            LambdaUpdater *updater;
            updater = new LambdaUpdater(c, d, this);
            lambda = updater->calculateLambdaArgmaxLc();
/*            cout<<"lambda zeta ";
cout<<lambda<<' '<<zeta<<' ';
            cout<<endl;
            cout<<"eta ";
            for(int k=0; k<_K; k++) cout<<_eta[c][k]<<' ';
            cout<<endl;
            cout<<"d ";
            for(int k=0; k<_K; k++) cout<<d[k]<<' ';
            cout<<endl;
            cout<<endl;*/
            if(lambda < 0.00001)break;
            for(int k=0; k<_K; k++) _eta[c][k] = _eta[c][k] + lambda * d[k];
            prevNablaLc = thisNablaLc;
            this->updateExpEtaZ(c);
            thisNablaLc = this->calculateNablaL(c);
            double denominator, numerator = 0;
            for(int k=0; k<_K; k++){
                numerator += thisNablaLc[k] * (thisNablaLc[k] - prevNablaLc[k]);
                denominator += prevNablaLc[k] * prevNablaLc[k];
            }
            zeta = numerator / denominator;
            for(int k=0; k<_K; k++) d[k] = thisNablaLc[k] + zeta * d[k];
            delete(updater);
        }
    }
}//}}}

double VariationalBayesEstimatorOnSLDA::calculateVariationalLowerBound()const{//{{{
    double LDAVariationalLowerBound = VariationalBayesEstimatorOnLDA::calculateVariationalLowerBound();
    double logisticTerm = 0;
    for(int d=0; d<_qz.size(); d++){
        double etaydZd = 0;
        for(int l=0; l<_qz[d].size(); l++){
            double etaydZdi = 0;
            for(int k=0; k<_K; k++){
                etaydZdi += _eta[_y[d]][k] * _qz[d][l][k];
            }
            unsigned int v = _docVoca[d][l];
            etaydZd += etaydZdi * _frequencyMatrix[d][v];
        }
        etaydZd /= _nd[d];
        logisticTerm += etaydZd - log(_xi[d]);
    }
    // _variationalLowerBoundOfQzがこちら側の使われているかチェック
    double variationalLowerBound = LDAVariationalLowerBound + logisticTerm;
    return(variationalLowerBound);
}//}}}

void VariationalBayesEstimatorOnSLDA::writeParameter(string thetaFilename, string phiFilename, string etaFilename, string xiFilename, string expEtaZFilename, string alphaFilename, string betaFilename)const{//{{{
    VariationalBayesEstimatorOnLDA::writeParameter(thetaFilename, phiFilename, alphaFilename, betaFilename);
    outputVector(_eta, etaFilename);
    outputVector(_xi, xiFilename);
    outputVector(_expEtaZ, expEtaZFilename);
}//}}}

void VariationalBayesEstimatorOnSLDA::runIteraions(){//{{{
    double prevVariationalLowerBound = 0;
    double thisVariationalLowerBound = 0;
    unsigned int count = 0;
    while(1){
        prevVariationalLowerBound = thisVariationalLowerBound;
        this->updateQz();
        this->updateNEx();
        this->updateXi();
        thisVariationalLowerBound = this->calculateVariationalLowerBound();
        // cout<<"VLB"<<thisVariationalLowerBound<<endl;
        _VLBTimeSeries.push_back(thisVariationalLowerBound);
        count++;
        if(std::isnan(thisVariationalLowerBound))break;
        if(count<2){
        }else if((thisVariationalLowerBound - prevVariationalLowerBound) / abs(thisVariationalLowerBound) < _convergenceDiterminationRate){
            _variationalLowerBound = thisVariationalLowerBound;
            break;
        }
        this->updateEta();
        this->updateHyperParameters();
    }
    this->calculateEx();
    // cout<<"iter:"<<count<<endl;
}//}}}
