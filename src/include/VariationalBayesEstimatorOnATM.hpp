//
// VariationalBayesEstimatorOnATM.hpp
//
// Copyright (c) 2017 Shion Hosoda
//
// This software is released under the MIT License.
// http://opensource.org/licenses/mit-license.php
//

#ifndef VBATM
#define VBATM

#include<stdlib.h>
#include<math.h>
#include<cmath>
#include<iostream>
#include<vector>
#include<numeric>
#include<memory>
#include<random>
#include<iomanip>
#include<fstream>
#include <boost/math/special_functions/digamma.hpp>
#include <boost/math/special_functions/gamma.hpp>
#include"BOWFileParser.hpp"
#include"AuthorFileParser.hpp"
#include"utils.hpp"


enum BetaUpdateManner{
    SYMMETRY, ASYMMETRY
};

class VariationalBayesEstimatorOnATM{
protected:
    const std::vector<std::vector<unsigned int> > &_frequencyMatrix;
    const std::vector<std::vector<unsigned int> > &_docVoca;
    const unsigned int _K, _V, _D, _M;
    const std::vector<std::vector<unsigned int> > &_docAuth;
    double _convergenceDiterminationRate;
    std::vector<std::vector<std::vector<std::vector<double> > > > _qzy;
    std::vector<std::vector<double> > _thetaEx, _phiEx;
    std::vector<double> _alpha, _beta;
    double _alphaSum, _betaSum;
    std::vector<std::vector<double> > _alphaTimeSeries, _betaTimeSeries;
    std::vector<double> _nm, _nk;
    std::vector<std::vector<double> > _nmk, _nkv;
    double _variationalLowerBound;
    std::vector<double> _VLBTimeSeries;
public:
    VariationalBayesEstimatorOnATM(const BOWFileParser &parser, const AuthorFileParser &aParser, const unsigned int K, const double convergenceDiterminationRate);
    virtual ~VariationalBayesEstimatorOnATM();
    virtual void initializeParam();
    virtual void initializeHyperParam();
    virtual void calculateEx();
    virtual void calculateHyperParamSum();
    virtual void updateQzy();
    virtual void updateNEx();
    virtual void updateBeta(BetaUpdateManner manner);
    virtual void updateHyperParameters();
    virtual double calculateVariationalLowerBound()const;
    virtual void writeParameter(std::string thetaFilename, std::string phiFilename, std::string alphaFilename, std::string betaFilename)const;
    virtual void writeVariationalLowerBound(std::string VLBFilename, std::string VLBTimeSeriesFilename)const;
    virtual void runIteraions();
};

#endif
