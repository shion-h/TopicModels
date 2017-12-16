//
// VariationalBayesEstimatorOnLDA.hpp
//
// Copyright (c) 2017 Shion Hosoda
//
// This software is released under the MIT License.
// http://opensource.org/licenses/mit-license.php
//

#ifndef INCLUDED_VB
#define INCLUDED_VB

#include<stdlib.h>
#include<math.h>
#include<iostream>
#include<vector>
#include<numeric>
#include<memory>
#include<random>
#include<iomanip>
#include<fstream>
#include <boost/math/special_functions/digamma.hpp>
#include <boost/math/special_functions/gamma.hpp>


class VariationalBayesEstimatorOnLDA{
protected:
    const std::vector<std::vector<unsigned int> > &_bagOfWordsNum;
    const unsigned int _K, _V;
    double _convergenceDiterminationRate;
    std::vector<std::vector<double> > _thetaEx, _phiEx;
    std::vector<double> _alpha, _beta;
    double _alphaSum, _betaSum;
    std::vector<std::vector<double> > _alphaTimeSeries, _betaTimeSeries;
    std::vector<double> _nd, _nk;
    std::vector<std::vector<double> > _ndk, _nkv;
    const std::vector<std::string> &_wordList;
    double _variationalLowerBoundOfQz;
    double _variationalLowerBound;
    std::vector<double> _VLBTimeSeries;
public:
    VariationalBayesEstimatorOnLDA(const std::vector<std::vector<unsigned int > > &bagOfWordsNum, const std::vector<std::string > &wordList, const unsigned int K, const unsigned int V, const double convergenceDiterminationRate);
    virtual ~VariationalBayesEstimatorOnLDA();
    virtual void initializeHyperParam();
    virtual void initializeParam();
    virtual void calculateEx();
    virtual void calculateHyperParamSum();
    virtual std::vector<double> calculateQz(unsigned int d, unsigned int i);
    virtual void nExUpdate();
    virtual void hyperParamUpdate();
    virtual double calculateVariationalLowerBound();
    virtual void printHyperParameter()const;
    virtual void printTopFactor(int numOfTopFactor)const;
    virtual void printThetaEx()const;
    virtual void printPhiEx()const;
    virtual void printNum()const;
    virtual void writeParameter(std::string thetaFilename, std::string phiFilename, std::string alphaFilename, std::string betaFilename)const;
    virtual void writeVariationalLowerBound(std::string VLBFilename, std::string VLBTimeSeriesFilename)const;
    virtual void runIteraions();
};

#endif
