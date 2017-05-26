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

using namespace std;


class VariationalBayesEstimatorOnLDA{
protected:
    const vector<vector<unsigned int> > &bagOfWordsNum;
    const unsigned int K,V;
    vector<vector<double> > thetaEx,phiEx;
    vector<double> alpha,beta;
    double alphaSum,betaSum;
    vector<double> nd, nk;
    vector<vector<double> > ndk, nkv;
    const vector<string> &_wordList;
    double _variationalLowerBoundOfQz;
    double _variationalLowerBound;
public:
    VariationalBayesEstimatorOnLDA(const vector<vector<unsigned int > > &bagOfWordsNum,const vector<string > &wordList,const unsigned int K,const unsigned int V);
    virtual ~VariationalBayesEstimatorOnLDA();
    virtual void initializeHyperParam();
    virtual void initializeParam();
    virtual void calculateEx();
    virtual void calculateHyperParamSum();
    virtual vector<double> calculateQz(unsigned int d, unsigned int i);
    virtual void nExUpdate();
    virtual void hyperParamUpdate();
    virtual double calculateVariationalLowerBound();
    virtual void printHyperParameter()const;
    virtual void printTopFactor(int numOfTopFactor)const;
    virtual void printThetaEx()const;
    virtual void printPhiEx()const;
    virtual void printNum()const;
    virtual void writeParameter(string thetaFilename,string phiFilename,string alphaFilename,string betaFilename)const;
    virtual void writeThetaEx(string thetaFilename)const;
    virtual void writePhiEx(string phiFilename)const;
    virtual void writeVariationalLowerBound(string VLBFilename)const;
    virtual void runIteraions();
};

#endif
