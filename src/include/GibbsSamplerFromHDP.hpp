//
// GibbsSamplerFromHDP.hpp
//
// Copyright (c) 2018 Shion Hosoda
//
// This software is released under the MIT License.
// http://opensource.org/licenses/mit-license.php
//

#ifndef HDP
#define HDP

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
#include<limits>
#include <boost/math/special_functions/digamma.hpp>
#include <boost/math/special_functions/gamma.hpp>
#include"BOWFileParser.hpp"
#include"utils.hpp"


class GibbsSamplerFromHDP{
protected:
    const std::vector<std::vector<unsigned int> > &_frequencyMatrix;
    const std::vector<std::vector<unsigned int> > &_docVoca;
    const unsigned int _J, _V, _iterationNumber, _burnIn, _samplingInterval;
    const double _alpha, _beta, _gamma;
    unsigned int _K, _samplingCount;
    std::vector<unsigned int> _Tj;
    std::vector<std::vector<unsigned int> >  _xji, _tji, _kjt;
    std::vector<std::vector<double> > _thetaEx, _phiEx;
    std::vector<unsigned int> _nj, _nk, _mk;
    std::vector<double> _sumOfNkSampled;
    std::vector<std::vector<unsigned int> > _njt, _nkv;
    std::vector<std::vector<double> > _sumOfNjkSampled, _sumOfNkvSampled;
    std::vector<std::vector<std::vector<unsigned int> > > _njtv;
    std::vector<double> _logLikelihood;
public:
    GibbsSamplerFromHDP(const BOWFileParser &parser, double alpha, double beta, double gamma, unsigned int iterationNumber, unsigned int burnIn, unsigned int samplingInterval);
    ~GibbsSamplerFromHDP();
    virtual void initializeParameters();
    virtual void sampleTji();
    virtual void sampleKjt(unsigned int j, unsigned int t, bool isNewT=false);
    virtual void countN();
    virtual void addT(unsigned int j, unsigned int i);
    virtual void addK();
    virtual void deleteT(unsigned int j, unsigned int t);
    virtual void deleteK(unsigned int k);
    virtual void calculateExpectations(bool isUseOnlyLastSample=false);
    virtual void calculateLogLikelihood();
    virtual void checkNumbers();
    virtual void writeParameters(std::string thetaFilename, std::string phiFilename)const;
    virtual void writeLogLikelihood(std::string logLikelihoodFilename)const;
    virtual void storeSamples();
    virtual void runIteraions();
};

#endif
