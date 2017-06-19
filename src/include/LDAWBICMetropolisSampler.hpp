//
// LDAWBICMetropolisSampler.hpp
//
// Copyright (c) 2017 Shion Hosoda
//
// This software is released under the MIT License.
// http://opensource.org/licenses/mit-license.php
//

#include<random>
#include<math.h>
#include <boost/math/special_functions/gamma.hpp>
#include "include/BOWFileParser.hpp"
// #include <boost/python.hpp>
// #include <boost/python/suite/indexing/vector_indexing_suite.hpp>

using namespace std;
// typedef std::vector<std::vector<double> > dmatrix;


vector<vector<double> > prod(const vector<vector<double> > &theta, const vector<vector<double> > &phi);


double calculateLogLikelihood(const vector<vector<double> > &theta, const vector<vector<double> > &phi, const vector<vector<unsigned int> > &BOW);


double calculateDirichletLogPdf(vector<double> probVar, vector<double> param);


double calculateLogPriorValue(const vector<vector<double> > &theta, const vector<vector<double> > &phi, const vector<double> alpha, const vector<double> beta);


double calculateTargetDistributionValue(const double &logLikelihood, const double &logPriorValue, const unsigned int &n);


vector<vector<double> > samplingParamFromDirichlet(const vector<vector<double> > &param, const unsigned int &N);


double runWBICMetropolis(const vector<vector<unsigned int> > &BOW, const vector<double> &alpha, const vector<double> &beta, unsigned int n);

vector<double> readHyperParam(string hyperParamFilename);
