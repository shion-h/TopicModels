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
#include<boost/math/special_functions/gamma.hpp>
#include"BOWFileParser.hpp"

using namespace std;


vector<vector<double> > prod(const vector<vector<double> > &theta, const vector<vector<double> > &phi);

double calculateLogLikelihood(const vector<vector<double> > &theta, const vector<vector<double> > &phi, const vector<vector<unsigned int> > &BOW);

double calculateDirichletLogPdf(vector<double> probVar, vector<double> param);

double calculateLogPriorValue(const vector<vector<double> > &theta, const vector<vector<double> > &phi, const vector<double> alpha, const vector<double> beta);

double calculateTargetDistributionValue(const double &logLikelihood, const double &logPriorValue, const unsigned int &n);

vector<vector<double> > samplingParamFromDirichlet(const vector<vector<double> > &param, const unsigned int &N);

double runWBICMetropolis(const vector<vector<unsigned int> > &BOW, const vector<double> &alpha, const vector<double> &beta, unsigned int n=0);

vector<double> readHyperParam(string hyperParamFilename);
