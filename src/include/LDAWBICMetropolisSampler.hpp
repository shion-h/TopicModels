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


std::vector<std::vector<double> > prod(const std::vector<std::vector<double> > &theta, const std::vector<std::vector<double> > &phi);

double calculateLogLikelihood(const std::vector<std::vector<double> > &theta, const std::vector<std::vector<double> > &phi, const std::vector<std::vector<unsigned int> > &BOW);

double calculateDirichletLogPdf(std::vector<double> probVar, std::vector<double> param);

double calculateLogPriorValue(const std::vector<std::vector<double> > &theta, const std::vector<std::vector<double> > &phi, const std::vector<double> alpha, const std::vector<double> beta);

double calculateTargetDistributionValue(const double &logLikelihood, const double &logPriorValue, const unsigned int &n);

std::vector<std::vector<double> > samplingParamFromDirichlet(const std::vector<std::vector<double> > &param, const unsigned int &N);

double runWBICMetropolis(const std::vector<std::vector<unsigned int> > &frequencyMatrix, const std::vector<std::vector<unsigned int> > &docVoca, const std::vector<double> &alpha, const std::vector<double> &beta, unsigned int n=0);

std::vector<double> readHyperParam(std::string hyperParamFilename);
