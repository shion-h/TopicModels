//
// LDAWBICMetropolisSampler.cpp
//
// Copyright (c) 2017 Shion Hosoda
//
// This software is released under the MIT License.
// http://opensource.org/licenses/mit-license.php
//

#include "include/LDAWBICMetropolisSampler.hpp"


vector<vector<double> > prod(const vector<vector<double> > &theta, const vector<vector<double> > &phi){
    unsigned int D=theta.size(), K=phi.size(), V=phi[0].size();
    vector<vector<double> > ret(D);
    for(int d=0; d<ret.size(); d++)ret[d].assign(V,0);
    for(int d=0; d<D; d++){
        for(int k=0; k<K; k++){
            double thetadk = theta[d][k];
            for(int v=0; v<V; v++){
                ret[d][v] += thetadk * phi[k][v];
            }
        }
    }
    return ret;
}


double calculateLogLikelihood(const vector<vector<double> > &theta, const vector<vector<double> > &phi, const vector<vector<unsigned int> > &BOW){
    double logLikelihood = 0.0;
    unsigned int D=BOW.size(), V=phi[0].size();
    vector<vector<double> > wordProbMat(D);
    for(int d=0;d<wordProbMat.size();d++)wordProbMat[d].assign(V,0);
    wordProbMat = prod(theta, phi);
    for(int d=0; d<D; d++){
        double documentLogLikelihood = 0.0;
        for(int i=0; i<BOW[d].size(); i++){
            documentLogLikelihood += log(wordProbMat[d][BOW[d][i]]);
        }
        logLikelihood += documentLogLikelihood;
    }
    return logLikelihood;
}


double calculateDirichletLogPdf(vector<double> probVar, vector<double> param){
    double pdf = 0;
    unsigned int p = param.size();
    double paramSum = 0.0;
    double term2 = 0.0, term3 = 0.0;
    for(int i=0; i<p; i++){
        paramSum += param[i];
        term2 += boost::math::lgamma(param[i]);
        term3 += (param[i]-1) * log(probVar[i]);
    }
    // cout<<param[0]<<' '<<boost::math::lgamma(param[0])<<' ';
    // cout<<(param[0]-1) * log(probVar[0])<<endl;
    pdf += boost::math::lgamma(paramSum);
    pdf += -term2;
    pdf += term3;
    return pdf;
}


double calculateLogPriorDistributionValue(const vector<vector<double> > &theta, const vector<vector<double> > &phi, const vector<double> alpha, const vector<double> beta){
    double logPriorDistributionValue = 0.0;
    for(int d=0; d<theta.size(); d++){
        logPriorDistributionValue += calculateDirichletLogPdf(theta[d], alpha);
    }
        // cout<<"thetapdf:"<<calculateDirichletLogPdf(theta[0], alpha)<<endl;
    for(int k=0; k<phi.size(); k++){
        logPriorDistributionValue += calculateDirichletLogPdf(phi[k], beta);
    }
    return logPriorDistributionValue;
}


double calculateLogProposalDistributionValue(const vector<vector<double> > &thetaDestination, const vector<vector<double> > &phiDestination, const vector<vector<double>> alphaSource, const vector<vector<double>> betaSource){
    double logPriorDistributionValue = 0.0;
    for(int d=0; d<thetaDestination.size(); d++){
        logPriorDistributionValue += calculateDirichletLogPdf(thetaDestination[d], alphaSource[d]);
    }
    for(int k=0; k<phiDestination.size(); k++){
        logPriorDistributionValue += calculateDirichletLogPdf(phiDestination[k], betaSource[k]);
    }
    return logPriorDistributionValue;
}


double calculateLogTargetDistributionValue(const double &logLikelihood, const double &logPriorDistributionValue, const unsigned int &n){
    double inverseTemperature = 1.0 / log(n);
    double logTargetDistributionValue = inverseTemperature * logLikelihood + logPriorDistributionValue;
    return logTargetDistributionValue;
}


vector<vector<double> > samplingParamFromDirichlet(const vector<vector<double> > &param, const unsigned int &N){
    //TODO:
    // unsigned int N = param.size();
    unsigned int p = param[0].size();
    vector<vector<double> > ret(N);
    for(int n=0;n<ret.size();n++)ret[n].assign(p,0);
    random_device rnd;
    for(int n=0; n<N; n++){
        vector<double> sample;
        for(int i=0; i<p; i++){
            gamma_distribution<double> gamma(param[n][i],1);
            sample.push_back(gamma(rnd));
        }
        double sampleSum=accumulate(sample.begin(),sample.end(), 0.0);
        for(int i=0; i<p; i++){
            ret[n][i] = sample[i] / sampleSum;
        }
    }
    return ret;
}


double runWBICMetropolis(const vector<vector<unsigned int> > &BOW, const vector<double> &alpha, const vector<double> &beta, unsigned int n){
    unsigned int D = BOW.size(), K = alpha.size(), V = beta.size();
    if(n==0){
        for(int d=0; d<D; d++){
            n += BOW[d].size();
        }
    }
    vector<vector<double> > alphaMatrix(D, alpha);
    vector<vector<double> > betaMatrix(K, beta);
    vector<vector<double> > alphaMatrixCandidate(D, alpha);
    vector<vector<double> > betaMatrixCandidate(K, beta);
    vector<vector<double> > alphaMatrixUniform(D, alpha);
    vector<vector<double> > betaMatrixUniform(K, beta);
    for(int d=0; d<D; d++){
        for(int k=0; k<K; k++){
            alphaMatrixUniform[d][k] = 1;
        }
    }
    for(int k=0; k<K; k++){
        for(int v=0; v<V; v++){
            betaMatrixUniform[k][v] = 1;
        }
    }
    unsigned int iteration = 60000;
    unsigned int burnIn = 50000;
    unsigned int samplingInterval = 100;
    const unsigned int RATE = 1000;
    unsigned int samplingTimes = (iteration - burnIn) / samplingInterval;
    unsigned int acceptanceTimes = 0;
    double logLikelihoodAverage = 0.0;
    vector<vector<double> > theta = samplingParamFromDirichlet(alphaMatrixUniform, D);
    vector<vector<double> > phi = samplingParamFromDirichlet(betaMatrixUniform, K);
    double logLikelihood = calculateLogLikelihood(theta, phi, BOW);
    double logPriorDistributionValue = calculateLogPriorDistributionValue(theta, phi, alpha, beta);
    double logTargetDistributionValue = calculateLogTargetDistributionValue(logLikelihood, logPriorDistributionValue, n);
    // double logTargetDistributionValue = 0;
    random_device rnd;
    mt19937 mt(rnd());
    uniform_real_distribution<double> randN(0,1);
    const unsigned int alphaSum = RATE * K;
    const unsigned int betaSum = RATE * V;
    for(int d=0; d<D; d++){
        for(int k=0; k<K; k++){
            alphaMatrix[d][k] = theta[d][k] * (alphaSum - K) + 1;
        }
    }
    for(int k=0; k<K; k++){
        for(int v=0; v<V; v++){
            betaMatrix[k][v] = phi[k][v] * (betaSum - V) + 1;
        }
    }
    for(int i=0; i<iteration; i++){
        vector<vector<double> > thetaCandidate = samplingParamFromDirichlet(alphaMatrix, D);
        vector<vector<double> > phiCandidate = samplingParamFromDirichlet(betaMatrix, K);

        for(int d=0; d<D; d++){
            for(int k=0; k<K; k++){
                alphaMatrixCandidate[d][k] = thetaCandidate[d][k] * (alphaSum - K) + 1;
            }
        }
        for(int k=0; k<K; k++){
            for(int v=0; v<V; v++){
                betaMatrixCandidate[k][v] = phiCandidate[k][v] * (betaSum - V) + 1;
            }
        }
        double logLikelihoodCandidate = calculateLogLikelihood(thetaCandidate, phiCandidate, BOW);
        double logPriorDistributionValueCandidate = calculateLogPriorDistributionValue(thetaCandidate, phiCandidate, alpha, beta);
        double logTargetDistributionValueCandidate = calculateLogTargetDistributionValue(logLikelihoodCandidate, logPriorDistributionValueCandidate, n);
        double logProposalDistributionValueToCandidate = calculateLogProposalDistributionValue(thetaCandidate, phiCandidate, alphaMatrix, betaMatrix);
        double logProposalDistributionValueToSample = calculateLogProposalDistributionValue(theta, phi, alphaMatrixCandidate, betaMatrixCandidate);

        double acceptanceProb = exp(logTargetDistributionValueCandidate - logTargetDistributionValue + logProposalDistributionValueToSample - logProposalDistributionValueToCandidate);
        double randomValue = randN(mt);
        if(acceptanceProb > randomValue){
            acceptanceTimes += 1;
            theta = thetaCandidate;
            phi = phiCandidate;
            alphaMatrix = alphaMatrixCandidate;
            betaMatrix = betaMatrixCandidate;
            logLikelihood = logLikelihoodCandidate;
            logPriorDistributionValue = logPriorDistributionValueCandidate;
            logTargetDistributionValue = logTargetDistributionValueCandidate;
        }
        if((i > burnIn-1) && (i-burnIn+1)%samplingInterval == 0){
            logLikelihoodAverage += logLikelihood;
        }
    }
    logLikelihoodAverage /= samplingTimes;
    cout<<"acceptanceTimes: "<<acceptanceTimes<<endl;
    return logLikelihoodAverage;
}


vector<double> readHyperParam(string hyperParamFilename){
    ifstream ifs(hyperParamFilename);
    vector<double> param;
    if(!ifs){
        cout<<"error";
        exit(1);
    }
    string str;
    while(getline(ifs,str)){
        string token;
        istringstream stream(str);
        while(getline(stream,token,'\n')){
            double temp=stof(token);
            param.push_back(temp);
        }
    }
    return param;
}
