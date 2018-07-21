//
// VariationalBayesEstimatorOnSLDA.hpp
//
// Copyright (c) 2017 Shion Hosoda
//
// This software is released under the MIT License.
// http://opensource.org/licenses/mit-license.php
//

#ifndef VBSLDA
#define VBSLDA

#include"VariationalBayesEstimatorOnLDA.hpp"
#include"LabelFileParser.hpp"

class VariationalBayesEstimatorOnSLDA:public VariationalBayesEstimatorOnLDA{
protected:
    const unsigned int _C;
    const std::vector<unsigned int> _y;
    std::vector<double> _xi;
    std::vector<std::vector<double> > _eta;
    std::vector<std::vector<double> > _expEtaZ;
public:
    VariationalBayesEstimatorOnSLDA(const BOWFileParser &parser, const LabelFileParser &labelFileParser, const unsigned int K, const double convergenceDiterminationRate);
    ~VariationalBayesEstimatorOnSLDA();
    virtual void initializeParam();
    virtual void updateExpEtaZ(int c);
    virtual void updateExpEtaZ();
    virtual void updateXi();
    virtual void updateQz();
    virtual std::vector<double> calculateNablaL(unsigned int c)const;
    virtual void updateEta();
    virtual double calculateVariationalLowerBound()const;
    virtual void writeParameter(std::string thetaFilename, std::string phiFilename, std::string etaFilename, std::string xiFilename, std::string expEtaZFilename, std::string alphaFilename, std::string betaFilename)const;
    virtual void runIteraions();
};
#endif
