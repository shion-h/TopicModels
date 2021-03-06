//
// sLDA.cpp
//
// Copyright (c) 2017 Shion Hosoda
//
// This software is released under the MIT License.
// http://opensource.org/licenses/mit-license.php
//

//include{{{
#include <stdlib.h>
#include <boost/program_options.hpp>
#include "include/BOWFileParser.hpp"
#include "include/LabelFileParser.hpp"
#include "include/MultiLabelFileParser.hpp"
#include "include/VariationalBayesEstimatorOnSLDA.hpp"
//}}}

using namespace std;

int main(int argc, char *argv[]){

    //Options{{{
    boost::program_options::options_description opt("Options");
    opt.add_options()
    ("help,h", "show help")
    ("nmtp,k", boost::program_options::value<unsigned int>()->default_value(10), "number of topics")
    ("cdrt,d", boost::program_options::value<double>()->default_value(0.001), "convergence ditermination rate")
    ("nmsh,f", boost::program_options::value<unsigned int>()->default_value(5), "number of factors with high probability to show")
    ("otpt,o", boost::program_options::value<string>()->default_value("./"), "directory name for output")
    ;

    boost::program_options::positional_options_description pd;
    pd.add("BOWfile", 1);
    pd.add("labelfile", 2);

    boost::program_options::options_description hidden("hidden");
    hidden.add_options()
        ("BOWfile", boost::program_options::value<string>(), "hidden")
        ("labelfile", boost::program_options::value<string>(), "hidden")
        ;
    boost::program_options::options_description cmdline_options;
    cmdline_options.add(opt).add(hidden);

    boost::program_options::variables_map vm;
    try{
        boost::program_options::store(boost::program_options::command_line_parser(argc, argv).options(cmdline_options).positional(pd).run(), vm);
    }catch(const boost::program_options::error_with_option_name& e){
        cout<<e.what()<<endl;
    }
    boost::program_options::notify(vm);
    string BOWFilename;
    string labelFilename;
    unsigned int K;
    double convergenceDiterminationRate;
    unsigned int numOfTopFactor;
    string outputDirectory;
    string thetaFilename;
    string phiFilename;
    string etaFilename;
    string xiFilename;
    string expEtaZFilename;
    string alphaFilename;
    string betaFilename;
    string VLBFilename;
    string VLBTimeSeriesFilename;
    string wordListFilename;
    string labelListFilename;
    if (vm.count("help") || !vm.count("BOWfile") || !vm.count("labelfile")){
        cout<<"Usage:\n LDA [BOW file] [label file] [-options] "<<endl;
        cout<<endl;
        cout<<opt<<endl;
        exit(1);
    }else{
        BOWFilename = vm["BOWfile"].as<std::string>();
        labelFilename = vm["labelfile"].as<std::string>();
        if(vm.count("nmtp"))K = vm["nmtp"].as<unsigned int>();
        if(vm.count("cdrt"))convergenceDiterminationRate = vm["cdrt"].as<double>();
        if(vm.count("nmsh"))numOfTopFactor = vm["nmsh"].as<unsigned int>();
        if(vm.count("otpt"))outputDirectory = vm["otpt"].as<std::string>();
        if(outputDirectory[outputDirectory.size()-1] != '/')outputDirectory.push_back('/');
        thetaFilename = outputDirectory + "theta.csv";
        phiFilename = outputDirectory + "phi.csv";
        etaFilename = outputDirectory + "eta.csv";
        xiFilename = outputDirectory + "xi.csv";
        expEtaZFilename = outputDirectory + "expEtaZ.csv";
        alphaFilename = outputDirectory + "alpha.csv";
        betaFilename = outputDirectory + "beta.csv";
        VLBFilename = outputDirectory + "variationalLowerBound.csv";
        VLBTimeSeriesFilename = outputDirectory + "VLBTimeSeries.csv";
        wordListFilename = outputDirectory + "wordList.csv";
        labelListFilename = outputDirectory + "labelList.csv";
    }

    BOWFileParser parser(BOWFilename);
    parser.writeWordList(wordListFilename);

    LabelFileParser labelFileParser(labelFilename);
    labelFileParser.writeLabelList(labelListFilename);
    //}}}

//estimation{{{
    VariationalBayesEstimatorOnSLDA *estimator;
    estimator = new VariationalBayesEstimatorOnSLDA(parser, labelFileParser, K, convergenceDiterminationRate);
    estimator->runIteraions();
    estimator->writeParameter(thetaFilename, phiFilename, etaFilename, xiFilename, expEtaZFilename, alphaFilename, betaFilename);
    estimator->writeVariationalLowerBound(VLBFilename, VLBTimeSeriesFilename);
    delete estimator;
//}}}
    return 0;
}
