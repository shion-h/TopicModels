//
// LDA.cpp
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
#include "include/VariationalBayesEstimatorOnLDA.hpp"
//}}}

using namespace std;

int main(int argc, char *argv[]){

    //Options{{{
    boost::program_options::options_description opt("Options");
    opt.add_options()
    ("help,h", "show help")
    ("nmtp,k", boost::program_options::value<unsigned int>()->default_value(10), "number of topics")
    ("cdrt,d", boost::program_options::value<double>()->default_value(0.001), "convergence ditermination rate")
    ("otpt,o", boost::program_options::value<string>()->default_value("./"), "directory name for output")
    ;

    boost::program_options::positional_options_description pd;
    pd.add("BOWfile", 1);

    boost::program_options::options_description hidden("hidden");
    hidden.add_options()
        ("BOWfile", boost::program_options::value<string>(), "hidden")
        ;
    boost::program_options::options_description cmdline_options;
    cmdline_options.add(opt).add(hidden);

    boost::program_options::variables_map vm;
    try{
        boost::program_options::store(boost::program_options::command_line_parser(argc, argv).options(cmdline_options).positional(pd).run(), vm);
    // }catch(const boost::program_options::error_with_option_name& e){
    //     cout<<e.what()<<endl;
    }catch(...){
        cout<<"Option error."<<endl;
    }
    boost::program_options::notify(vm);
    string BOWFilename;
    unsigned int K;
    double convergenceDiterminationRate;
    string outputDirectory;
    string thetaFilename;
    string phiFilename;
    string alphaFilename;
    string betaFilename;
    string VLBFilename;
    string VLBTimeSeriesFilename;
    string wordListFilename;
    if (vm.count("help") || !vm.count("BOWfile")){
        cout<<"Usage:\n LDA [BOW file] [-options] "<<endl;
        cout<<endl;
        cout<<opt<<endl;
        exit(1);
    }else{
        BOWFilename = vm["BOWfile"].as<std::string>();
        if(vm.count("nmtp"))K = vm["nmtp"].as<unsigned int>();
        if(vm.count("cdrt"))convergenceDiterminationRate = vm["cdrt"].as<double>();
        if(vm.count("otpt"))outputDirectory = vm["otpt"].as<std::string>();
        if(outputDirectory[outputDirectory.size()-1] != '/')outputDirectory.push_back('/');
        thetaFilename = outputDirectory + "theta.csv";
        phiFilename = outputDirectory + "phi.csv";
        alphaFilename = outputDirectory + "alpha.csv";
        betaFilename = outputDirectory + "beta.csv";
        VLBFilename = outputDirectory + "variationalLowerBound.csv";
        VLBTimeSeriesFilename = outputDirectory + "VLBTimeSeries.csv";
        wordListFilename = outputDirectory + "wordList.csv";
    }

    BOWFileParser parser(BOWFilename);
    parser.writeWordList(wordListFilename);
    //}}}

//estimation{{{
    VariationalBayesEstimatorOnLDA *estimator;
    estimator = new VariationalBayesEstimatorOnLDA(parser, K, convergenceDiterminationRate);
    estimator->runIteraions();
    estimator->writeParameter(thetaFilename, phiFilename, alphaFilename, betaFilename);
    estimator->writeVariationalLowerBound(VLBFilename, VLBTimeSeriesFilename);
    delete estimator;
//}}}
    return 0;
}
