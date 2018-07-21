//
// HDP.cpp
//
// Copyright (c) 2018 Shion Hosoda
//
// This software is released under the MIT License.
// http://opensource.org/licenses/mit-license.php
//

//include{{{
#include <stdlib.h>
#include <boost/program_options.hpp>
#include "include/BOWFileParser.hpp"
#include "include/GibbsSamplerFromHDP.hpp"
//}}}

using namespace std;

int main(int argc, char *argv[]){

    //Options{{{
    boost::program_options::options_description opt("Options");
    opt.add_options()
    ("help,h", "show help")
    ("otpt,o", boost::program_options::value<string>()->default_value("./"), "directory name for output")
    ("itnm,n", boost::program_options::value<unsigned int>()->default_value(1000), "the number of iteration")
    ("intr,i", boost::program_options::value<unsigned int>()->default_value(5), "sampling interval")
    ("bnin,b", boost::program_options::value<unsigned int>()->default_value(500), "burn-in term")
    ("alph,a", boost::program_options::value<double>()->default_value(100), "alpha value")
    ("beta,c", boost::program_options::value<double>()->default_value(100), "beta value")
    ("gamm,g", boost::program_options::value<double>()->default_value(100), "gamma value")
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
    double alpha, beta, gamma;
    unsigned int iterationNumber;
    unsigned int samplingInterval;
    unsigned int burnIn;
    string outputDirectory;
    string thetaFilename;
    string phiFilename;
    string logLikelihoodFilename;
    string wordListFilename;
    if (vm.count("help") || !vm.count("BOWfile")){
        cout<<"Usage:\n HDP [BOW file] [-options] "<<endl;
        cout<<endl;
        cout<<opt<<endl;
        exit(1);
    }else{
        BOWFilename = vm["BOWfile"].as<std::string>();
        if(vm.count("alph"))alpha = vm["alph"].as<double>();
        if(vm.count("beta"))beta = vm["beta"].as<double>();
        if(vm.count("gamm"))gamma = vm["gamm"].as<double>();
        if(vm.count("itnm"))iterationNumber = vm["itnm"].as<unsigned int>();
        if(vm.count("intr"))samplingInterval = vm["intr"].as<unsigned int>();
        if(vm.count("bnin"))burnIn = vm["bnin"].as<unsigned int>();
        if(vm.count("otpt"))outputDirectory = vm["otpt"].as<std::string>();
        if(outputDirectory[outputDirectory.size()-1] != '/')outputDirectory.push_back('/');
        thetaFilename = outputDirectory + "theta.csv";
        phiFilename = outputDirectory + "phi.csv";
        logLikelihoodFilename = outputDirectory + "logLikelihood.csv";
        wordListFilename = outputDirectory + "wordList.csv";
    }

    BOWFileParser parser(BOWFilename);
    parser.writeWordList(wordListFilename);
    //}}}

//estimation{{{
    GibbsSamplerFromHDP *sampler;
    sampler = new GibbsSamplerFromHDP(parser, alpha, beta, gamma, iterationNumber, burnIn, samplingInterval);
    sampler->runIteraions();
    // sampler->writeParameters(thetaFilename, phiFilename);
    // sampler->writeLogLikelihood(logLikelihoodFilename);
    // delete sampler;
//}}}
    return 0;
}
