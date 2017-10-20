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
#include "include/OrdinaryGibbsSamplerFromLDA.hpp"
#include "include/CollapsedGibbsSamplerFromLDA.hpp"
#include "include/CollapsedGibbsSamplerFromATM.hpp"
#include "include/VariationalBayesEstimatorOnLDA.hpp"
//}}}

using namespace std;

int main(int argc, char *argv[]){

    //Options{{{
    boost::program_options::options_description opt("Options");
    opt.add_options()
    ("help,h", "show help")
    ("nmtp,k", boost::program_options::value<unsigned int>()->default_value(10), "number of topics")
    ("nmit,s", boost::program_options::value<unsigned int>()->default_value(1000), "number of iterations(only used by sampling algorythm)")
    ("bnin,b", boost::program_options::value<unsigned int>(), "burn in period(only used by sampling algorythm)")
    ("intv,i", boost::program_options::value<unsigned int>()->default_value(10), "sampling interval(only used by sampling algorythm)")
    ("cdrt,d", boost::program_options::value<double>()->default_value(0.001), "convergence ditermination rate(only used by VB algorythm)")
    ("lrna,l", boost::program_options::value<unsigned int>()->default_value(1), "learning algorythm(0:Gibbs sampling 1:Collapsed gibbs sampling 2:Variational Bayes)")
    ("nmsh,f", boost::program_options::value<unsigned int>()->default_value(5), "number of factors with high probability to show")
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
    }catch(const boost::program_options::error_with_option_name& e){
        cout<<e.what()<<endl;
    }
    boost::program_options::notify(vm);
    string BOWFilename;
    unsigned int K;
    unsigned int S;
    unsigned int burnIn;
    unsigned int samplingInterval;
    double convergenceDiterminationRate;
    unsigned int learningAlgorythmFlag;
    unsigned int numOfTopFactor;
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
        if(vm.count("nmit"))S = vm["nmit"].as<unsigned int>();
        if(vm.count("bnin"))burnIn = vm["bnin"].as<unsigned int>();
        else burnIn = static_cast<int>(S*0.8);
        if(vm.count("intv"))samplingInterval = vm["intv"].as<unsigned int>();
        if(vm.count("cdrt"))convergenceDiterminationRate = vm["cdrt"].as<double>();
        if(vm.count("lrna"))learningAlgorythmFlag = vm["lrna"].as<unsigned int>();
        if(vm.count("nmsh"))numOfTopFactor = vm["nmsh"].as<unsigned int>();
        if(vm.count("otpt"))outputDirectory = vm["otpt"].as<std::string>();
        cout<<outputDirectory<<endl;
        thetaFilename = outputDirectory + "theta.csv";
        phiFilename = outputDirectory + "phi.csv";
        alphaFilename = outputDirectory + "alpha.csv";
        betaFilename = outputDirectory + "beta.csv";
        VLBFilename = outputDirectory + "variationalLowerBound.csv";
        VLBTimeSeriesFilename = outputDirectory + "VLBTimeSeries.csv";
        wordListFilename = outputDirectory + "wordList.csv";
    }

    BOWFileParser parser(BOWFilename);
    parser.readBOWFile();
    parser.makeBagOfWordsNum();
    parser.writeWordList(wordListFilename);
    unsigned int V = parser.getV();
    vector<vector<unsigned int> > bagOfWordsNum = parser.getBagOfWordsNum();

    cout<<"filename = "<<BOWFilename<<' ';
    cout<<"K = "<<K<<' ';
    cout<<"V = "<<V<<' ';
    cout<<"learningAlgorythmFlag = "<<learningAlgorythmFlag<<' ';
    if(learningAlgorythmFlag == 0 || learningAlgorythmFlag == 1){
        cout<<"iterationTimes = "<<S<<' ';
        cout<<"burnIn = "<<burnIn<<' ';
        cout<<"samplingInterval = "<<samplingInterval<<' ';
    }else if(learningAlgorythmFlag == 2){
        cout<<"convergenceDiterminationRate = "<<convergenceDiterminationRate<<' ';
    }
    cout<<"thetaFilename = "<<thetaFilename<<' ';
    cout<<"phiFilename = "<<phiFilename<<' ';
    cout<<"alphaFilename = "<<alphaFilename<<' ';
    cout<<"betaFilename = "<<betaFilename<<' ';
    cout<<"wordListFilename = "<<wordListFilename<<endl;
    //}}}

//estimation{{{
    if(learningAlgorythmFlag == 2){
        VariationalBayesEstimatorOnLDA *estimator;
        estimator = new VariationalBayesEstimatorOnLDA(bagOfWordsNum, parser.getWordList(), K, V, convergenceDiterminationRate);
        estimator->runIteraions();
        estimator->writeParameter(thetaFilename, phiFilename, alphaFilename, betaFilename);
        estimator->writeVariationalLowerBound(VLBFilename, VLBTimeSeriesFilename);
        estimator->printTopFactor(numOfTopFactor);
        delete estimator;
    }else{
        vector<unsigned int> nThIterationSample;
        for(int i=0;i < ((S-burnIn)/samplingInterval);i++){
            nThIterationSample.push_back(burnIn + samplingInterval*i);
        }
        // shared_ptr<GibbsSamplerFromLDA> gibbsSampler;
        GibbsSamplerFromLDA *gibbsSampler;
        switch(learningAlgorythmFlag){
            case 0:
                // gibbsSampler = make_shared<OrdinaryGibbsSamplerFromLDA>(bagOfWordsNum, K, V, nThIterationSample);
                gibbsSampler = new OrdinaryGibbsSamplerFromLDA(bagOfWordsNum, parser.getWordList(), K, V, nThIterationSample);
                break;
            case 1:
                // gibbsSampler = make_shared<CollapsedGibbsSamplerFromLDA>(bagOfWordsNum, K, V, nThIterationSample);
                gibbsSampler = new CollapsedGibbsSamplerFromLDA(bagOfWordsNum, parser.getWordList(), K, V, nThIterationSample);
                break;
        }
        for(int s=0;s<S;s++){
            gibbsSampler->runIteraion();
        }
        gibbsSampler->computeParameter();
        gibbsSampler->writeParameter(thetaFilename, phiFilename, alphaFilename, betaFilename);
        gibbsSampler->printTopFactor(numOfTopFactor);
        delete gibbsSampler;
    }
//}}}
    return 0;
}
