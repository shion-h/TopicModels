//
// ATM.cpp
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
#include "include/AuthorFileParser.hpp"
#include "include/CollapsedGibbsSamplerFromATM.hpp"
//}}}

using namespace std;

int main(int argc, char *argv[]){
    //オプション読み取り{{{
    boost::program_options::options_description opt("Options");
    opt.add_options()
    ("help,h", "show help")
    ("nmtp,k", boost::program_options::value<unsigned int>()->default_value(10), "number of topics")
    ("nmit,s", boost::program_options::value<unsigned int>()->default_value(1000), "number of iterations")
    ("bnin,b", boost::program_options::value<unsigned int>(), "burn in period")
    ("intv,i", boost::program_options::value<unsigned int>()->default_value(10), "sampling interval")
    ("nmsh,f", boost::program_options::value<unsigned int>()->default_value(5), "number of factors with high probability to show")
    ("otpt,o", boost::program_options::value<string>()->default_value("./"), "directory name for output")
    ;

    boost::program_options::positional_options_description pd;
    pd.add("BOWfile",1);
    pd.add("authorfile",2);

    boost::program_options::options_description hidden("処理用");
    hidden.add_options()
        ("BOWfile",boost::program_options::value<string>(), "hidden")
        ("authorfile",boost::program_options::value<string>(), "hidden")
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
    string authorFilename;
    unsigned int K;
    unsigned int S;
    unsigned int burnIn;
    unsigned int samplingInterval;
    unsigned int numOfTopFactor;
    string outputDirectory;
    string thetaFilename;
    string phiFilename;
    string alphaFilename;
    string betaFilename;
    string wordListFilename;
    string authorListFilename;
    if (vm.count("help") || !vm.count("BOWfile") || !vm.count("authorfile")){
        cout<<"Usage:\n ATM [BOW file] [Author file] [-options] "<<endl;
        cout<<endl;
        cout<<opt<<endl;
        exit(1);
    }else{
        BOWFilename=vm["BOWfile"].as<std::string>();
        authorFilename=vm["authorfile"].as<std::string>();
        if(vm.count("nmtp"))K=vm["nmtp"].as<unsigned int>();
        if(vm.count("nmit"))S=vm["nmit"].as<unsigned int>();
        if(vm.count("bnin"))burnIn=vm["bnin"].as<unsigned int>();
        else burnIn=static_cast<int>(S*0.8);
        if(vm.count("intv"))samplingInterval=vm["intv"].as<unsigned int>();
        if(vm.count("nmsh"))numOfTopFactor=vm["nmsh"].as<unsigned int>();
        if(vm.count("otpt"))outputDirectory=vm["otpt"].as<std::string>();
        if(outputDirectory[outputDirectory.size()-1] != '/')outputDirectory.push_back('/')
        thetaFilename=outputDirectory+"theta.csv";
        phiFilename=outputDirectory+"phi.csv";
        alphaFilename=outputDirectory+"alpha.csv";
        betaFilename=outputDirectory+"beta.csv";
        wordListFilename=outputDirectory+"wordList.csv";
        authorListFilename=outputDirectory+"authorList.csv";
    }

    BOWFileParser parser(BOWFilename);
    parser.readBOWFile();
    parser.makeBagOfWordsNum();
    parser.writeWordList(wordListFilename);
    unsigned int V=parser.getV();
    vector<vector<unsigned int> > bagOfWordsNum=parser.getBagOfWordsNum();

    cout<<"filename="<<BOWFilename<<' ';
    cout<<"author filename="<<authorFilename<<' ';
    cout<<"K="<<K<<' ';
    cout<<"V="<<V<<' ';
    cout<<"S="<<S<<' ';
    cout<<"burnIn="<<burnIn<<' ';
    cout<<"samplingInterval="<<samplingInterval<<' ';
    cout<<"thetaFilename="<<thetaFilename<<' ';
    cout<<"phiFilename="<<phiFilename<<' ';
    cout<<"wordListFilename="<<wordListFilename<<' ';
    cout<<"authorListFilename="<<authorListFilename<<endl;
    //}}}

//サンプリング{{{

    vector<unsigned int> nThIterationSample;
    for(int i=0;i<((S-burnIn)/samplingInterval);i++){
        nThIterationSample.push_back(burnIn+samplingInterval*i);
    }
    AuthorFileParser authorFileParser(authorFilename);
    authorFileParser.readAuthorFile();
    authorFileParser.makeBagOfAuthorsNum();
    authorFileParser.writeAuthorList(authorListFilename);
    unsigned int A=authorFileParser.getY();
    vector<vector<unsigned int> > bagOfAuthorsNum=authorFileParser.getBagOfAuthorsNum();
    // shared_ptr<GibbsSamplerFromATM> gibbsSampler;
    // gibbsSampler=make_shared<CollapsedGibbsSamplerFromATM>(bagOfWordsNum,bagOfAuthorsNum,K,V,A,nThIterationSample);
    GibbsSamplerFromATM *gibbsSampler=new CollapsedGibbsSamplerFromATM(bagOfWordsNum,bagOfAuthorsNum,parser.getWordList(),authorFileParser.getAuthorList(),K,V,A,nThIterationSample);
    for(int s=0;s<S;s++){
        gibbsSampler->runIteraion();
    }
    gibbsSampler->computeParameter();
    gibbsSampler->writeParameter(thetaFilename,phiFilename,alphaFilename,betaFilename);
    gibbsSampler->printTopFactor(numOfTopFactor);
    delete gibbsSampler;
//}}}
    return 0;
}
