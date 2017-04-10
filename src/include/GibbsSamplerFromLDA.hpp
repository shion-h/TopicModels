//
// GibbsSamplerFromLDA.hpp
//
// Copyright (c) 2017 Shion Hosoda
//
// This software is released under the MIT License.
// http://opensource.org/licenses/mit-license.php
//

#ifndef INCLUDED_GSTEMP
#define INCLUDED_GSTEMP

#include<stdlib.h>
#include<iostream>
#include<vector>
#include<numeric>
#include<memory>
#include<random>
#include<iomanip>
#include<fstream>
#include <boost/math/special_functions/digamma.hpp>

using namespace std;

struct LDA{
    const vector<vector<unsigned int> > &bagOfWordsNum;
    const unsigned int K,V;
    vector<vector<unsigned int > > z;
    vector<vector<double> > theta,phi;
    vector<vector<double> > thetaEx,phiEx;
    vector<double> alpha,beta;
    double alphaSum,betaSum;
    vector<unsigned int > topicCounts;
    vector<vector<unsigned int > > topicCountsInDoc,vocaCountsFromTopic;

    LDA(const vector<vector<unsigned int> > &bagOfWordsNum,const unsigned int K,const unsigned int V);
    virtual void initializeHyperParam();
    virtual void initializeParam();
    virtual void calculateEx();
    virtual void calculateHyperParamSum();
    virtual void nUpdate();
};

template<typename MODEL>
class BaseGibbsSampler{
protected:
    shared_ptr<MODEL> _model;
    vector<vector<double> > _thetaAve,_phiAve;
    unsigned int _iteration;
    const vector<unsigned int> _nThIterationSample;
    const vector<string> &_wordList;
public:
    BaseGibbsSampler(const vector<vector<unsigned int > > &bagOfWordsNum,const vector<string > &wordList,const unsigned int K,const unsigned int V,const vector<unsigned int> nThIterationSample);
    BaseGibbsSampler(const vector<string> &wordList,const vector<unsigned int> nThIterationSample);
    virtual ~BaseGibbsSampler()=0;
    virtual void hyperParamUpdate();
    virtual void computeParameter();
    virtual void printTopFactor(int numOfTopFactor)const;
    virtual void printTheta()const;
    virtual void printPhi()const;
    virtual void printThetaEx()const;
    virtual void printPhiEx()const;
    virtual void printNum()const;
    virtual void printHyperParameter()const;
    virtual void writeParameter(string thetaFilename,string phiFilename,string alphaFilename,string betaFilename)const;
    virtual void writeTheta(string thetaFilename)const;
    virtual void writePhi(string phiFilename)const;
    virtual void samplingParameter();
    virtual void runIteraion()=0;
};

template<typename MODEL>
BaseGibbsSampler<MODEL>::BaseGibbsSampler(const vector<string> &wordList,const vector<unsigned int> nThIterationSample)//{{{
    :_iteration(0),
     _wordList(wordList),
     _nThIterationSample(nThIterationSample){
}//}}}
template<typename MODEL>
BaseGibbsSampler<MODEL>::BaseGibbsSampler(const vector<vector<unsigned int> > &bagOfWordsNum,const vector<string> &wordList,const unsigned int K,const unsigned int V,const vector<unsigned int> nThIterationSample)//{{{
    :_thetaAve(bagOfWordsNum.size()),
     _phiAve(K),
     _iteration(0),
     _wordList(wordList),
     _nThIterationSample(nThIterationSample){
    _model=make_shared<MODEL>(bagOfWordsNum,K,V);
    for(int d=0;d<_thetaAve.size();d++)_thetaAve[d].assign(_model->K,0);
    for(int k=0;k<_phiAve.size();k++)_phiAve[k].assign(_model->V,0);
}//}}}
template<typename MODEL>
BaseGibbsSampler<MODEL>::~BaseGibbsSampler(){//{{{

}//}}}
template<typename MODEL>
void BaseGibbsSampler<MODEL>::hyperParamUpdate(){//{{{
    _model->calculateEx();
    double numerator=0,denominator=0;

    for(int k=0;k<_model->K;k++){
        for(int d=0;d<_model->z.size();d++){
            numerator+=(boost::math::digamma(_model->thetaEx[d][k]*_model->z[d].size()+_model->alpha[k])-boost::math::digamma(_model->alpha[k]))*_model->alpha[k];
            denominator+=boost::math::digamma(_model->z[d].size()+_model->alphaSum)-boost::math::digamma(_model->alphaSum);
        }
        _model->alpha[k]=numerator/denominator;
    }
    double betaUpdated;
    numerator=0;
    denominator=0;
    for(int k=0;k<_model->K;k++){
        double numOfTopicsExpected=0;
        for(int v=0;v<_model->V;v++){
            double numOfVocasfromTopicExpected=0;
            for(int d=0;d<_model->z.size();d++){
                numOfVocasfromTopicExpected+=_model->thetaEx[d][k]*_model->phiEx[k][v]*_model->z[d].size();
            }
            numerator+=(boost::math::digamma(numOfVocasfromTopicExpected+_model->beta[v])-boost::math::digamma(_model->beta[v]))*_model->beta[v];
            numOfTopicsExpected+=numOfVocasfromTopicExpected;
        }
        denominator+=boost::math::digamma(numOfTopicsExpected+_model->betaSum)-boost::math::digamma(_model->betaSum);
    }
    betaUpdated=numerator/denominator/_model->V;
    _model->beta.assign(_model->V,betaUpdated);
    _model->calculateHyperParamSum();
}//}}}
template<typename MODEL>
void BaseGibbsSampler<MODEL>::printTheta()const{//{{{
    cout<<fixed<<setprecision(5);
    cout<<"theta:"<<endl;
    for(int i=0;i<_model->theta.size();i++){
        for(int j=0;j<_model->theta[i].size();j++){
            cout<<_model->theta[i][j]<<' ';
        }
        cout<<endl;
    }
}//}}}
template<typename MODEL>
void BaseGibbsSampler<MODEL>::printPhi()const{//{{{
    cout<<fixed<<setprecision(5);
    cout<<"phi:"<<endl;
    for(int i=0;i<_model->phi.size();i++){
        for(int j=0;j<_model->phi[i].size();j++){
            cout<<_model->phi[i][j]<<' ';
        }
        cout<<endl;
    }
}//}}}
template<typename MODEL>
void BaseGibbsSampler<MODEL>::printThetaEx()const{//{{{
    cout<<fixed<<setprecision(5);
    cout<<"thetaEx:"<<endl;
    for(int i=0;i<_model->thetaEx.size();i++){
        for(int j=0;j<_model->thetaEx[i].size();j++){
            cout<<_model->thetaEx[i][j]<<' ';
        }
        cout<<endl;
    }
}//}}}
template<typename MODEL>
void BaseGibbsSampler<MODEL>::printPhiEx()const{//{{{
    cout<<fixed<<setprecision(5);
    cout<<"phiEx:"<<endl;
    for(int i=0;i<_model->phiEx.size();i++){
        for(int j=0;j<_model->phiEx[i].size();j++){
            cout<<_model->phiEx[i][j]<<' ';
        }
        cout<<endl;
    }
}//}}}
template<typename MODEL>
void BaseGibbsSampler<MODEL>::printNum()const{//{{{
    cout<<fixed<<setprecision(5);
    cout<<"numOfOccurencesOfTopic:"<<endl;
    for(int i=0;i<_model->topicCounts.size();i++){
        cout<<_model->topicCounts[i]<<' ';
    }
    cout<<endl;

    cout<<"numOfOccurencesOfTopicInDoc:"<<endl;
    for(int i=0;i<_model->topicCountsInDoc.size();i++){
        for(int j=0;j<_model->topicCountsInDoc[i].size();j++){
            cout<<_model->topicCountsInDoc[i][j]<<' ';
        }
        cout<<endl;
    }

    cout<<"numOfOccurencesOfVocaFromTopic:"<<endl;
    for(int i=0;i<_model->vocaCountsFromTopic.size();i++){
        for(int j=0;j<_model->vocaCountsFromTopic[i].size();j++){
            cout<<_model->vocaCountsFromTopic[i][j]<<' ';
        }
        cout<<endl;
    }
}//}}}
template<typename MODEL>
void BaseGibbsSampler<MODEL>::printHyperParameter()const{//{{{
    cout<<fixed<<setprecision(5);
    cout<<"alpha:"<<endl;
    for(int i=0;i<_model->alpha.size();i++){
        cout<<_model->alpha[i]<<' ';
    }
    cout<<endl;
    cout<<"beta:"<<endl;
    for(int i=0;i<_model->beta.size();i++){
        cout<<_model->beta[i]<<' ';
    }
    cout<<endl;
}//}}}
template<typename MODEL>
void BaseGibbsSampler<MODEL>::printTopFactor(int numOfTopFactor)const{//{{{
    vector<vector<double> > topFactorOfTheta;
    vector<vector<int> >  topFactorOfThetaIndex;
    vector<vector<double> > topFactorOfPhi;
    vector<vector<int> >  topFactorOfPhiIndex;
    vector<double> maxValue;
    vector<int> maxValueIndex;

    for(int d=0;d<_thetaAve.size();d++){
        for(int k=0;k<_thetaAve[d].size();k++){
            if(k<numOfTopFactor){
                maxValue.push_back(_thetaAve[d][k]);
                maxValueIndex.push_back(k);
            }else{
                double minInMaxValue=maxValue[0];
                int minInMaxValueIndex=0;
                for(int i=1;i<maxValue.size();i++){//find min in maxValue
                    if(minInMaxValue>maxValue[i]){
                        minInMaxValue=maxValue[i];
                        minInMaxValueIndex=i;
                    }
                }
                if(minInMaxValue<_thetaAve[d][k]){
                    // cout<<"min:"<<minInMaxValue<<"_thetaAve:"<<_thetaAve[d][k]<<endl;
                    maxValue[minInMaxValueIndex]=_thetaAve[d][k];
                    maxValueIndex[minInMaxValueIndex]=k;
                }
            }
        }
        topFactorOfTheta.push_back(maxValue);
        topFactorOfThetaIndex.push_back(maxValueIndex);
        maxValue.clear();
        maxValueIndex.clear();
    }
    for(int k=0;k<_phiAve.size();k++){
        for(int v=0;v<_phiAve[k].size();v++){
            if(v<numOfTopFactor){
                maxValue.push_back(_phiAve[k][v]);
                maxValueIndex.push_back(v);
            }else{
                double minInMaxValue=maxValue[0];
                int minInMaxValueIndex=0;
                for(int i=1;i<maxValue.size();i++){//find min in maxValue
                    if(minInMaxValue>maxValue[i]){
                        minInMaxValue=maxValue[i];
                        minInMaxValueIndex=i;
                    }
                }
                if(minInMaxValue<_phiAve[k][v]){
                    // cout<<"min:"<<minInMaxValue<<"_phiAve:"<<_phiAve[k][v]<<endl;
                    maxValue[minInMaxValueIndex]=_phiAve[k][v];
                    maxValueIndex[minInMaxValueIndex]=v;
                }
            }
        }
        topFactorOfPhi.push_back(maxValue);
        topFactorOfPhiIndex.push_back(maxValueIndex);
        maxValue.clear();
        maxValueIndex.clear();
    }
    for(int d=0;d<topFactorOfTheta.size();d++){
        cout<<"doc"<<d<<':'<<endl;
        for(int n=0;n<topFactorOfTheta[d].size();n++){
            cout<<"topic"<<topFactorOfThetaIndex[d][n]<<": "<<topFactorOfTheta[d][n]<<endl;
        }
    }
    for(int k=0;k<topFactorOfPhi.size();k++){
        cout<<"topic"<<k<<':'<<endl;
        for(int n=0;n<topFactorOfPhi[k].size();n++){
            cout<<_wordList[topFactorOfPhiIndex[k][n]]<<": "<<topFactorOfPhi[k][n]<<endl;
        }
    }
    cout<<"alpha"<<endl;
    for(int i=0;i<_model->alpha.size();i++){
        cout<<_model->alpha[i]<<',';
    }
    cout<<endl;

    cout<<"beta"<<endl;
    cout<<_model->beta[0];
    cout<<endl;
}//}}}
template<typename MODEL>
void BaseGibbsSampler<MODEL>::computeParameter(){//{{{
    for(int d=0;d<_thetaAve.size();d++){
        for(int k=0;k<_model->K;k++){
            _thetaAve[d][k]=_thetaAve[d][k]/static_cast<double>(_nThIterationSample.size());
        }
    }
    for(int k=0;k<_model->K;k++){
        for(int v=0;v<_model->V;v++){
            _phiAve[k][v]=_phiAve[k][v]/static_cast<double>(_nThIterationSample.size());
        }
    }

}//}}}
template<typename MODEL>
void BaseGibbsSampler<MODEL>::writeParameter(string thetaFilename,string phiFilename,string alphaFilename,string betaFilename)const{//{{{
    ofstream thetaOutput;
    ofstream phiOutput;
    ofstream alphaOutput;
    ofstream betaOutput;
    thetaOutput.open(thetaFilename,ios::out);
    phiOutput.open(phiFilename,ios::out);
    alphaOutput.open(alphaFilename,ios::out);
    betaOutput.open(betaFilename,ios::out);

    for(int i=0;i<_thetaAve.size();i++){
        for(int j=0;j<_thetaAve[i].size();j++){
            thetaOutput<<_thetaAve[i][j];
            if(j!=(_thetaAve[i].size()-1)){
                thetaOutput<<',';
            }
        }
        thetaOutput<<endl;
    }

    for(int i=0;i<_phiAve.size();i++){
        for(int j=0;j<_phiAve[i].size();j++){
            phiOutput<<_phiAve[i][j];
            if(j!=(_phiAve[i].size()-1)){
                phiOutput<<',';
            }
        }
        phiOutput<<endl;
    }
    for(int i=0;i<_model->alpha.size();i++){
        alphaOutput<<_model->alpha[i];
        alphaOutput<<endl;
    }
    for(int i=0;i<_model->beta.size();i++){
        betaOutput<<_model->beta[i];
        betaOutput<<endl;
    }
    thetaOutput.close();
    phiOutput.close();
    alphaOutput.close();
    betaOutput.close();
}//}}}
template<typename MODEL>
void BaseGibbsSampler<MODEL>::writeTheta(string thetaFilename)const{//{{{
    ofstream thetaOutput;
    thetaOutput.open(thetaFilename,ios::out);

    for(int i=0;i<_model->theta.size();i++){
        for(int j=0;j<_model->theta[i].size();j++){
            thetaOutput<<_model->theta[i][j];
            if(j!=(_model->theta[i].size()-1)){
                thetaOutput<<',';
            }
        }
        thetaOutput<<endl;
    }
}//}}}
template<typename MODEL>
void BaseGibbsSampler<MODEL>::writePhi(string phiFilename)const{//{{{
    ofstream phiOutput;
    phiOutput.open(phiFilename,ios::out);

    for(int i=0;i<_model->phi.size();i++){
        for(int j=0;j<_model->phi[i].size();j++){
            phiOutput<<_model->phi[i][j];
            if(j!=(_model->phi[i].size()-1)){
                phiOutput<<',';
            }
        }
        phiOutput<<endl;
    }
}//}}}
template<typename MODEL>
void BaseGibbsSampler<MODEL>::samplingParameter(){//{{{

    for(int d=0;d<_thetaAve.size();d++){
        for(int k=0;k<_model->K;k++){
            _thetaAve[d][k]+=_model->thetaEx[d][k];
        }
    }
    for(int k=0;k<_model->K;k++){
        for(int v=0;v<_model->V;v++){
            _phiAve[k][v]+=_model->phiEx[k][v];
        }
    }
}//}}}

typedef BaseGibbsSampler<LDA> GibbsSamplerFromLDA;

#endif
