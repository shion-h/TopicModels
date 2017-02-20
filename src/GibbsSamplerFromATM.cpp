//
// GibbsSamplerFromATM.cpp
//
// Copyright (c) 2017 Shion Hosoda
//
// This software is released under the MIT License.
// http://opensource.org/licenses/mit-license.php
//

#include"include/GibbsSamplerFromATM.hpp"

using namespace std;

ATM::ATM(const vector<vector<unsigned int > > &bagOfWordsNum,const vector<vector<unsigned int> > &bagOfAuthorsNum,const unsigned int K,const unsigned int V,const unsigned int A)//{{{
    :LDA(bagOfWordsNum,K,V),
     bagOfAuthorsNum(bagOfAuthorsNum),
     A(A),
     y(bagOfWordsNum.size()),
     wordCountsInAuthor(A),
     topicCountsInAuthor(A){
    initializeParam();
}//}}}

void ATM::initializeParam(){//{{{
    for(int d=0;d<y.size();d++)y[d].assign(bagOfWordsNum[d].size(),0);
    theta.assign(A,*(new vector<double>));
    for(int a=0;a<theta.size();a++)theta[a].assign(K,0);
    thetaEx.assign(A,*(new vector<double>));
    for(int a=0;a<thetaEx.size();a++)thetaEx[a].assign(K,0);
    wordCountsInAuthor.assign(A,0);
    for(int a=0;a<topicCountsInAuthor.size();a++)topicCountsInAuthor[a].assign(K,0);
}//}}}

void ATM::calculateEx(){//{{{
    for(int a=0;a<A;a++){
        for(int k=0;k<K;k++){
            thetaEx[a][k]=(topicCountsInAuthor[a][k]+alpha[k])/(wordCountsInAuthor[a]+alphaSum);
        }
    }
    for(int k=0;k<K;k++){
        for(int v=0;v<V;v++){
            phiEx[k][v]=(vocaCountsFromTopic[k][v]+beta[v])/(topicCounts[k]+betaSum);
        }
    }

}//}}}

void ATM::nUpdate(){//{{{
    for(int d=0;d<topicCountsInAuthor.size();d++)topicCountsInAuthor[d].assign(K,0);
    for(int k=0;k<vocaCountsFromTopic.size();k++)vocaCountsFromTopic[k].assign(V,0);
    topicCounts.assign(K,0);
    wordCountsInAuthor.assign(A,0);

    for(int d=0;d<z.size();d++){
        for(int i=0;i<z[d].size();i++){
            topicCountsInAuthor[y[d][i]][z[d][i]]++;
        }
    }
    for(int d=0;d<z.size();d++){
        for(int i=0;i<z[d].size();i++){
                vocaCountsFromTopic[z[d][i]][bagOfWordsNum[d][i]]++;
        }
    }
    for(int k=0;k<K;k++){
        for(int v=0;v<V;v++){
            topicCounts[k]+=vocaCountsFromTopic[k][v];
        }
    }
    for(int d=0;d<z.size();d++){
        for(int i=0;i<z[d].size();i++){
            wordCountsInAuthor[y[d][i]]++;
        }
    }

}//}}}

GibbsSamplerFromATM::GibbsSamplerFromATM(const vector<vector<unsigned int > > &bagOfWordsNum,const vector<vector<unsigned int> > &bagOfAuthorsNum,const vector<string> &wordList,const vector<string> &authorList,const unsigned int K,const unsigned int V,unsigned int A,vector<unsigned int> nThIterationSample)//{{{
    :BaseGibbsSampler(wordList,nThIterationSample),
     _authorList(authorList){
    _thetaAve.assign(A,*(new vector<double>));
    _phiAve.assign(K,*(new vector<double>));
    _model=make_shared<ATM>(bagOfWordsNum,bagOfAuthorsNum,K,V,A);
    for(int d=0;d<_thetaAve.size();d++)_thetaAve[d].assign(_model->K,0);
    for(int k=0;k<_phiAve.size();k++)_phiAve[k].assign(_model->V,0);
}//}}}

GibbsSamplerFromATM::~GibbsSamplerFromATM(){//{{{

}//}}}

void GibbsSamplerFromATM::yUpdate(){//{{{
    random_device rnd;
    mt19937 mt(rnd());
    for(int d=0;d<_model->y.size();d++){
        uniform_int_distribution<int> randN(0,_model->bagOfAuthorsNum[d].size()-1);
        for(int i=0;i<_model->y[d].size();i++){
            _model->y[d][i]=_model->bagOfAuthorsNum[d][randN(mt)];
        }
    }
}//}}}

void GibbsSamplerFromATM::hyperParamUpdate(){//{{{
    _model->calculateEx();
    double numerator=0,denominator=0;
    for(int k=0;k<_model->K;k++){
        for(int a=0;a<_model->A;a++){
            double wordCountsInAuthorExpected=0;
            for(int d=0;d<_model->z.size();d++){
                if(find(_model->bagOfAuthorsNum[d].begin(),_model->bagOfAuthorsNum[d].begin(),a)!=_model->bagOfAuthorsNum[d].end()){
                    wordCountsInAuthorExpected+=1.0/_model->bagOfAuthorsNum[d].size()*_model->z[d].size();
                }
            }
            numerator+=(boost::math::digamma(_model->thetaEx[a][k]*wordCountsInAuthorExpected+_model->alpha[k])-boost::math::digamma(_model->alpha[k]))*_model->alpha[k];
            denominator+=boost::math::digamma(wordCountsInAuthorExpected+_model->alphaSum)-boost::math::digamma(_model->alphaSum);
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
            for(int a=0;a<_model->A;a++){
                double wordCountsInAuthorExpected=0;
                for(int d=0;d<_model->z.size();d++){
                    if(find(_model->bagOfAuthorsNum[d].begin(),_model->bagOfAuthorsNum[d].begin(),a)!=_model->bagOfAuthorsNum[d].end()){
                        wordCountsInAuthorExpected+=1.0/_model->bagOfAuthorsNum[d].size()*_model->z[d].size();
                    }
                }
                numOfVocasfromTopicExpected+=_model->thetaEx[a][k]*_model->phiEx[k][v]*wordCountsInAuthorExpected;
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

void GibbsSamplerFromATM::printTopFactor(int numOfTopFactor)const{//{{{
    vector<vector<double> > topFactorOfTheta;
    vector<vector<int> >  topFactorOfThetaIndex;
    vector<vector<double> > topFactorOfPhi;
    vector<vector<int> >  topFactorOfPhiIndex;
    vector<double> maxValue;
    vector<int> maxValueIndex;

    for(int a=0;a<_thetaAve.size();a++){
        for(int k=0;k<_thetaAve[a].size();k++){
            if(k<numOfTopFactor){
                maxValue.push_back(_thetaAve[a][k]);
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
                if(minInMaxValue<_thetaAve[a][k]){
                    maxValue[minInMaxValueIndex]=_thetaAve[a][k];
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
    for(int a=0;a<topFactorOfTheta.size();a++){
        cout<<_authorList[a]<<':'<<endl;
        for(int n=0;n<topFactorOfTheta[a].size();n++){
            cout<<"topic"<<topFactorOfThetaIndex[a][n]<<": "<<topFactorOfTheta[a][n]<<endl;
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
