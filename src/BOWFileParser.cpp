//
// BOWFileParser.hpp
//
// Copyright (c) 2017 Shion Hosoda
//
// This software is released under the MIT License.
// http://opensource.org/licenses/mit-license.php
//

#include"include/BOWFileParser.hpp"
#include<sstream>

using namespace std;

BOWFileParser::BOWFileParser(string filename):_inputText(filename),_V(0){//{{{

}//}}}

BOWFileParser::~BOWFileParser(){//{{{
    _inputText.close();
}//}}}

void BOWFileParser::readBOWFile(){//{{{
    if(!_inputText){
        cout<<"Can't open BOWfile";
        exit(1);
    }
    string str;
    for(int i=0;getline(_inputText,str);i++){
        vector<unsigned int > intVec;
        string token;
        istringstream stream(str);

        if(i==0){
            for(int j=0;getline(stream,token,',');j++){
                if(j==0)continue;
                _wordList.push_back(token);
            }
            continue;
        }
        for(int j=0;getline(stream,token,',');j++){
            if(j==0)continue;
            intVec.push_back(stoi(token));
        }
        _frequencyMatrix.push_back(intVec);
    }


    _V=_wordList.size();
}//}}}

void BOWFileParser::makeBagOfWordsNum(){//{{{
    for(int i=0;i<_frequencyMatrix.size();i++){
        vector<unsigned int> intVec;
        for(int j=0;j<_frequencyMatrix[i].size();j++){
            for(int k=0;k<_frequencyMatrix[i][j];k++){
                intVec.push_back(j);
            }
        }
        _bagOfWordsNum.push_back(intVec);
    }
}//}}}

void BOWFileParser::writeWordList(string wordListFilename)const{//{{{
    ofstream wordListFp;
    wordListFp.open(wordListFilename,ios::out);
    for(int v=0;v<_wordList.size();v++){
        wordListFp<<_wordList[v];
        if(v!=(_wordList.size()-1))wordListFp<<endl;
    }
    wordListFp.close();

}//}}}
