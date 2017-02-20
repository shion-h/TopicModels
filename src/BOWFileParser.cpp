//
// BOWFileParser.hpp
//
// Copyright (c) 2017 Shion Hosoda
//
// This software is released under the MIT License.
// http://opensource.org/licenses/mit-license.php
//

#include"include/BOWFileParser.hpp"

using namespace std;

BOWFileParser::BOWFileParser(string filename):_inputText(filename),_V(0){//{{{

}//}}}

BOWFileParser::~BOWFileParser(){//{{{
    _inputText.close();
}//}}}

void BOWFileParser::readBOWFile(){//{{{
    vector<string > vbuf;
    string buf;

    if(!_inputText){
        cout<<"erRor";
        exit(1);
    }
    char c;
    while(_inputText.get(c)){
        if(c==','){
            // buf+='\0';
            vbuf.push_back(buf);
            buf.clear();
        }else if(c=='\n'){
            // buf+='\0';
            vbuf.push_back(buf);
            buf.clear();
            _bagOfWords.push_back(vbuf);
            vbuf.clear();
        }else{
            buf+=c;
        }
    }

}//}}}

void BOWFileParser::makeWordList(){//{{{
    unsigned int flag=0;

    for(int d=0;d<_bagOfWords.size();d++){
        for(int i=0;i<_bagOfWords[d].size();i++){
            for(int v=0;v<_V;v++){
                if(_wordList[v]==_bagOfWords[d][i]){
                    flag=1;
                    break;
                }
            }
            if(flag==0){//new word
                _wordList.push_back(_bagOfWords[d][i]);
                _V++;
            }
            flag=0;
        }
    }

}//}}}

void BOWFileParser::writeWordList(string wordListFilename)const{//{{{
    ofstream _wordListFp;
    _wordListFp.open(wordListFilename,ios::out);
    for(int v=0;v<_V;v++){
        _wordListFp<<_wordList[v];
        if(v!=(_V-1))_wordListFp<<endl;
    }
    _wordListFp.close();

}//}}}

void BOWFileParser::assignNumbersToWords(){//{{{
    for(int d=0;d<_bagOfWords.size();d++){
        vector<unsigned int> ibuf;
        for(int i=0;i<_bagOfWords[d].size();i++){
            for(int v=0;v<_V;v++){
                if(_bagOfWords[d][i]==_wordList[v]){
                    ibuf.push_back(v);
                    break;
                }
            }
        }
        _bagOfWordsNum.push_back(ibuf);
    }

}//}}}
