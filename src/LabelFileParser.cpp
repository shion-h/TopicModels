//
// LabelFileParser.hpp
//
// Copyright (c) 2017 Shion Hosoda
//
// This software is released under the MIT License.
// http://opensource.org/licenses/mit-license.php
//

#include"include/LabelFileParser.hpp"
#include<sstream>

using namespace std;

LabelFileParser::LabelFileParser(string filename):_inputText(filename),_C(0){//{{{
    this->readLabelFile();
    this->makeLabelList();
}//}}}

LabelFileParser::~LabelFileParser(){//{{{
    _inputText.close();
}//}}}

void LabelFileParser::readLabelFile(){//{{{
    if(!_inputText){
        cout<<"Can't open labelfile";
        exit(1);
    }
    string str;
    for(int i=0;getline(_inputText,str);i++){
        string token;
        istringstream stream(str);
        if(i == 0) continue;
        for(int j=0; getline(stream,token,','); j++){
            if(j == 0)continue;
            if(j > 1){
                cout<<"Label file must be one dimension.";
                exit(1);
            }
            _label.push_back(token);
        }
    }
}//}}}

void LabelFileParser::makeLabelList(){//{{{
    auto p = _label.begin();
    unsigned int count = 0;
    for(; p!=_label.end(); p++){
        vector<string>::iterator sIter = find(_label.begin(), p, *p);
        if(sIter != p){
            vector<string>::iterator lIter = find(_labelList.begin(), _labelList.end(), *sIter);
            unsigned int index = distance(_labelList.begin(), lIter);
            _y.push_back(index);
        }else{
            _y.push_back(count);
            _labelList.push_back(*p);
            count++;
        }
    }
    if(_labelList.size() < 2){
        cout<<"The number of different labels must be two or more."<<endl<<flush;
        exit(1);
    }
    _C = _labelList.size();
}//}}}

void LabelFileParser::writeLabelList(string labelListFilename)const{//{{{
    ofstream labelListFp;
    labelListFp.open(labelListFilename,ios::out);
    for(int c=0; c<_labelList.size(); c++){
        labelListFp<<_labelList[c];
        if(c != (_labelList.size()-1))labelListFp<<endl;
    }
    labelListFp.close();
}//}}}
