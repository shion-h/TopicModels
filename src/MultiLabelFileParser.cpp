//
// MultiLabelFileParser.cpp
//
// Copyright (c) 2017 Shion Hosoda
//
// This software is released under the MIT License.
// http://opensource.org/licenses/mit-license.php
//

#include"include/MultiLabelFileParser.hpp"
#include<sstream>

using namespace std;

MultiLabelFileParser::MultiLabelFileParser(string filename):_inputText(filename){//{{{
    this->readLabelFile();
    this->makeLabelList();
}//}}}

MultiLabelFileParser::~MultiLabelFileParser(){//{{{
    _inputText.close();
}//}}}

void MultiLabelFileParser::readLabelFile(){//{{{
    if(!_inputText){
        cout<<"Can't open labelfile";
        exit(1);
    }
    string line;
    vector<vector<string> > labelMatrixBuf;
    for(int i=0;getline(_inputText, line);i++){
        Tokenizer tok(line);
        vector<string> labelBuf;
        if(i == 0){
            _labelColumnList.assign(tok.begin(),tok.end());
            _labelColumnList.erase(_labelColumnList.begin(), _labelColumnList.begin());
        }else{
            labelBuf.assign(tok.begin(),tok.end());
            labelBuf.erase(labelBuf.begin(), labelBuf.begin()+1);
            labelMatrixBuf.push_back(labelBuf);
        }
    }
    _labelMatrix = vector<vector<string> >(labelMatrixBuf[0].size(), vector<string>(labelMatrixBuf.size(), ""));
    for(int i=0; i<_labelMatrix.size(); i++){
        for(int j=0; j<_labelMatrix[0].size(); j++){
            _labelMatrix[i][j] = labelMatrixBuf[j][i];
        }
    }
    _L = _labelMatrix.size();
}//}}}

void MultiLabelFileParser::makeLabelList(){//{{{
    for(int l=0; l<_labelMatrix.size(); l++){
        auto p = _labelMatrix[l].begin();
        unsigned int count = 0;
        vector<unsigned int> yl;
        vector<string> labelListl;
        for(; p!=_labelMatrix[l].end(); p++){
            vector<string>::iterator sIter = find(_labelMatrix[l].begin(), p, *p);
            if(sIter != p){
                vector<string>::iterator lIter = find(labelListl.begin(), labelListl.end(), *sIter);
                unsigned int index = distance(labelListl.begin(), lIter);
                yl.push_back(index);
            }else{
                yl.push_back(count);
                labelListl.push_back(*p);
                count++;
            }
        }
        _y.push_back(yl);
        _labelList.push_back(labelListl);
        _C.push_back(labelListl.size());
    }
    cout<<_L<<endl;
    outputVectorForDebug(_labelList[0], "l0");
    outputVectorForDebug(_labelList[1], "l1");
}//}}}

void MultiLabelFileParser::writeLabelList(string labelListFilename)const{//{{{
    ofstream labelListFp;
    labelListFp.open(labelListFilename,ios::out);
    for(int l=0; l<_L; l++){
        for(int c=0; c<_C[l]; c++){
            labelListFp<<_labelList[l][c];
            if(c != (_C[l]-1))labelListFp<<',';
        }
        if(l != (_L-1))labelListFp<<endl;
    }
    labelListFp.close();
}//}}}
