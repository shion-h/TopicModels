//
// BOWFileParser.hpp
//
// Copyright (c) 2017 Shion Hosoda
//
// This software is released under the MIT License.
// http://opensource.org/licenses/mit-license.php
//

#ifndef BOWFILE
#define BOWFILE
#include<iostream>
#include<string>
#include<vector>
#include<fstream>

using namespace std;

class BOWFileParser{
protected:
    ifstream _inputText;
    vector<string> _wordList;
    unsigned int _V;
    vector<vector<unsigned int> > _frequencyMatrix;
    vector<vector<unsigned int> > _bagOfWordsNum;
public:
    BOWFileParser(string filename);
    ~BOWFileParser();
    void readBOWFile();
    void makeBagOfWordsNum();
    void writeWordList(string wordListFilename)const;
    const vector<vector<unsigned int> > &getBagOfWordsNum()const{
        return _bagOfWordsNum;
    }
    const vector<string> &getWordList()const{
        return _wordList;
    }
    const unsigned int &getV()const{
        return _V;
    }
};

#endif
