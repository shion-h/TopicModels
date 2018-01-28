//
// BOWFileParser.hpp
//
// Copyright (c) 2017 Shion Hosoda
//
// This software is released under the MIT License.
// http://opensource.org/licenses/mit-license.php
//

#ifndef BOWFILEPASER
#define BOWFILEPASER
#include<iostream>
#include<string>
#include<vector>
#include<fstream>

class BOWFileParser{
protected:
    std::ifstream _inputText;
    std::vector<std::string> _wordList;
    unsigned int _V;
    std::vector<std::vector<unsigned int> > _frequencyMatrix;
    std::vector<std::vector<unsigned int> > _docVoca;
public:
    BOWFileParser(std::string filename);
    ~BOWFileParser();
    void readBOWFile();
    void makeDocVoca();
    void writeWordList(std::string wordListFilename)const;
    const std::vector<std::vector<unsigned int> > &getFrequencyMatrix()const{
        return _frequencyMatrix;
    }
    const std::vector<std::vector<unsigned int> > &getDocVoca()const{
        return _docVoca;
    }
    const std::vector<std::string> &getWordList()const{
        return _wordList;
    }
    const unsigned int &getV()const{
        return _V;
    }
};

#endif
