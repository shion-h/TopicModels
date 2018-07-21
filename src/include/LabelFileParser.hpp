//
// LabelFileParser.hpp
//
// Copyright (c) 2017 Shion Hosoda
//
// This software is released under the MIT License.
// http://opensource.org/licenses/mit-license.php
//
//
#ifndef LABELFILEPASER
#define LABELFILEPASER

#include<iostream>
#include<string>
#include<vector>
#include<algorithm>
#include<fstream>

class LabelFileParser{
protected:
    std::ifstream _inputText;
    std::vector<std::string> _labelList;
    unsigned int _C;
    unsigned int _L;
    std::vector<unsigned int> _y;
    std::vector<std::string> _label;
public:
    LabelFileParser(std::string filename);
    ~LabelFileParser();
    void readLabelFile();
    void makeLabelList();
    void writeLabelList(std::string labelListFilename)const;
    const std::vector<unsigned int> &getY()const{
        return _y;
    }
    const std::vector<std::string> &getLabelList()const{
        return _labelList;
    }
    const unsigned int &getC()const{
        return _C;
    }
};
#endif
