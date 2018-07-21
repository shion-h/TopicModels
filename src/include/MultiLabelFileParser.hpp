//
// MultiLabelFileParser.hpp
//
// Copyright (c) 2017 Shion Hosoda
//
// This software is released under the MIT License.
// http://opensource.org/licenses/mit-license.php
//
//
#ifndef MULTILABELFILEPASER
#define MULTILABELFILEPASER

#include<iostream>
#include<string>
#include<vector>
#include<algorithm>
#include<fstream>
#include <boost/tokenizer.hpp>
#include"utils.hpp"

typedef boost::tokenizer<boost::escaped_list_separator<char> > Tokenizer;
class MultiLabelFileParser{
protected:
    std::ifstream _inputText;
    std::vector<std::string> _labelColumnList;
    unsigned int _L;
    std::vector<std::vector<std::string> > _labelList;
    std::vector<unsigned int> _C;
    std::vector<std::vector<unsigned int> > _y;
    std::vector<std::vector<std::string> > _labelMatrix;
public:
    MultiLabelFileParser(std::string filename);
    ~MultiLabelFileParser();
    void readLabelFile();
    void makeLabelList();
    void writeLabelList(std::string labelListFilename)const;
    const std::vector<std::vector<unsigned int> > &getY()const{
        return _y;
    }
    const std::vector<std::string> &getLabelColumnList()const{
        return _labelColumnList;
    }
    const std::vector<std::vector<std::string> > &getLabelList()const{
        return _labelList;
    }
    const std::vector<unsigned int> &getC()const{
        return _C;
    }
    const unsigned int &getL()const{
        return _L;
    }
};
#endif
