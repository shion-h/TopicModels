//
// AuthorFileParser.hpp
//
// Copyright (c) 2017 Shion Hosoda
//
// This software is released under the MIT License.
// http://opensource.org/licenses/mit-license.php
//

#ifndef AUTHORFILEPASER
#define AUTHORFILEPASER
#include<iostream>
#include<string>
#include<vector>
#include<fstream>

class AuthorFileParser:private BOWFileParser{
public:
    AuthorFileParser(std::string filename):BOWFileParser(filename){};
    ~AuthorFileParser(){};
    const std::vector<std::vector<unsigned int> > &getDocAuth()const{
        return _docVoca;
    }
    const std::vector<std::string> &getAuthorList()const{
        return _wordList;
    }
    const unsigned int &getM()const{
        return _V;
    }
    void writeAuthorList(std::string authorListFilename)const{
        writeWordList(authorListFilename);
    };
};

#endif
