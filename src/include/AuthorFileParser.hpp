//
// AuthorFileParser.hpp
//
// Copyright (c) 2017 Shion Hosoda
//
// This software is released under the MIT License.
// http://opensource.org/licenses/mit-license.php
//

#include"BOWFileParser.hpp"

using namespace std;

class AuthorFileParser:private BOWFileParser{
public:
    AuthorFileParser(string filename);
    void readAuthorFile();
    void makeAuthorList();
    void writeAuthorList(string authorListFilename)const;
    void assignNumbersToAuthors();
    const vector<vector<unsigned int> > &getBagOfAuthorsNum()const{
        return _bagOfWordsNum;
    }
    const vector<string> &getAuthorList()const{
        return _wordList;
    }
    const unsigned int &getY()const{
        return _V;
    }
};
