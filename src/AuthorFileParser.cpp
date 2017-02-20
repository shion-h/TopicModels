//
// AuthorFileParser.cpp
//
// Copyright (c) 2017 Shion Hosoda
//
// This software is released under the MIT License.
// http://opensource.org/licenses/mit-license.php
//

#include"include/AuthorFileParser.hpp"

using namespace std;

AuthorFileParser::AuthorFileParser(string filename):BOWFileParser(filename){//{{{

}//}}}

void AuthorFileParser::readAuthorFile(){//{{{
    readBOWFile();
}//}}}

void AuthorFileParser::makeAuthorList(){//{{{
    makeWordList();
}//}}}

void AuthorFileParser::writeAuthorList(string authorListFilename)const{//{{{
    writeWordList(authorListFilename);
}//}}}

void AuthorFileParser::assignNumbersToAuthors(){//{{{
    assignNumbersToWords();
}//}}}
