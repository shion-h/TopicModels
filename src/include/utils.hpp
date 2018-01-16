//
// utils.hpp
//
// Copyright (c) 2018 Shion Hosoda
//
// This software is released under the MIT License.
// http://opensource.org/licenses/mit-license.php
//


#include<vector>
#include<string>
#include<iostream>
#include<fstream>


template<typename T>
void outputVector(const std::vector<std::vector<T> > &vector, std::string filename){
    std::ofstream stream;
    stream.open(filename, std::ios::out);

    for(int i=0;i<vector.size();i++){
        for(int j=0;j<vector[i].size();j++){
            stream<<vector[i][j];
            if(j!=(vector[i].size()-1)){
                stream<<',';
            }
        }
        stream<<std::endl;
    }
    stream.close();
}

template<typename T>
void outputVector(const std::vector<T> &vector, std::string filename){
    std::ofstream stream;
    stream.open(filename, std::ios::out);
    for(int i=0;i<vector.size();i++){
        stream<<vector[i];
        stream<<std::endl;
    }
    stream.close();
}
