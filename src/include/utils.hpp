//
// utils.hpp
//
// Copyright (c) 2018 Shion Hosoda
//
// This software is released under the MIT License.
// http://opensource.org/licenses/mit-license.php
//

#ifndef UTILS
#define UTILS

#include<stdlib.h>
#include<math.h>
#include<cmath>
#include<iostream>
#include<vector>
#include<numeric>
#include<memory>
#include<random>
#include<iomanip>
#include<fstream>
#include<limits>




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

template<typename T>
void outputVectorForDebug(const std::vector<T> &vector, std::string msg){
    std::ofstream stream;
    stream.open("./check", std::ios::app);
    stream<<msg<<std::endl;
    for(int i=0; i<vector.size(); i++){
        stream<<vector[i]<<',';
    }
    stream<<std::endl;
    stream.close();
}

template<typename T>
void printVector(const std::vector<std::vector<T> > &vector){
    for(int i=0;i<vector.size();i++){
        for(int j=0;j<vector[i].size();j++){
            std::cout<<vector[i][j];
            if(j!=(vector[i].size()-1)){
                std::cout<<',';
            }
        }
        std::cout<<std::endl;
    }
}

template<typename T>
void printVector(const std::vector<T> &vector){
    for(int i=0;i<vector.size();i++){
        std::cout<<vector[i]<<',';
    }
    std::cout<<std::endl;
}

template<typename T>
unsigned int sampleDiscreteValues(std::vector<T> &discreteDistribution, bool isLogValue = false){
    if(isLogValue == true){
        std::vector<T> discreteNotLogDistribution(discreteDistribution);
        double bottomValue = *std::max_element(discreteDistribution.begin(), discreteDistribution.end());
        for(int i=0; i<discreteNotLogDistribution.size(); i++){
            discreteNotLogDistribution[i] -= bottomValue;
            discreteNotLogDistribution[i] = exp(discreteNotLogDistribution[i]);
            if(std::isinf(discreteNotLogDistribution[i])){
                discreteNotLogDistribution[i] = std::numeric_limits<T>::max();
            }
        }
        return sampleDiscreteValues(discreteNotLogDistribution, false);
    }
    double normalizationConstant = accumulate(discreteDistribution.begin(), discreteDistribution.end(), 0.0);
    std::random_device rnd;
    std::mt19937 mt(rnd());
    std::uniform_real_distribution<> rand12(0.0, normalizationConstant);
    double probability = 0;
    double randomValue = rand12(mt);
    unsigned int i;
    for(i=0; i<discreteDistribution.size(); i++){
        probability += discreteDistribution[i];
        if(probability > randomValue){
            break;
        }
    }
    return i;
}

template<typename T>
void fillTensor(std::vector<std::vector<std::vector<T> > > &tensor, T value){
    for(int i=0; i<tensor.size(); i++){
        for(int j=0; j<tensor[i].size(); j++){
            for(int k=0; k<tensor[i][j].size(); k++){
                tensor[i][j][k] = value;
            }
        }
    }
}

template<typename T>
void fillMatrix(std::vector<std::vector<T> > &matrix, T value){
    for(int i=0; i<matrix.size(); i++){
        for(int j=0; j<matrix[i].size(); j++){
            matrix[i][j] = value;
        }
    }
}

template<typename T>
T sum(std::vector<std::vector<std::vector<T> > > tensor){
    T result(0);
    for(int i=0; i<tensor.size(); i++){
        for(int j=0; j<tensor[i].size(); j++){
            result += accumulate(tensor[i][j].begin(), tensor[i][j].end(), 0.0);
        }
    }
    return result;
}

template<typename T>
T sum(std::vector<std::vector<T> > matrix){
    T result(0);
    for(int i=0; i<matrix.size(); i++){
        result += accumulate(matrix[i].begin(), matrix[i].end(), 0.0);
    }
    return result;
}

template<typename T>
T sum(std::vector<T> vector){
    return accumulate(vector.begin(), vector.end(), 0.0);
}

template<typename T>
std::vector<std::vector<T> > prod(std::vector<std::vector<T> > matrix1, std::vector<std::vector<T> > matrix2){
    std::vector<std::vector<T> > resultMatrix(matrix1.size(), std::vector<T>(matrix2[0].size(), 0));
    for(int i=0; i<matrix1.size(); i++){
        for(int j=0; j<matrix1[0].size(); j++){
            for(int k=0; k<matrix2[0].size(); k++){
                resultMatrix[i][k] += matrix1[i][j] * matrix2[j][k];
            }
        }
    }
    return resultMatrix;
}
#endif
