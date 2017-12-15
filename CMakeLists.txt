cmake_minimum_required(VERSION 3.5)
set (CMAKE_CXX_STANDARD 11)
# INCLUDE_DIRECTORIES(/usr/include/python2.7)
add_executable(LDA src/LDA.cpp src/VariationalBayesEstimatorOnLDA.cpp src/BOWFileParser.cpp)
add_executable(LDAWBIC src/LDAWBIC.cpp src/LDAWBICMetropolisSampler.cpp
    src/BOWFileParser.cpp)
TARGET_LINK_LIBRARIES(LDA boost_program_options)
