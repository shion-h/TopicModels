cmake_minimum_required(VERSION 3.5)
set (CMAKE_CXX_STANDARD 11)
add_executable(LDA src/LDA.cpp src/VariationalBayesEstimatorOnLDA.cpp src/BOWFileParser.cpp)
add_executable(sLDA src/sLDA.cpp src/VariationalBayesEstimatorOnSLDA.cpp src/VariationalBayesEstimatorOnLDA.cpp
    src/BOWFileParser.cpp src/LabelFileParser.cpp src/MultiLabelFileParser.cpp)
add_executable(ATM src/ATM.cpp src/VariationalBayesEstimatorOnATM.cpp src/BOWFileParser.cpp)
add_executable(HDP src/HDP.cpp src/GibbsSamplerFromHDP.cpp src/BOWFileParser.cpp)
add_executable(LDAWBIC src/LDAWBIC.cpp src/LDAWBICMetropolisSampler.cpp
    src/BOWFileParser.cpp)
TARGET_LINK_LIBRARIES(LDA boost_program_options)
TARGET_LINK_LIBRARIES(sLDA boost_program_options)
TARGET_LINK_LIBRARIES(ATM boost_program_options)
TARGET_LINK_LIBRARIES(HDP boost_program_options)
