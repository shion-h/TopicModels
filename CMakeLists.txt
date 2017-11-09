cmake_minimum_required(VERSION 3.5)
set (CMAKE_CXX_STANDARD 11)
# INCLUDE_DIRECTORIES(/usr/include/python2.7)
add_executable(LDA src/LDA.cpp src/GibbsSamplerFromLDA.cpp
    src/CollapsedGibbsSamplerFromLDA.cpp src/OrdinaryGibbsSamplerFromLDA.cpp
    src/VariationalBayesEstimatorOnLDA.cpp src/BOWFileParser.cpp)
add_executable(ATM src/ATM.cpp src/GibbsSamplerFromLDA.cpp
    src/CollapsedGibbsSamplerFromLDA.cpp
    src/GibbsSamplerFromATM.cpp
    src/CollapsedGibbsSamplerFromATM.cpp src/BOWFileParser.cpp
    src/AuthorFileParser.cpp)
add_executable(LDAWBIC src/LDAWBIC.cpp src/LDAWBICMetropolisSampler.cpp
    src/BOWFileParser.cpp)
TARGET_LINK_LIBRARIES(LDA boost_program_options)
TARGET_LINK_LIBRARIES(ATM boost_program_options)
