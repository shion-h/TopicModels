bin_PROGRAMS = LDA ATM LDAWBIC
LDA_SOURCES = BOWFileParser.cpp include/BOWFileParser.hpp CollapsedGibbsSamplerFromLDA.cpp include/CollapsedGibbsSamplerFromLDA.hpp GibbsSamplerFromLDA.cpp include/GibbsSamplerFromLDA.hpp LDA.cpp OrdinaryGibbsSamplerFromLDA.cpp include/OrdinaryGibbsSamplerFromLDA.hpp VariationalBayesEstimatorOnLDA.cpp include/VariationalBayesEstimatorOnLDA.hpp
LDA_CXXFLAGS = -std=c++11
LDA_LDFLAGS = -lboost_program_options
ATM_SOURCES = CollapsedGibbsSamplerFromATM.cpp include/CollapsedGibbsSamplerFromATM.hpp GibbsSamplerFromATM.cpp include/GibbsSamplerFromATM.hpp ATM.cpp AuthorFileParser.cpp include/AuthorFileParser.hpp
ATM_CXXFLAGS = -std=c++11
ATM_LDFLAGS = -lboost_program_options
ATM_LDADD = LDA-BOWFileParser.o LDA-GibbsSamplerFromLDA.o
LDAWBIC_SOURCES = LDAWBICMetropolisSampler.cpp include/LDAWBICMetropolisSampler.hpp
LDAWBIC_CXXFLAGS = -std=c++11
LDAWBIC_LDADD = LDA-BOWFileParser.o
