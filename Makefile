COMPILER=g++
BIN=./bin/
SRC=./src/
INCLUDE=./src/include/
LDA_OBJECT=$(BIN)LDA.o $(BIN)GibbsSamplerFromLDA.o $(BIN)OrdinaryGibbsSamplerFromLDA.o $(BIN)CollapsedGibbsSamplerFromLDA.o $(BIN)BOWFileParser.o $(BIN)VariationalBayesEstimatorOnLDA.o
ATM_OBJECT=$(BIN)ATM.o $(BIN)GibbsSamplerFromLDA.o $(BIN)GibbsSamplerFromATM.o $(BIN)CollapsedGibbsSamplerFromATM.o $(BIN)BOWFileParser.o $(BIN)AuthorFileParser.o
LDA_HEADER=$(INCLUDE)GibbsSamplerFromLDA.hpp $(INCLUDE)OrdinaryGibbsSamplerFromLDA.hpp $(INCLUDE)CollapsedGibbsSamplerFromLDA.hpp $(INCLUDE)BOWFileParser.hpp $(INCLUDE)VariationalBayesEstimatorOnLDA.hpp
ATM_HEADER=$(INCLUDE)GibbsSamplerFromLDA.hpp $(INCLUDE)GibbsSamplerFromATM.hpp $(INCLUDE)CollapsedGibbsSamplerFromATM.hpp $(INCLUDE)BOWFileParser.hpp $(INCLUDE)AuthorFileParser.hpp

all:LDA ATM

LDA: $(LDA_OBJECT)
	$(COMPILER) -o $(BIN)LDA $(LDA_OBJECT) -O2 -lboost_program_options
$(BIN)LDA.o: $(SRC)LDA.cpp
	$(COMPILER) -c $(SRC)LDA.cpp -o $(BIN)LDA.o -std=c++11 -O2
$(BIN)LDA.o: $(LDA_HEADER)

ATM: $(ATM_OBJECT)
	$(COMPILER) -o $(BIN)ATM $(ATM_OBJECT) -O2 -lboost_program_options
$(BIN)ATM.o: $(SRC)ATM.cpp
	$(COMPILER) -c $(SRC)ATM.cpp -o $(BIN)ATM.o -std=c++11 -O2
$(BIN)ATM.o: $(ATM_HEADER)

$(BIN)GibbsSamplerFromLDA.o: $(SRC)GibbsSamplerFromLDA.cpp
	$(COMPILER) -c $(SRC)GibbsSamplerFromLDA.cpp -o $(BIN)GibbsSamplerFromLDA.o -std=c++11 -O2
$(BIN)GibbsSamplerFromLDA.o: $(INCLUDE)GibbsSamplerFromLDA.hpp

$(BIN)OrdinaryGibbsSamplerFromLDA.o: $(SRC)OrdinaryGibbsSamplerFromLDA.cpp
	$(COMPILER) -c $(SRC)OrdinaryGibbsSamplerFromLDA.cpp -o $(BIN)OrdinaryGibbsSamplerFromLDA.o -std=c++11 -O2
$(BIN)OrdinaryGibbsSamplerFromLDA.o: $(INCLUDE)GibbsSamplerFromLDA.hpp $(INCLUDE)OrdinaryGibbsSamplerFromLDA.hpp

$(BIN)CollapsedGibbsSamplerFromLDA.o: $(SRC)CollapsedGibbsSamplerFromLDA.cpp
	$(COMPILER) -c $(SRC)CollapsedGibbsSamplerFromLDA.cpp -o $(BIN)CollapsedGibbsSamplerFromLDA.o -std=c++11 -O2
$(BIN)CollapsedGibbsSamplerFromLDA.o: $(INCLUDE)GibbsSamplerFromLDA.hpp $(INCLUDE)CollapsedGibbsSamplerFromLDA.hpp

$(BIN)BOWFileParser.o: $(SRC)BOWFileParser.cpp
	$(COMPILER) -c $(SRC)BOWFileParser.cpp -o $(BIN)BOWFileParser.o -std=c++11 -O2
$(BIN)BOWFileParser.o: $(INCLUDE)BOWFileParser.hpp

$(BIN)AuthorFileParser.o: $(SRC)AuthorFileParser.cpp
	$(COMPILER) -c $(SRC)AuthorFileParser.cpp -o $(BIN)AuthorFileParser.o -std=c++11 -O2
$(BIN)AuthorFileParser.o: $(INCLUDE)AuthorFileParser.hpp

$(BIN)GibbsSamplerFromATM.o: $(SRC)GibbsSamplerFromATM.cpp
	$(COMPILER) -c $(SRC)GibbsSamplerFromATM.cpp -o $(BIN)GibbsSamplerFromATM.o -std=c++11 -O2
$(BIN)GibbsSamplerFromATM.o: $(INCLUDE)GibbsSamplerFromLDA.hpp $(INCLUDE)GibbsSamplerFromATM.hpp

$(BIN)CollapsedGibbsSamplerFromATM.o: $(SRC)CollapsedGibbsSamplerFromATM.cpp
	$(COMPILER) -c $(SRC)CollapsedGibbsSamplerFromATM.cpp -o $(BIN)CollapsedGibbsSamplerFromATM.o -std=c++11 -O2
$(BIN)CollapsedGibbsSamplerFromATM.o: $(INCLUDE)GibbsSamplerFromATM.hpp $(INCLUDE)CollapsedGibbsSamplerFromATM.hpp
$(BIN)VariationalBayesEstimatorOnLDA.o: $(SRC)VariationalBayesEstimatorOnLDA.cpp
	$(COMPILER) -c $(SRC)VariationalBayesEstimatorOnLDA.cpp -o $(BIN)VariationalBayesEstimatorOnLDA.o -std=c++11 -O2
$(BIN)VariationalBayesEstimatorOnLDA.o: $(INCLUDE)VariationalBayesEstimatorOnLDA.hpp $(INCLUDE)VariationalBayesEstimatorOnLDA.hpp
