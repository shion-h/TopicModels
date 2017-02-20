SRC="./src/"
INCLUDE="./src/include/"
gvim -p ${SRC}LDA.cpp ${INCLUDE}GibbsSamplerFromLDA.hpp ${SRC}GibbsSamplerFromLDA.cpp ${INCLUDE}OrdinaryGibbsSamplerFromLDA.hpp ${SRC}OrdinaryGibbsSamplerFromLDA.cpp ${INCLUDE}CollapsedGibbsSamplerFromLDA.hpp ${SRC}CollapsedGibbsSamplerFromLDA.cpp ${SRC}ATM.cpp ${INCLUDE}GibbsSamplerFromATM.hpp ${SRC}GibbsSamplerFromATM.cpp ${INCLUDE}CollapsedGibbsSamplerFromATM.hpp ${SRC}CollapsedGibbsSamplerFromATM.cpp ${INCLUDE}BOWFileParser.hpp ${SRC}BOWFileParser.cpp ${INCLUDE}AuthorFileParser.hpp ${SRC}AuthorFileParser.cpp
