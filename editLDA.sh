SRC="./src/"
INCLUDE="./src/include/"
gvim -p ${SRC}LDA.cpp ${INCLUDE}GibbsSamplerFromLDA.hpp ${SRC}GibbsSamplerFromLDA.cpp ${INCLUDE}OrdinaryGibbsSamplerFromLDA.hpp ${SRC}OrdinaryGibbsSamplerFromLDA.cpp ${INCLUDE}CollapsedGibbsSamplerFromLDA.hpp ${SRC}CollapsedGibbsSamplerFromLDA.cpp ${INCLUDE}BOWFileParser.hpp ${SRC}BOWFileParser.cpp
