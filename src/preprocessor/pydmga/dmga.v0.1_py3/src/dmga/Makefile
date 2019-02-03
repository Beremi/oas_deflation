#CONFIGURE THIS
# C++ compiler
CXX=g++

#NO CONFIGURE BELOW THIS LINE

BASE_DIR=$(CURDIR)

#where the 3rd party source code lies
PARTY3RD_INC=$(BASE_DIR)/3rd

#where the 3rd party source code lies
SRC_INC=$(BASE_DIR)

#all include paths
INC = -I$(PARTY3RD_INC) -I$(SRC_INC)

# Flags for the C++ compiler
CFLAGS=-Wall -ansi -pedantic $(INC) -std=c++0x -ggdb -Wsign-compare -DEIGEN_MPL2_ONLY
# todo: -Wsign-compare wyrzucić i poprawić interfejs na styku moj program / stl

#where to put executables/libraries
OUT_DIR=$(BASE_DIR)/../../bin
OUT_EXT=

$(OUT_DIR): 
	@mkdir -p $(OUT_DIR)

# Tell make that these are phony targets
.PHONY: all clean default voropp test dppc sasa_test

#default target for eclipse
default: all

#build voro++ if necessary
voropp:
	$(MAKE) -C 3rd/voro all

#build tests
test:
	g++ $(CFLAGS) test.cpp -L3rd/voro -lvoro++ -o $(OUT_DIR)/test/run_tests$(OUT_EXT)
	
#build sasa tmp test
sasa_test:
	g++ $(CFLAGS) sasa_test.cpp -L3rd/voro -lvoro++ -o $(OUT_DIR)/test/sasa_test$(OUT_EXT)
	
#build current program for publication
dppc:
	g++ $(CFLAGS) dppc.cpp -L3rd/voro -lvoro++ -o $(OUT_DIR)/dppc/dppc$(OUT_EXT)
	
#build my basic application
dmgalpha:
	g++ $(CFLAGS) dmgalpha.cpp -L3rd/voro -lvoro++ -o $(OUT_DIR)/dmgalpha$(OUT_EXT)	

# Build all of the executable files
all: voropp dmgalpha
	
	
dmgaclean:
	rm -f $(OUT_DIR)/dmgalpha$(OUT_EXT)
	rm -f $(OUT_DIR)/test/run_tests$(OUT_EXT)
	rm -f $(OUT_DIR)/dppc/dppc$(OUT_EXT)
	
voroclean:
	$(MAKE) -C 3rd/voro clean

# Clean up the executable files
clean: voroclean dmgaclean

	
