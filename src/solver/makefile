# @Author: jose
# @Date:   2019-04-05T15:17:03+02:00
# @Last modified by:   jose
# @Last modified time: 2019-04-11T11:10:29+02:00



# Include Petsc-defined variables
#include ${PETSC_DIR}/conf/variables
#include ${PETSC_DIR}/conf/rules

#STATIC LINKING: -static-libgcc -static -static-libstdc++
ARCH := 64

CPP_FLAGS_32 := -m32
CPP_FLAGS_64 := -m64

CCOMP := gcc
CPPCOMP := g++

# CPPFLAGS = $(CPP_FLAGS_$(ARCH)) -ggdb -c -O3 -w -std=c++11
# CFLAGS = $(CPP_FLAGS_$(ARCH)) -ggdb -c -O3 -w -fno-strict-aliasing

CPPFLAGS = $(CPP_FLAGS_$(ARCH)) -c -O3 -w -std=c++11
CFLAGS = $(CPP_FLAGS_$(ARCH)) -c -O3 -w -fno-strict-aliasing

PR_FOLDER = reader/
PR_SOURCES := $(wildcard $(PR_FOLDER)*.cpp)
PR_HEADERS := $(wildcard $(PR_FOLDER)*.h)
PR_OBJECTS := $(patsubst %.cpp,%.o,$(PR_SOURCES))

CO_FOLDER = common/
CO_SOURCES := $(wildcard $(CO_FOLDER)*.cpp)
CO_HEADERS := $(wildcard $(CO_FOLDER)*.h)
CO_OBJECTS := $(patsubst %.cpp,%.o,$(CO_SOURCES))

LO_FOLDER = reader/loader/
LO_SOURCES := $(wildcard $(LO_FOLDER)*.cpp)
LO_HEADERS := $(wildcard $(LO_FOLDER)*.h)
LO_OBJECTS := $(patsubst %.cpp,%.o,$(LO_SOURCES))

EL_FOLDER = reader/elements/
EL_SOURCES := $(wildcard $(EL_FOLDER)*.cpp)
EL_HEADERS := $(wildcard $(EL_FOLDER)*.h)
EL_OBJECTS := $(patsubst %.cpp,%.o,$(EL_SOURCES))

PA_FOLDER = reader/particles/
PA_SOURCES := $(wildcard $(PA_FOLDER)*.cpp)
PA_HEADERS := $(wildcard $(PA_FOLDER)*.h)
PA_OBJECTS := $(patsubst %.cpp,%.o,$(PA_SOURCES))

SO_FOLDER = solver/
SO_SOURCES := $(wildcard $(SO_FOLDER)*.cpp)
SO_HEADERS := $(wildcard $(SO_FOLDER)*.h)
SO_OBJECTS := $(patsubst %.cpp,%.o,$(SO_SOURCES))

SY_FOLDER = solver/system_of_eq/
SY_SOURCES := $(wildcard $(SY_FOLDER)*.cpp)
SY_HEADERS := $(wildcard $(SY_FOLDER)*.h)
SY_OBJECTS := $(patsubst %.cpp,%.o,$(SY_SOURCES))

SP_FOLDER = solver/sparse/
SP_SOURCES := $(wildcard $(SP_FOLDER)*.cpp)
SP_HEADERS := $(wildcard $(SP_FOLDER)*.h)
SP_OBJECTS := $(patsubst %.cpp,%.o,$(SP_SOURCES))

PO_FOLDER = export/
PO_SOURCES := $(wildcard $(PO_FOLDER)*.cpp)
PO_HEADERS := $(wildcard $(PO_FOLDER)*.h)
PO_OBJECTS := $(patsubst %.cpp,%.o,$(PO_SOURCES))

DiscreteModel: main.o $(PR_OBJECTS) $(CO_OBJECTS) $(LO_OBJECTS) $(EL_OBJECTS) $(PA_OBJECTS) $(SO_OBJECTS) $(SY_OBJECTS) $(SP_OBJECTS) $(PO_OBJECTS)
	$(CPPCOMP) $(PR_OBJECTS) $(CO_OBJECTS) $(LO_OBJECTS) $(EL_OBJECTS) $(PA_OBJECTS) $(SO_OBJECTS) $(SY_OBJECTS) $(SP_OBJECTS) $(PO_OBJECTS) main.o -m${ARCH}  -o DiscreteModel
	cp DiscreteModel ../
	@echo "***********************  B U I L D    D O N E  ************************"

main.o: main.cpp
	$(CPPCOMP) ${CPPFLAGS} $<

%.o: $(PR_FOLDER)%.cpp
	$(CPPCOMP) ${CPPFLAGS} $@ $<

%.o: $(CO_FOLDER)%.cpp
	$(CPPCOMP) ${CPPFLAGS} $@ $<

%.o: $(LO_FOLDER)%.cpp
	$(CPPCOMP) ${CPPFLAGS} $@ $<

%.o: $(EL_FOLDER)%.cpp
	$(CPPCOMP) ${CPPFLAGS} $@ $<

%.o: $(PA_FOLDER)%.cpp
	$(CPPCOMP) ${CPPFLAGS} $@ $<

%.o: $(SO_FOLDER)%.cpp
	$(CPPCOMP) ${CPPFLAGS} $@ $<

%.o: $(SY_FOLDER)%.cpp
	$(CPPCOMP) ${CPPFLAGS} $@ $<

%.o: $(SP_FOLDER)%.cpp
	$(CPPCOMP) ${CPPFLAGS} $@ $<

%.o: $(PO_FOLDER)%.cpp
	$(CPPCOMP) ${CPPFLAGS} $@ $<

clean_all:
	# rm -rf *o
    # find DiscreteModel -type f -exec rm {} + || echo " "
	find . $(PR_OBJECTS) $(CO_OBJECTS) $(LO_OBJECTS) $(EL_OBJECTS) $(PA_OBJECTS) $(SO_OBJECTS) $(SY_OBJECTS) $(SP_OBJECTS) $(PO_OBJECTS) -type f -name '*.o' -exec rm {} + || echo " "
	@echo "***********************  C L E A N    D O N E  ************************"
