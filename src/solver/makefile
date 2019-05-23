#STATIC LINKING: -static-libgcc -static -static-libstdc++
ARCH := 64

CPP_FLAGS_32 := -m32
CPP_FLAGS_64 := -m64

CCOMP := gcc
CPPCOMP := g++

CPPFLAGS = $(CPP_FLAGS_$(ARCH)) -c -O3 -w -std=c++11
CFLAGS = $(CPP_FLAGS_$(ARCH)) -c -O3 -w -fno-strict-aliasing

SRC_FOLDER = src/
SRC_SOURCES := $(wildcard $(SRC_FOLDER)*.cpp)
SRC_HEADERS := $(wildcard $(SRC_FOLDER)*.h)
SRC_OBJECTS := $(patsubst %.cpp,%.o,$(SRC_SOURCES))

DiscreteModel: $(SRC_OBJECTS)
	$(CPPCOMP) $(SRC_OBJECTS) -m${ARCH}  -o DiscreteModel
	cp DiscreteModel ../
	@echo "***********************  B U I L D    D O N E  ************************"

%.o: $(SRC_FOLDER)%.cpp
	$(CPPCOMP) ${CPPFLAGS} $@ $<

clean_all:
	rm -f $(SRC_OBJECTS)
	@echo "***********************  C L E A N    D O N E  ************************"
