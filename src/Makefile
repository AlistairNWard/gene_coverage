CXX=g++
C=gcc

# Compiler flags
CFLAGS=-O3 -D_FILE_OFFSET_BITS=64 -g

BAMTOOLS_ROOT=../bamtools

LIBS=-L./ -L$(BAMTOOLS_ROOT)/lib -lz -lm
INCLUDE=-I$(BAMTOOLS_ROOT)/src

all: ../bin/coverage
debug: ../bin/coverage

# builds bamtools static lib, and copies into root
$(BAMTOOLS_ROOT)/lib/libbamtools.a:
	cd $(BAMTOOLS_ROOT) && mkdir -p build && cd build && cmake .. && $(MAKE)

# Objects
OBJECTS=dataProcessing.o \
	mean.o \
	pileup.o \
	$(BAMTOOLS_ROOT)/lib/libbamtools.a

# Executables
coverage ../bin/coverage: coverage.o mean.o pileup.o dataProcessing.o $(OBJECTS)
	$(CXX) $(CFLAGS) $(INCLUDE) coverage.o $(OBJECTS) -o ../bin/coverage $(LIBS)

# Objects
dataProcessing.o: dataProcessing.cpp
	$(CXX) $(CFLAGS) $(INCLUDE) -c dataProcessing.cpp

coverage.o: coverage.cpp $(BAMTOOLS_ROOT)/lib/libbamtools.a
	$(CXX) $(CFLAGS) $(INCLUDE) -c coverage.cpp

pileup.o: pileup.cpp $(BAMTOOLS_ROOT)/lib/libbamtools.a
	$(CXX) $(CFLAGS) $(INCLUDE) -c pileup.cpp

mean.o: mean.cpp
	$(CXX) $(CFLAGS) $(INCLUDE) -c mean.cpp

clean:
	-@rm *.o
	-@rm ../bin/*
