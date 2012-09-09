SHELL=/bin/bash

ifndef BUILD_TYPE
  BUILD_TYPE=testing
endif

ifdef USE_ICC
  CXX=/opt/intel/cce/10.1.015/bin/icc
endif

ifndef CXXFLAGS
  ifeq ($(BUILD_TYPE), production)
    CXXFLAGS=-O3 -DNDEBUG
  else
    ifeq ($(BUILD_TYPE), testing)
      CXXFLAGS=-g -O3
    else
      CXXFLAGS=-g
    endif
  endif
endif

ifdef USE_ICC
  CXXFLAGS+=-I. -Kc++ -wd383,981,1572 -axP -ipo -openmp -static-intel
else
  CXXFLAGS+=-I. -mmmx -msse -msse2 -fopenmp -Wall -Wextra -Wno-deprecated
endif

#GIT_VERSION=$(shell ./get_git_version)
#override CXXFLAGS+=-DGIT_VERSION=$(GIT_VERSION)

LD=$(CXX)

UNAME := $(shell uname)
ifeq ($(UNAME), Darwin)
LDFLAGS=-lm -lz
else
LDFLAGS=-lm -lz -lrt
endif

LN=ln

all: include-boost get-regions inalign get-coverage zc add-extra-sam-flags move-bins

include-boost:
	if [[ ! -e boost && -d "$(BOOST)" ]]; then ln -s "$(BOOST)" boost; fi

get-regions: get-regions.o Read.o Clone.o Mapping.o DNASequence.o \
	SamMapping.o SamMappingSetGen.o CloneGen.o Pairing.o Fasta.o common.o globals.o \
	inalign_core.o util.o ScoreMatrix.o deep_size.o Range.o Cigar.o gzstream/gzstream.o
	$(LD) $(CXXFLAGS) -o $@ $+ $(LDFLAGS)

inalign: inalign.o inalign_core.o \
	Clone.o Read.o Mapping.o Pairing.o DNASequence.o \
	Fasta.o common.o globals.o ScoreMatrix.o deep_size.o gzstream/gzstream.o
	$(LD) $(CXXFLAGS) -o $@ $+ $(LDFLAGS)

add-extra-sam-flags: add-extra-sam-flags.o globals.o util.o deep_size.o Pairing.o \
	DNASequence.o Read.o Mapping.o CloneGen.o SamMapping.o SamMappingSetGen.o common.o \
	gzstream/gzstream.o Cigar.o Clone.o
	$(LD) $(CXXFLAGS) -o $@ $+ $(LDFLAGS)

get-coverage: get-coverage.o DNASequence.o SamMapping.o Cigar.o globals.o Range.o \
	gzstream/gzstream.o
	$(LD) $(CXXFLAGS) -o $@ $+ $(LDFLAGS)

get-regions.o: get-regions.cpp globals.hpp \
	Read.hpp Clone.hpp Mapping.hpp SamMapping.hpp SamMappingSetGen.hpp \
	CloneGen.hpp Fasta.hpp common.hpp inalign_core.hpp deep_size.hpp Range.hpp Cigar.hpp
	$(CXX) $(CXXFLAGS) -c -o $@ $<

inalign.o: inalign.c inalign_core.hpp Clone.hpp DNASequence.hpp Fasta.hpp \
	common.hpp gzstream/gzstream.o
	$(CXX) $(CXXFLAGS) -c -o $@ $<

inalign_core.o: inalign_core.?pp globals.hpp Clone.hpp DNASequence.hpp \
	Fasta.hpp common.hpp ScoreMatrix.hpp
	$(CXX) $(CXXFLAGS) -c -o $@ $<

add-extra-sam-flags.o: add-extra-sam-flags.cpp globals.hpp \
	Read.hpp Mapping.hpp SamMapping.hpp SamMappingSetGen.hpp Mapping.hpp CloneGen.hpp \
	common.hpp Cigar.hpp
	$(CXX) $(CXXFLAGS) -c -o $@ $<

gzstream/gzstream.o: gzstream/gzstream.C gzstream/gzstream.h
	$(CXX) $(shell echo $(CXXFLAGS) | sed 's/I./&\/gzstream/g') -c -o $@ $<

gflags/gflags.o: gflags/gflags.cc gflags/gflags.h
	$(CXX) $(CXXFLAGS) -c -o $@ $<

Read.o: Read.?pp Mapping.hpp globals.hpp deep_size.hpp
	$(CXX) $(CXXFLAGS) -c -o $@ $<

Clone.o: Clone.?pp Read.hpp DNASequence.hpp Interval.hpp deep_size.hpp
	$(CXX) $(CXXFLAGS) -c -o $@ $<

Mapping.o: Mapping.?pp Read.hpp DNASequence.hpp deep_size.hpp
	$(CXX) $(CXXFLAGS) -c -o $@ $<

DNASequence.o: DNASequence.?pp
	$(CXX) $(CXXFLAGS) -c -o $@ $<

SamMapping.o: SamMapping.?pp DNASequence.hpp
	$(CXX) $(CXXFLAGS) -c -o $@ $<

SamMappingSetGen.o: SamMappingSetGen.?pp SamMapping.hpp
	$(CXX) $(CXXFLAGS) -c -o $@ $<

Pairing.o: Pairing.?pp Mapping.hpp Read.hpp DNASequence.hpp Interval.hpp
	$(CXX) $(CXXFLAGS) -c -o $@ $<

Interval.o: Interval.?pp
	$(CXX) $(CXXFLAGS) -c -o $@ $<

CloneGen.o: CloneGen.?pp Clone.hpp SamMapping.hpp Cigar.hpp Interval.hpp
	$(CXX) $(CXXFLAGS) -c -o $@ $<

Fasta.o: Fasta.?pp DNASequence.hpp globals.hpp
	$(CXX) $(CXXFLAGS) -c -o $@ $<

common.o: common.?pp Clone.hpp
	$(CXX) $(CXXFLAGS) -c -o $@ $<

util.o: util.c util.h
	$(CXX) $(CXXFLAGS) -c -o $@ $<

globals.o: globals.?pp
	$(CXX) $(CXXFLAGS) -c -o $@ $<

ScoreMatrix.o: ScoreMatrix.?pp
	$(CXX) $(CXXFLAGS) -c -o $@ $<

deep_size.o: deep_size.?pp
	$(CXX) $(CXXFLAGS) -c -o $@ $<

Range.o: Range.?pp
	$(CXX) $(CXXFLAGS) -c -o $@ $<

Cigar.o: Cigar.?pp Interval.hpp
	$(CXX) $(CXXFLAGS) -c -o $@ $<

get-coverage.o: get-coverage.?pp Range.hpp
	$(CXX) $(CXXFLAGS) -c -o $@ $<

zc: zc.c
	$(LD) $(CXXFLAGS) -o $@ $+ $(LDFLAGS)

move-bins:
	mv get-regions inalign get-coverage zc add-extra-sam-flags ../bin 

clean:
	rm -f *.o gzstream/gzstream.o get-regions inalign get-coverage zc boost ../bin/get-regions ../bin/inalign ../bin/get-coverage ../bin/zc ../bin/add-extra-sam-flags
