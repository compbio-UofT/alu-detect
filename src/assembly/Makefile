CXX=g++

CXXFLAGS=-g -O3 -Wall -Wextra -pedantic-errors -m64 -Werror# -mmmx -msse -msse2

all: get_reads merge filter_reads

get_reads: get_reads.cpp
	$(CXX) $(CXXFLAGS) -o ../../bin/assembly/$@ $<

filter_reads: filter_reads.cpp
	$(CXX) $(CXXFLAGS) -o ../../bin/assembly/$@ $<

Coord.o: Coord.cpp
	$(CXX) $(CXXFLAGS) -c -o $@ $<

AluCall.o: AluCall.cpp AluCall.h
	$(CXX) $(CXXFLAGS) -c -o $@ $<

Call.o: Call.cpp Call.h Coord.o AluCall.o
	$(CXX) $(CXXFLAGS) -c -o $@ $<

CallTree.o: CallTree.cpp CallTree.h Call.o
	$(CXX) $(CXXFLAGS) -c -o $@ $<

merge: merge.cpp CallTree.o
	$(CXX) $(CXXFLAGS) -o ../../bin/assembly/$@ $< CallTree.o Call.o Coord.o AluCall.o

clean_o:
	rm *.o

clean: clean_o
	rm ../../bin/assembly/get_reads
	rm ../../bin/assembly/merge
	rm ../../bin/assembly/filter_reads
