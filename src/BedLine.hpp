#ifndef BedLine_hpp_
#define BedLine_hpp_

using namespace std;

#include <ostream>
#include <string>

#include "DNASequence.hpp"
#include "Interval.hpp"


class BedLine
{
public:
  string line;
  Contig * db;
  Interval<long long> pos;
  int support;
  int skip_left;
  int len_left;
  int skip_right;
  int len_right;

  BedLine() : db(NULL), support(0), skip_left(0), len_left(0), skip_right(0), len_right(0) {}
  BedLine(const string&, SQDict*);

};


#endif
