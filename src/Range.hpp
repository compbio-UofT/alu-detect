#ifndef Range_hpp_
#define Range_hpp_

#include <string>

#include "DNASequence.hpp"


class Range
{
public:
  string dbName;
  const Contig* db;
  long long int start;
  long long int end;

  Range(const string&);

  void parse(const SQDict&);
};


#endif
