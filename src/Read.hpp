#ifndef Read_hpp_
#define Read_hpp_

using namespace std;

#include <ostream>
#include <string>
#include <vector>

#include "Mapping.hpp"


typedef bool RepeatEvidence[2];

class Read
{
public:
  string name;
  string seq;
  string qvString;
  Read* mp;
  int len;
  int nip;
  int st;
  vector<Mapping> mapping;

  RepeatEvidence mappedToRepeatSt;

  static int min_read_len;

  Read() { mappedToRepeatSt[0] = false; mappedToRepeatSt[1] = false; }

  bool usable() { return len >= min_read_len; }
  bool mapped() { return mapping.size() > 0; }
  int getRefMappingIdx();
  void flip_st();
  bool fully_mapped() {
    return mapping.size() > 0
      and mapping[0].qrPos[0] == 0 and mapping[mapping.size() - 1].qrPos[1] == len - 1;
  }
};

ostream& operator <<(ostream& ostr, const Read&);

size_t size_below(const Read&);


#endif
