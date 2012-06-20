#ifndef Clone_hpp_
#define Clone_hpp_

using namespace std;

#include <ostream>
#include <string>
#include <vector>

#include "Read.hpp"
#include "DNASequence.hpp"
#include "Interval.hpp"
#include "Pairing.hpp"


class BP
{
public:
  Interval<long long int> pos;
  bool solid;
  bool refAnchor[2];
};


class Clone
{
public:
  string name;
  Read read[2];
  Contig* ref;
  Pairing* pairing;
  Interval<long long int> fragPos;
  vector<BP> bp;
  RepeatEvidence mappedToRepeatSt;

  bool use;

  Clone(const string& _name = string())
    : name(_name), ref(NULL), pairing(NULL) {
    read[0].st = 0;
    read[1].st = 0;
    read[0].nip = 0;
    read[1].nip = 1;
    read[0].mp = &read[1];
    read[1].mp = &read[0];
    mappedToRepeatSt[0] = false;
    mappedToRepeatSt[1] = false;
  }

  void computePosition(bool, bool);
  bool fully_mapped() {
    return (not read[0].usable() or read[0].fully_mapped())
      and (not read[1].usable() or read[1].fully_mapped());
  }
};

ostream& operator <<(ostream&, const Clone&);
pair<int,Clone>* readFastq(istream&, void (*fullNameParser)(const string&, Clone&, int&));
vector<Clone> readAllFastq(istream&,
			   string (*cloneNameParser)(const string&),
			   void (*fullNameParser)(const string&, Clone&, int&));

size_t size_below(const Clone&);


#endif
