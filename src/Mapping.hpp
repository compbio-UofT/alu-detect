#ifndef Mapping_hpp_
#define Mapping_hpp_

using namespace std;

#include <ostream>

#include "Interval.hpp"


class Read;
class Contig;

class Mapping
{
public:
  Read* qr;
  Contig* db;
  Interval<int> qrPos;
  Interval<long long int> dbPos;
  int st;
  int mqv;
  bool is_ref;
};

ostream& operator <<(ostream&, const Mapping&);

size_t size_below(const Mapping&);


#endif
