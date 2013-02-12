#ifndef Pairing_hpp_
#define Pairing_hpp_

using namespace std;

#include <string>
#include <vector>
#include <map>

#include "Interval.hpp"
#include "Mapping.hpp"


class Pairing
{
public:
  bool paired;
  int st_diff;
  int min;
  int max;
  int mean;
  int stddev;
  int r_len[2];

  Pairing() { paired = false; }
  Pairing(const string&);

  void parse_token(const string&);
  int get_mp_st(const Mapping&, int) const;
  bool is_mp_downstream(const Mapping&, int) const;
  vector<Interval<long long int> > get_mp_pos(const Mapping&, int) const;
  int get_t_len(const Mapping&, int, const Mapping&, int) const;
  bool pair_concordant(const Mapping&, int, const Mapping&, int) const;
};

typedef map<string,Pairing> RGDict;
typedef map<string,string> RGRGDict;

ostream& operator <<(ostream&, const Pairing&);
void load_pairing(istream&, RGDict&, RGDict&, RGRGDict&);


#endif
