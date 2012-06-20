#include "Pairing.hpp"

#include "strtk/strtk.hpp"
#include "Read.hpp"

void
Pairing::parse_token(const string& t)
{
  strtk::std_string::token_list_type token_list;
  strtk::split("=", t, back_inserter(token_list));
  strtk::std_string::token_list_type::iterator itr = token_list.begin();
  if (itr == token_list.end()) {
    cerr << "could not parse pairing token: " << t << endl;
    exit(1);
  }
  string key(itr->first, itr->second);
  ++itr;
  if (itr == token_list.end()) {
    cerr << "could not parse pairing token: " << t << endl;
    exit(1);
  }
  int val=atoi(itr->first);
  if (!key.compare("paired")) {
    paired = val;
  } else if (!key.compare("st_diff")) {
    st_diff = val;
  } else if (!key.compare("min")) {
    min = val;
  } else if (!key.compare("max")) {
    max = val;
  } else if (!key.compare("mean")) {
    mean = val;
  } else if (!key.compare("stddev")) {
    stddev = val;
  } else if (!key.compare("r1_len")) {
    r_len[0] = val;
  } else if (!key.compare("r2_len")) {
    r_len[1] = val;
  } else {
    cerr << "could not parse pairing token: " << t << endl;
    exit(1);
  }    
}

Pairing::Pairing(const string& s)
{
  strtk::std_string::token_list_type token_list;
  strtk::split(",", s, back_inserter(token_list));
  strtk::std_string::token_list_type::iterator itr = token_list.begin();

  while (itr != token_list.end()) {
    parse_token(string(itr->first, itr->second));
    ++itr;
  }
}

int
Pairing::get_mp_st(const Mapping& mapping, int r_st) const
{
  int st = (r_st + mapping.st) % 2;
  return (st + st_diff) % 2;  
}

bool
Pairing::is_mp_downstream(const Mapping& mapping, int r_st) const
{
  assert(mapping.qr != NULL);
  int st = (r_st + mapping.st) % 2;
  int nip = mapping.qr->nip;
  int mp_st = get_mp_st(mapping, r_st);
  int sign_5p_diff;

  if (nip == 0) {
    sign_5p_diff = (st == 0? 1 : -1);
  } else {
    sign_5p_diff = (mp_st == 0? -1 : 1);
  }

  return mean * sign_5p_diff > 0;
}

vector<Interval<long long int> >
Pairing::get_mp_pos(const Mapping& mapping, int r_st) const
{
  vector<Interval<long long int> > result;

  assert(mapping.qr != NULL);
  assert(mapping.qr->mp != NULL);
  int st = (r_st + mapping.st) % 2;
  int nip = mapping.qr->nip;
  int mp_st = get_mp_st(mapping, r_st);
  int mp_rlen = mapping.qr->mp->len;

  long long int pos_5p = mapping.dbPos[st];
  int sign_5p_diff;
  if (nip == 0) {
    sign_5p_diff = (st == 0? 1 : -1);
  } else {
    sign_5p_diff = (mp_st == 0? -1 : 1);
  }

  Interval<long long int> mp_pos_5p;
  mp_pos_5p[0] = pos_5p + sign_5p_diff * min;
  mp_pos_5p[1] = pos_5p + sign_5p_diff * max;
  mp_pos_5p.sort();

  int sign_mp_st = (mp_st == 0? 1 : -1);
  Interval<long long int> tmp;
  tmp[0] = mp_pos_5p[0];
  tmp[1] = mp_pos_5p[0] + sign_mp_st * (mp_rlen - 1);
  tmp.sort();
  result.push_back(tmp);
  tmp[0] = mp_pos_5p[1];
  tmp[1] = mp_pos_5p[1] + sign_mp_st * (mp_rlen - 1);
  tmp.sort();
  result.push_back(tmp);
    
  return result;
}

bool
Pairing::pair_concordant(const Mapping& mapping0, int r0_st,
			 const Mapping& mapping1, int r1_st) const
{
  //int st = (r0_st + mapping0.st) % 2;
  int mp_st = get_mp_st(mapping0, r0_st);
  vector<Interval<long long int> > mp_pos = get_mp_pos(mapping0, r0_st);

  return (r1_st + mapping1.st) % 2 == mp_st and
    mapping1.dbPos[0] >= mp_pos[0][0] and mapping1.dbPos[0] <= mp_pos[1][0];
}

ostream&
operator <<(ostream& ostr, const Pairing& pairing)
{
  ostr << "paired=" << (pairing.paired? 1 : 0) <<
    ",st_diff=" << pairing.st_diff <<
    ",min=" << pairing.min <<
    ",max=" << pairing.max <<
    ",mean=" << pairing.mean <<
    ",stddev=" << pairing.stddev;
  return ostr;
}
