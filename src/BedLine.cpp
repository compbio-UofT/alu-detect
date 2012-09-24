#include "BedLine.hpp"

#include <cstdlib>
#include <limits>

#include "strtk/strtk.hpp"


BedLine::BedLine(const string& s, SQDict* dict)
  : support(0), skip_left(0), len_left(0), skip_right(0), len_right(0)
{
  line = s;

  strtk::std_string::token_list_type token_list;
  strtk::split("\t", s, back_inserter(token_list));
  strtk::std_string::token_list_type::iterator itr = token_list.begin();

  if (itr == token_list.end()) { cerr << "invalid BED line: " << s << endl; exit(1); }
  string tmp(itr->first, itr->second);
  db = &(*dict)[tmp];
  if (db->name.length() == 0) {
    cerr << "warning: adding placeholder for missing contig [" << tmp << "]" << endl;
    db->name = tmp;
    db->len = numeric_limits<long long>::max();
    db->idx = dict->size() - 1;
  }
  ++itr;

  if (itr == token_list.end()) { cerr << "invalid BED line: " << s << endl; exit(1); }
  pos[0] = atoll(itr->first);
  if (pos[0] < 0 or pos[0] >= db->len) { cerr << "error: got position " << pos[0] << " in contig " << db->name << " which has length " << db->len << endl; exit(1); }
  ++itr;

  if (itr == token_list.end()) { cerr << "invalid BED line: " << s << endl; exit(1); }
  pos[1] = atoll(itr->first) - 1;
  if (pos[1] < 0 or pos[1] >= db->len) { cerr << "error: got position " << pos[1] << " in contig " << db->name << " which has length " << db->len << endl; exit(1); }
  if (pos[0] > pos[1]) { cerr << "invalid BED line: " << s << endl; exit(1); }
}
