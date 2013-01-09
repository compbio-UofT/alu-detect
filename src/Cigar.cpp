#include "Cigar.hpp"

#include <iostream>

#include "globals.hpp"


vector<pair<char,int> >
getCigarOps(const string& cigar)
{
  vector<pair<char,int> > result;
  char* s = (char *)cigar.c_str();
  char* c;
  int l;
  while (true) {
    l = strtol(s, &c, 10);
    if (c == s)
      break;
    result.push_back(pair<char,int>(*c, l));
    s = c + 1;
  }

  //cerr << "getCigarOps: cigar=" << cigar << " =>";
  //for (vector<pair<char,int> >::iterator it = result.begin(); it != result.end(); ++it) {
  //  cerr << " ['" << it->first << "'," << it->second << "]";
  ///}
  //cerr << endl;

  return result;
}


void
parseCigar(const string& cigar, long long int dbPos_start,
	   Interval<int>& qrPos, Interval<long long int>& dbPos)
{
  vector<pair<char,int> > cigar_ops = getCigarOps(cigar);

  int x = 0;
  long long int y = dbPos_start - 1;
  vector<pair<char,int> >::iterator it = cigar_ops.begin();

  while (it != cigar_ops.end()) {
    bool done = false;
    switch (it->first) {
    case 'M':
    case '=':
    case 'X':
      if (it->second >= global::min_tail_len) {
	done = true;
      } else {
	x += it->second;
	y += it->second;
      }
      break;

    case 'I':
    case 'S':
    case 'H':
      x += it->second;
      break;

    case 'D':
    case 'N':
    case 'P':
      y += it->second;
      break;
    }
    if (done)
      break;
    ++it;
  }
  qrPos[0] = x;
  dbPos[0] = y;
  if (it == cigar_ops.end()) {
    cerr << "warning: cigar string doesn't contain enough matches: " << cigar << endl;
    qrPos[1] = x;
    dbPos[1] = y;
    return;
    //exit(1);
  }
  x += it->second;
  y += it->second;
  ++it;
  int tmp_x = 0;
  long long int tmp_y = 0;
  while (it != cigar_ops.end()) {
    switch (it->first) {
    case 'M':
    case '=':
    case 'X':
      tmp_x += it->second;
      tmp_y += it->second;
      if (it->second >= global::min_tail_len) {
	x += tmp_x;
	y += tmp_y;
	tmp_x = 0;
	tmp_y = 0;
      }
      break;

    case 'I':
    case 'S':
    case 'H':
      tmp_x += it->second;
      break;

    case 'D':
    case 'N':
    case 'P':
      tmp_y += it->second;
      break;
    }
    ++it;
  }
  qrPos[1] = x - 1;
  dbPos[1] = y - 1;
}


void
get_tail_insert_size_aux(const vector<pair<char,int> >& cigar_ops, int min_tail_match_len,
			 int start_idx, int end_idx, int step, int& dest)
{
  dest = 0;
  int i = start_idx;
  while (i != end_idx) {
    if (cigar_ops[i].first == 'H'
	or cigar_ops[i].first == 'S'
	or cigar_ops[i].first == 'I') {
      dest += cigar_ops[i].second;
      i += step;
    } else if (cigar_ops[i].first == 'M'
	       or cigar_ops[i].first == '='
	       or cigar_ops[i].first == 'X') {
      int tmp = cigar_ops[i].second;
      int j = i + step;
      while (j != end_idx
	     and (cigar_ops[j].first == 'M'
		  or cigar_ops[j].first == '='
		  or cigar_ops[j].first == 'X')) {
	tmp += cigar_ops[j].second;
	j += step;
      }
      if (tmp < min_tail_match_len) {
	dest += tmp;
	i = j;
      } else { // longer match
	break;
      }
    } else { // DN ops
      break;
    }
  }
}


void
get_tail_insert_size(const string& cigar, int min_tail_match_len, vector<int>& tails)
{
  vector<pair<char,int> > cigar_ops = getCigarOps(cigar);

  get_tail_insert_size_aux(cigar_ops, min_tail_match_len, 0, cigar_ops.size(), 1, tails[0]);
  get_tail_insert_size_aux(cigar_ops, min_tail_match_len, cigar_ops.size() - 1, -1, -1, tails[1]);
}
