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
    cerr << "warning: cigar string contains no matches: " << cigar << endl;
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
