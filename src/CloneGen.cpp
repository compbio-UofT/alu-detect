#include "CloneGen.hpp"

#include <iostream>

#include "globals.hpp"
#include "Cigar.hpp"


/*
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
  if (it == cigar_ops.end()) {
    cerr << "cigar string contains no matches: " << cigar << endl;
    exit(1);
  }
  qrPos[0] = x;
  dbPos[0] = y;
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
*/


Mapping
convert_SamMapping_to_Mapping(const SamMapping& samMapping)
{
  Mapping result;

  result.db = samMapping.db;
  result.st = samMapping.st;
  result.mqv = samMapping.mqv;
  parseCigar(samMapping.cigar, samMapping.dbPos, result.qrPos, result.dbPos);

  //cerr << "parsing cigar=" << samMapping.cigar << ", start_dbPos=" << samMapping.dbPos
  //     << " => qrPos=" << result.qrPos << ", dbPos=" << result.dbPos << endl;

  return result;
}


Clone*
CloneGen::get_next()
{
  pair<string,vector<SamMapping> >* next_ref = refGen_.get_next();
  if (next_rep_ == NULL) {
    next_rep_ = repGen_.get_next();
  }
  if (next_ref == NULL) {
    if (next_rep_ == NULL) {
      return NULL;
    } else {
      cerr << "error: clone [" << next_rep_->first
	   << "] has repeat mappings but no reference mappings" << endl;
      exit(1);
    }
  }

  // create clone object containing reference mappings
  Clone* result = new Clone();
  result->name = next_ref->first;
  for_iterable(vector<SamMapping>, next_ref->second, it) {
    int nip;
    fullNameParser_(it->name, *result, nip);
    if (it->mapped) {
      Mapping m = convert_SamMapping_to_Mapping(*it);
      m.qr = &result->read[nip];
      m.is_ref = true;
      result->read[nip].mapping.push_back(m);
    }
  }

  delete next_ref;
  next_ref = NULL;

  for (int nip = 0; nip < 2; nip++) {
    // check only 1 reference mapping per read
    if (result->read[nip].mapping.size() >= 2) {
      cerr << "error: clone [" << result->name
	   << "] has multiple reference mappings for nip [" << nip << "]" << endl;
      exit(1);
    }
    // if missing a read sequence of length 1, add a dummy
    // if missing a longer read sequence, crash
    if (result->read[nip].seq.length() == 0) {
      if (result->read[nip].len == 1) {
	result->read[nip].seq = "N";
	result->read[nip].qvString = "!";
      } else {
	cerr << "error: missing read sequence from clone [" << result->name
	     << "], nip [" << nip << "]" << endl;
	exit(1);
      }
    }
  }

  if (next_rep_ != NULL && !result->name.compare(next_rep_->first)) {
    // merge with reference mappings
    for_iterable (vector<SamMapping>, next_rep_->second, it) {
      if (it->mapped) {
	int nip;
	fullNameParser_(it->name, *result, nip);
	result->read[nip].mappedToRepeatSt[it->st] = true;
      }
    }
    delete next_rep_;
    next_rep_ = NULL;
  }

  // compute fragment position and possible breakpoints
  //result->pairing = &global::pairing;
  result->computePosition(true, true);

  return result;
}
