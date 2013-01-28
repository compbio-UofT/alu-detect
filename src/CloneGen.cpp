#include "CloneGen.hpp"

#include <iostream>

#include "globals.hpp"
#include "Cigar.hpp"


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
    if (result->read[nip].seq.length() == 0) {
      if (not it->mapped or it->st == 0) {
	result->read[nip].seq = it->seq;
	result->read[nip].qvString = it->qvString;
      } else {
	result->read[nip].seq = reverseComplement(it->seq);
	result->read[nip].qvString = reverse(it->qvString);
      }
    }
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
