#ifndef SamMapping_hpp_
#define SamMapping_hpp_

using namespace std;

#include <ostream>
#include <string>
#include <vector>
#include <bitset>

#include "globals.hpp"
#include "DNASequence.hpp"


class ExtraSamField
{
public:
  string key;
  string type;
  string value;

  ExtraSamField(const string&);
};

ostream& operator <<(ostream&, const ExtraSamField&);


class SamMapping
{
public:
  string name;
  bitset<32> flags;
  Contig * db;
  long long int dbPos;
  int mqv;
  string cigar;
  Contig * mp_db;
  long long int mp_dbPos;
  long long int tLen;
  DNASequence seq;
  QVString qvString;
  vector<ExtraSamField> rest;

  int nip;
  int st;
  bool mapped;
  bool is_ref;

  SamMapping() {}
  SamMapping(const string&, SQDict*, bool);
};

ostream& operator <<(ostream&, const SamMapping&);
Pairing* get_pairing_from_SamMapping(const SamMapping&);


#endif
