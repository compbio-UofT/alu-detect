#include "SamMapping.hpp"

#include <cstdlib>
#include <iostream>
#include <map>

#include "strtk/strtk.hpp"
#include "globals.hpp"


ExtraSamField::ExtraSamField(const string& s)
{
  if (s.length() < 5 || s[2] != ':' || s[4] != ':') {
    cerr << "invalid SAM field: " << s << endl;
    exit(1);
  }
  key = s.substr(0, 2);
  type = s.substr(3, 1);
  value = s.substr(5);
}

ostream&
operator <<(ostream& ostr, const ExtraSamField& extraSamField)
{
  ostr << extraSamField.key << ":" << extraSamField.type << ":" << extraSamField.value;
  return ostr;
}


SamMapping::SamMapping(const string& s, SQDict* dict, bool add_to_dict)
{
  strtk::std_string::token_list_type token_list;
  strtk::split("\t", s, back_inserter(token_list));
  strtk::std_string::token_list_type::iterator itr = token_list.begin();
  
  if (itr == token_list.end()) { cerr << "invalid SAM line: " << s << endl; exit(1); }
  name = string(itr->first, itr->second);
  ++itr;

  if (itr == token_list.end()) { cerr << "invalid SAM line: " << s << endl; exit(1); }
  flags = bitset<32>(atol(itr->first));
  ++itr;

  if (itr == token_list.end()) { cerr << "invalid SAM line: " << s << endl; exit(1); }
  string tmp = string(itr->first, itr->second);
  //cerr << "got dbName=[" << tmp << "]" << endl;
  if (!tmp.compare("*")) {
    db = NULL;
  } else {
    db = &(*dict)[tmp];
    if (db->name.length() == 0) {
      if (add_to_dict) {
	if (global::verbosity > 0)
	  clog << "adding contig [" << tmp << "]\n";
	db->name = tmp;
	db->idx = dict->size() - 1;
      } else {
	cerr << "error: missing sequence for contig [" << tmp
	     << "] referred to in mapping [" << s << "]" << endl;
	exit(1);
      }
    }
    //cerr << "db now:" << *db << endl;
  }
  ++itr;

  if (itr == token_list.end()) { cerr << "invalid SAM line: " << s << endl; exit(1); }
  dbPos = atoll(itr->first);
  ++itr;

  if (itr == token_list.end()) { cerr << "invalid SAM line: " << s << endl; exit(1); }
  mqv = (int)atol(itr->first);
  ++itr;

  if (itr == token_list.end()) { cerr << "invalid SAM line: " << s << endl; exit(1); }
  cigar = string(itr->first, itr->second);
  ++itr;

  if (itr == token_list.end()) { cerr << "invalid SAM line: " << s << endl; exit(1); }
  tmp = string(itr->first, itr->second);
  if (!tmp.compare("*")) {
    mp_db = NULL;
  } else if (!tmp.compare("=")) {
    mp_db = db;
  } else {
    mp_db = &(*dict)[tmp];
    if (mp_db->name.length() == 0) {
      if (add_to_dict) {
	if (global::verbosity > 0)
	  clog << "adding contig [" << tmp << "]\n";
	mp_db->name = tmp;
	mp_db->idx = dict->size() - 1;
      } else {
	cerr << "error: missing sequence for contig [" << tmp
	     << "] referred to in mapping [" << s << "]" << endl;
	exit(1);
      }
    }
  }
  ++itr;

  if (itr == token_list.end()) { cerr << "invalid SAM line: " << s << endl; exit(1); }
  mp_dbPos = atoll(itr->first);
  ++itr;

  if (itr == token_list.end()) { cerr << "invalid SAM line: " << s << endl; exit(1); }
  tLen = atoll(itr->first);
  ++itr;

  if (itr == token_list.end()) { cerr << "invalid SAM line: " << s << endl; exit(1); }
  seq = DNASequence(itr->first, itr->second);
  ++itr;

  if (itr == token_list.end()) { cerr << "invalid SAM line: " << s << endl; exit(1); }
  qvString = QVString(itr->first, itr->second);
  ++itr;

  rest = vector<ExtraSamField>();
  while (itr != token_list.end()) {
    rest.push_back(ExtraSamField(string(itr->first, itr->second)));
    ++itr;
  }

  if (!flags[0]) {
    nip = 0;
  } else if (flags[6]) {
    nip = 0;
  } else {
    nip = 1;
  }

  mapped = (!flags[2]);

  if (mapped) {
    st = (flags[4]? 1 : 0);
  }
}

ostream&
operator <<(ostream& ostr, const SamMapping& samMapping)
{
  ostr << "name=[" << samMapping.name << "]"
       << " flags=[" << samMapping.flags.to_string() << "]";
  ostr << " db=[";
  if (samMapping.db == NULL)
    ostr << "*";
  else
    ostr << *samMapping.db;
  ostr << "]";
  ostr << " dbPos=[" << samMapping.dbPos << "]"
       << " mqv=[" << samMapping.mqv << "]"
       << " cigar=[" << samMapping.cigar << "]";
  ostr << " mp_db=[";
  if (samMapping.mp_db == NULL)
    ostr << "*";
  else
    ostr << *samMapping.mp_db;
  ostr << "]";
  ostr << " mp_dbPos=[" << samMapping.mp_dbPos << "]"
       << " tLen=[" << samMapping.tLen << "]"
       << " seq=[" << samMapping.seq << "]"
       << " qvString=[" << samMapping.qvString << "]";
  vector<ExtraSamField>::const_iterator it;
  for (it = samMapping.rest.begin(); it != samMapping.rest.end(); ++it) {
    ostr << " [" << *it << "]";
  }

  return ostr;
}

Pairing *
get_pairing_from_SamMapping(const SamMapping& m)
{
  RGDict::iterator it;
  for (size_t i = 0; i < m.rest.size(); ++i) {
    if (m.rest[i].key == "RG") {
      it = global::rg_dict.find(m.rest[i].value);
      if (it == global::rg_dict.end()) {
	cerr << "error: no pairing info for RG: " << m << endl;
	exit(1);
      }
      return &it->second;
    }
  }

  it = global::rg_dict.find(global::default_rg);
  if (it == global::rg_dict.end()) {
    cerr << "error: no pairing info for (default) RG: " << m << endl;
    exit(1);
  }
  return &it->second;
}
