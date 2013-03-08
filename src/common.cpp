#include "common.hpp"

#include <iostream>

#include "globals.hpp"


string
cloneNameParser(const string& name)
{
  size_t i = name.find_first_of(':') + 1;
  size_t j = name.find_first_of(':', i);
  return string(name.substr(i, j - i));
}

void
fullNameParser(const string& name, Clone& clone, int& nip)
{
  //cerr << "parsing read name: " << name << endl;

  int i = 0;
  int j = name.find(':');
  // first is the num_rgid
  if (global::num_rg_len > 0) {
    string tmp = name.substr(i, j);
    RGDict::iterator it = global::num_rg_dict.find(tmp);
    if (it == global::num_rg_dict.end()) {
      cerr << "error: no pairing info for RG of " << name << endl;
      exit(1);
    }
    clone.pairing = &it->second;
  }

  i = j + 1;
  j = name.find(':', i);
  // next, the clone name, but we already know it
  //cerr << "clone name: " << name.substr(i, j) << endl;

  i = j + 1;
  j = name.find(':', i);
  // this is nip
  nip = atoi(&(name.c_str()[i]));
  if (nip < 0 or nip > 2) {
    cerr << "error parsing read name: " << name << endl;
    exit(1);
  }
  if (nip > 0)
    --nip;
  //cerr << "nip: \"" << name.substr(i, j - i) << "\" -> " << nip << endl;

  i = j + 1;
  j = name.find(':', i);
  // this is len of first mate
  clone.read[0].len = atoi(&(name.c_str()[i]));
  //cerr << "len of read 0: \"" << name.substr(i, j - i) << "\" -> " << clone.read[0].len << endl;

  i = j + 1;
  j = name.find(':', i);
  // this is len of second mate
  clone.read[1].len = atoi(&(name.c_str()[i]));
  //cerr << "len of read 1: \"" << name.substr(i, j - i) << "\" -> " << clone.read[1].len << endl;
  int len = clone.read[nip].len;

  i = j + 1;
  j = name.find(':', i);
  // this is rc indicator
  clone.read[nip].st = atoi(&(name.c_str()[i]));

  i = j + 1;
  j = name.find(':', i);
  // this is seq/qvstring indicator
  if (!name.substr(i, j - i).compare("1")) {
    if (clone.read[nip].seq.length() == 0) {
      clone.read[nip].seq = name.substr(j + 1, len);
      clone.read[nip].qvString = name.substr(j + 1 + len + 1, len);
    }
    if (clone.read[nip].name.length() == 0) {
      clone.read[nip].name = name.substr(0, i) + "0" + name.substr(j + 1 + len + 1 + len);
    }
  } else {
    if (clone.read[nip].name.length() == 0) {
      clone.read[nip].name = name;
    }
  }
}
