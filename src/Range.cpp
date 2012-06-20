#include "Range.hpp"

#include <cstdlib>
#include <iostream>


Range::Range(const string& s)
 : db(NULL),
   start(-1),
   end(-1)
{
  size_t i = s.find(':');
  dbName = s.substr(0, i);
  if (i < s.length()) {
    size_t j = s.find('-', i + 1);
    if (j >= s.length()) {
      cerr << "invalid range specification: " << s << endl;
      exit(1);
    }
    start = atoll(&s.c_str()[i + 1]);
    if (start <= 0) {
      cerr << "invalid range specification: start<=0: " << s << endl;
      exit(1);
    }
    --start;
    end = atoll(&s.c_str()[j + 1]);
    if (end <= 0) {
      cerr << "invalid range specification: end<=0: " << s << endl;
      exit(1);
    }
    --end;
    if (start > end) {
      cerr << "invalid range specification: start>end: " << s << endl;
      exit(1);
    }
  }
}

void
Range::parse(const SQDict& dict)
{
  SQDict::const_iterator it = dict.find(dbName);
  if (it == dict.end()) {
    cerr << "invalid range restriction: " << dbName << " not found" << endl;
    exit(1);
  }
  db = &it->second;
  if (start < 0) {
    start = 0;
  } else if (start > db->len - 1) {
    cerr << "invalid range restriction: start=" << start + 1
	 << " larger than contig length=" << db->len << endl;
    exit(1);
  }
  if (end < 0) {
    end = db->len - 1;
  } else if (end > db->len - 1) {
    cerr << "invalid range restriction: end=" << end + 1
	 << " larger than contig length=" << db->len << endl;
    exit(1);
  }
}
