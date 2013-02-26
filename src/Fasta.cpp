#include "Fasta.hpp"

#include <iostream>


void
readFasta(istream& istr, SQDict& dict, bool parseSeqOffset)
{
  string name;
  vector<char> buffer;

  while (true) {
    string s;
    getline(istr, s);
    if (istr.bad()) {
      cerr << "error reading SAM file" << endl;
      exit(1);
    }
    if (istr.eof() or s[0] == '>') {
      if (name.length() > 0) {
	// add previous sequence
	long long int seqOffset = 0;
	if (parseSeqOffset) {
	  int j = name.find(':');
	  if (j < (int)name.length()) {
	    seqOffset = atoll(&name.c_str()[j + 1]) - 1;
	    name = name.substr(0, j);
	  }
	}
	Contig& c = dict[name];
	if (c.name.length() == 0) {
	  c.name = name;
	  c.len = buffer.size();
	  c.idx = dict.size() - 1;
	  c.seq[0] = DNASequence(buffer.begin(), buffer.end());
	  c.seqOffset[0] = seqOffset;
	  cerr << "added contig [" << c.name << "] of length [" << c.len << "]"
	       << " with start offset [" << c.seqOffset[0] << "]" << endl;
	}
      }
      if (istr.eof())
	break;
      int stop = s.find(" ")-1;
      if (stop > 0) {
      	name = s.substr(1,stop);
      } else {
      	name = s.substr(1);
      }
      buffer.clear();
    } else {
      buffer.insert(buffer.end(), s.begin(), s.end());
    }
  }
}
