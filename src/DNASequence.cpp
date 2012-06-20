#include "DNASequence.hpp"

#include <cstdlib>
#include <iostream>
#include <vector>


char
complement(char c)
{
  if (c == 'a') return 't';
  if (c == 'A') return 'T';
  if (c == 'c') return 'g';
  if (c == 'C') return 'G';
  if (c == 'g') return 'c';
  if (c == 'G') return 'C';
  if (c == 't') return 'a';
  if (c == 'T') return 'A';
  if (c == 'n') return 'n';
  if (c == 'N') return 'N';
  return '?';
}


ostream&
operator <<(ostream& ostr, const Contig& contig)
{
  ostr << "name=[" << contig.name << "] length=[" << contig.len << "] seq:" <<
    (contig.seq[0].length() > 0 ? "yes" : "no");
  return ostr;
}

DNASequence
reverseComplement(const DNASequence& s)
{
  vector<char> tmp;
  for (string::const_reverse_iterator it = s.rbegin(); it != s.rend(); ++it) {
    tmp.push_back(complement(*it));
  }
  return DNASequence(tmp.begin(), tmp.end());
}

string
reverse(const string& s)
{
  vector<char> tmp;
  for (string::const_reverse_iterator it = s.rbegin(); it != s.rend(); ++it) {
    tmp.push_back(*it);
  }
  return string(tmp.begin(), tmp.end());
}

void
addRCToDict(SQDict& dict)
{
  for (SQDict::iterator it = dict.begin(); it != dict.end(); ++it) {
    Contig& c = it->second;
    c.seq[1] = reverseComplement(c.seq[0]);
    c.seqOffset[1] = c.len - 1 - (c.seqOffset[0] + c.seq[1].length() - 1);
  }
}
