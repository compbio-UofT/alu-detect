#include "Read.hpp"

#include "globals.hpp"
#include "deep_size.hpp"


int Read::min_read_len = global::min_read_len;


int
Read::getRefMappingIdx()
{
  for (size_t i = 0; i < mapping.size(); ++i) {
    if (mapping[i].is_ref) {
      return i;
    }
  }
  return -1;
}

void
Read::flip_st()
{
  st = 1 - st;
  for (size_t i = 0; i < mapping.size(); ++i) {
    mapping[i].st = 1 - mapping[i].st;
  }
  swap(mappedToRepeatSt[0], mappedToRepeatSt[1]);
  seq = reverseComplement(seq);
  qvString = reverse(qvString);
}


ostream&
operator <<(ostream& ostr, const Read& read)
{
  ostr << "[name=" << read.name
       << ", seq=" << read.seq
       << ", qvString=" << read.qvString
       << ", len=" << read.len
       << ", nip=" << read.nip
       << ", st=" << read.st
       << ", mappedToRepeatSt=[" << read.mappedToRepeatSt[0] << "," << read.mappedToRepeatSt[1]
       << "]";
  for (size_t i = 0; i < read.mapping.size(); ++i) {
    ostr << ", " << read.mapping[i];
  }
  ostr << "]";
  return ostr;
}


size_t
size_below(const Read& r)
{
  return size_below(r.name)
    + size_below(r.seq)
    + size_below(r.qvString)
    + size_below(r.mapping);
}
