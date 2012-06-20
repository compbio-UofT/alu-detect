#include "Mapping.hpp"

#include "Read.hpp"
#include "DNASequence.hpp"


ostream&
operator <<(ostream& ostr, const Mapping& mapping)
{
  ostr << "[db=" << mapping.db->name
       << ", qrPos=" << mapping.qrPos
       << ", dbPos=" << mapping.dbPos
       << ", st=" << mapping.st << "]";
  return ostr;
}


size_t
size_below(const Mapping&)
{
  return 0;
}
