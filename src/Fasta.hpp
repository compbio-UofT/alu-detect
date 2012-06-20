#ifndef Fasta_hpp_
#define Fasta_hpp_

#include <istream>

#include "globals.hpp"
#include "DNASequence.hpp"


void readFasta(istream&, SQDict& dict, bool = false);


#endif
