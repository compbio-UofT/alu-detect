#ifndef Cigar_hpp_
#define Cigar_hpp_

using namespace std;

#include <cstdlib>
#include <string>
#include <vector>

#include "Interval.hpp"


void parseCigar(const string&, long long, Interval<int>&, Interval<long long>&);


#endif
