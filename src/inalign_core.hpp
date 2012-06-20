#ifndef inalign_core_h_
#define inalign_core_h_

using namespace std;

#include <ostream>
#include <string>
#include <vector>

#include "DNASequence.hpp"
#include "Clone.hpp"


void
run(const Contig&,
    long long int, long long int, long long int, long long int,
    const vector<pair<Contig*,int> >&,
    vector<Clone*>&,
    int,
    ostream&, ostream&);


#endif
