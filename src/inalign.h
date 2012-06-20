#ifndef inalign_h_
#define inalign_h_

using namespace std;

#include <string>
#include <vector>

#include "DNASequence.hpp"
#include "Clone.hpp"


void
run(const Contig&,
    long long int, long long int, long long int, long long int,
    const SQDict&,
    vector<pair<Contig*,int> >&,
    vector<Clone>&,
    int);


#endif
