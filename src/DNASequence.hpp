#ifndef DNASequence_hpp_
#define DNASequence_hpp_

using namespace std;

#include <ostream>
#include <string>
#include <map>


typedef string DNASequence;

typedef string QVString;

class Contig
{
public:
  string name;
  DNASequence seq[2];
  long long int len;
  long long int seqOffset[2];
  int idx;

  Contig() : len(0) { seqOffset[0] = 0; seqOffset[1] = 0; }
};

typedef map<string,Contig> SQDict;


ostream& operator <<(ostream&, const Contig&);

DNASequence reverseComplement(const DNASequence&);
string reverse(const string&);
void addRCToDict(SQDict&);


#endif
