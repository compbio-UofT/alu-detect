#ifndef CloneGen_hpp_
#define CloneGen_hpp_

#include <istream>
#include <string>
#include <vector>

#include "Clone.hpp"
#include "SamMapping.hpp"
#include "SamMappingSetGen.hpp"


class CloneGen
{
public:
  SamMappingSetGen refGen_;
  SamMappingSetGen repGen_;
  void (*fullNameParser_)(const string&, Clone&, int&);

  pair<string,vector<SamMapping> >* next_rep_;

  CloneGen(istream* ref_istr, istream* rep_istr,
	   string (*cloneNameParser)(const string&),
	   void (*fullNameParser)(const string&, Clone&, int&),
	   void (*ref_headerLineHook)(const string&),
	   void (*rep_headerLineHook)(const string&),
	   SQDict* refDict,
	   SQDict* repDict)
    : refGen_(ref_istr, cloneNameParser, ref_headerLineHook, refDict, false),
      repGen_(rep_istr, cloneNameParser, rep_headerLineHook, repDict, false),
      fullNameParser_(fullNameParser),
      next_rep_(NULL) {}

  Clone* get_next();
};

Mapping convert_SamMapping_to_Mapping(const SamMapping&);


#endif
