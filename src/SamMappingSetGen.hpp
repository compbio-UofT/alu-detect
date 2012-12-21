#ifndef SamMappingSetGen_hpp_
#define SamMappingSetGen_hpp_

using namespace std;

#include <istream>
#include <string>
#include <vector>

#include "SamMapping.hpp"


class SamMappingSetGen
{
public:
  istream *istr_;
  string (*cloneNameParser_)(const string&);
  void (*headerLineHook_)(const string&);
  SQDict *dict_;
  bool add_to_dict_;

  pair<string,vector<SamMapping> >* next_;

  SamMappingSetGen(istream* istr,
		   string (*cloneNameParser)(const string&),
		   void (*headerLineHook)(const string&),
		   SQDict* dict, bool add_to_dict)
    : istr_(istr),
      cloneNameParser_(cloneNameParser),
      headerLineHook_(headerLineHook),
      dict_(dict),
      add_to_dict_(add_to_dict),
      next_(NULL) {}

  pair<string,vector<SamMapping> >* get_next();
};


#endif
