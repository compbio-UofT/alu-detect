#include "SamMappingSetGen.hpp"

#include <cstdlib>
#include <cassert>
#include <iostream>


SamMapping*
getMapping(istream& istr, void (*headerLineHook)(const string&), SQDict* dict, bool add_to_dict)
{
  string line;
  while (getline(istr, line)) {
    if (line[0] != '@')
      return new SamMapping(line, dict, add_to_dict);
    if (headerLineHook != NULL)
      headerLineHook(line);
  }
  if (istr.bad()) {
    cerr << "error reading SAM mapping" << endl;
    exit(1);
  }
  return NULL;
}


pair<string,vector<SamMapping> >*
SamMappingSetGen::get_next()
{
  pair<string,vector<SamMapping> >* result = NULL;
  if (next_ != NULL) {
    result = next_;
    next_ = NULL;
  }
  if (istr_->eof()) {
    return result;
  }
  bool done = false;
  while (!done) {
    SamMapping* m = getMapping(*istr_, headerLineHook_, dict_, add_to_dict_);
    if (m == NULL)
      break;
    const string s = cloneNameParser_(m->name);
    if (result == NULL) {
      // start of new set; continue
      result = new pair<string,vector<SamMapping> >(s, vector<SamMapping>(1, *m));
    } else if (!result->first.compare(s)) {
      // same clone, add to existing set; continue
      result->second.push_back(*m);
    } else {
      // new set, save it, report the old one
      next_ = new pair<string,vector<SamMapping> >(s, vector<SamMapping>(1, *m));
      done = true;
    }
    delete m;
  }
  return result;
}
