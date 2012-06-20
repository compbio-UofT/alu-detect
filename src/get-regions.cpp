using namespace std;

#include <cstdlib>
#include <iostream>
#include <string>
#include <map>
#include <vector>
#include <omp.h>

#include "gzstream/gzstream.h"
#include "strtk/strtk.hpp"
#include "globals.hpp"
#include "SamMapping.hpp"
#include "Clone.hpp"
#include "Read.hpp"
#include "Mapping.hpp"
#include "SamMappingSetGen.hpp"
#include "CloneGen.hpp"
#include "Fasta.hpp"
#include "common.hpp"
#include "inalign_core.hpp"
#include "deep_size.hpp"
#include "Range.hpp"


int num_threads = 1;


class Chunk
{
public:
  long long int chunk_id;
  int thread_id;
  stringstream* out_str;
  stringstream* err_str;
};

class ChunkComparator
{
public:
  bool operator() (const Chunk& lhs, const Chunk& rhs) { return lhs.chunk_id > rhs.chunk_id; }
};


void
checkSQDict(const string& line, SQDict& dict)
{
  strtk::std_string::token_list_type token_list;
  strtk::split("\t", line, back_inserter(token_list));
  strtk::std_string::token_list_type::iterator itr = token_list.begin();

  string s = string(itr->first, itr->second);
  if (!s.compare("@SQ")) {
    string contig_name = string();
    long long int contig_len = 0;
    while (contig_name.length() == 0 || contig_len == 0) {
      ++itr;
      if (itr == token_list.end()) {
	cerr << "did not find contig name and len in @SQ line: " << line << endl;
	exit(1);
      }
      s = string(itr->first, itr->second);
      if (!s.substr(0, 2).compare("SN")) {
	contig_name = s.substr(3);
      } else if (!s.substr(0, 2).compare("LN")) {
	contig_len = atoll(s.substr(3).c_str());
      }
    }
    Contig* contig = &dict[contig_name];
    if (contig->name.length() == 0) {
      cerr << "error: missing sequence for contig [" << contig_name << "]";
      exit(1);
      //contig->name = contig_name;
      //contig->len = contig_len;
      //contig->idx = dict.size() - 1;
    } else if (contig->len != contig_len) {
      cerr << "error: contig [" << contig->name << "] has different lengths: "
	   << contig->len << " or " << contig_len << endl;
      exit(1);
    }
  }
}

void ref_checkSQDict(const string& line) { checkSQDict(line, global::refDict); }
void rep_checkSQDict(const string& line) { checkSQDict(line, global::repDict); }


void
collect_evidence(istream* refMapIn, istream* repMapIn,
		 vector<pair<int,Clone*> >& l,
		 const vector<Range>& range)
{
  cerr << "---------- collect_evidence start" << endl;

  CloneGen cloneGen(refMapIn, repMapIn, cloneNameParser, fullNameParser,
		    ref_checkSQDict, rep_checkSQDict, &global::refDict, &global::repDict);

  while (true) {
    Clone* c = cloneGen.get_next();
    if (c == NULL)
      break;

    if (c->ref == NULL or c->bp.size() == 0) {
      cerr << "dropping clone [" << c->name << "]" << endl;
      delete c;
      continue;
    }

    if (c->bp[0].pos[0] > c->bp[0].pos[1]) {
      cerr << "error: breakpoint positions in wrong order; this shouldn't happen" << endl;
      exit(1);
    }

    if (range.size() > 0) {
      bool pass = false;
      for (size_t i = 0; i < range.size(); ++i) {
	if (range[i].db == c->ref
	    and ((range[i].start <= c->fragPos[0] and c->fragPos[0] <= range[i].end)
		 or (c->fragPos[0] <= range[i].start and range[i].start <= c->fragPos[1]))) {
	  pass = true;
	  break;
	}
      }
      if (not pass) {
	cerr << "dropping clone [" << c->name << "]: not in range restriction" << endl;
	delete c;
	continue;
      }
    }

    l.push_back(pair<int,Clone*>(0, c));
    l.push_back(pair<int,Clone*>(1, c));
    cerr << *c << endl;
  }

  cerr << "---------- collect_evidence end" << endl;
}


bool
compareCloneEndpoints(const pair<int,Clone*>& p1, const pair<int,Clone*>& p2)
{
  // first order them by contigs
  if (p1.second->ref->idx != p2.second->ref->idx)
    return p1.second->ref->idx < p2.second->ref->idx;

  // next, by relevant breakpoint
  return p1.second->bp[0].pos[p1.first] < p2.second->bp[0].pos[p2.first];
}


void
sort_evidence(vector<pair<int,Clone*> >& l)
{
  cerr << "---------- sort_evidence start" << endl;

  sort(l.begin(), l.end(), compareCloneEndpoints);

  cerr << "---------- sort_evidence end" << endl;
}


void
print_evidence(ostream& ostr, const vector<pair<int,Clone*> >& l)
{
  cerr << "---------- print_evidence start" << endl;

  ostr << "len=" << l.size() << endl;
  ostr << "sizeof(l[0]): " << sizeof(l[0]) << endl;
  ostr << "sizeof(l[0].first): " << sizeof(l[0].first) << endl;
  ostr << "sizeof(l[0].second): " << sizeof(l[0].second) << endl;
  ostr << "sizeof(*(l[0].second)): " << sizeof(*(l[0].second)) << endl;
  ostr << "deep size(l): " << sizeof(l) + size_below(l) << endl;

  size_t tmp = 0;
  for (vector<pair<int,Clone*> >::const_iterator it = l.begin(); it != l.end(); ++it) {
    ostr << it->second->ref->name << " "
	 << it->second->bp[0].pos[it->first] << " " << it->second->bp[0].pos[1 - it->first] << " "
	 << it->second->mappedToRepeatSt[0] << " "
	 << it->second->mappedToRepeatSt[1] << " "
	 << it->second->name << endl;
    if (it->first == 0) {
      tmp += sizeof(*(it->second)) + size_below(*(it->second));
    }
  }

  cerr << "in-mem size of clones: " << tmp << endl;

  cerr << "---------- print_evidence end" << endl;
}


bool
duplicate_clones(Clone* c1, Clone* c2)
{
  return (c1->read[0].seq == c2->read[0].seq
	  and c1->read[1].seq == c2->read[1].seq)
    or (c1->read[0].seq == c2->read[1].seq
	and c1->read[1].seq == c2->read[0].seq);
}

void
filter_evidence(vector<pair<int,Clone*> >& l)
{
  cerr << "---------- filter_evidence start" << endl;

  vector<pair<int,Clone*> > new_l;

  for_iterable(CONCAT(vector<pair<int,Clone*> >), l, it) {
    //for (int i = 0; i < l.size(); ++i) {
    if (it->second->use) {
      new_l.push_back(*it);
    } else if (it->first == 1) {
      delete it->second;
    }
  }

  l = new_l;

  cerr << "---------- filter_evidence end" << endl;
}


void
remove_duplicate_evidence(vector<pair<int,Clone*> >& l)
{
  cerr << "---------- remove_duplicate_evidence start" << endl;

  for (size_t i = 0; i < l.size(); ++i) {
    if (l[i].first == 0) {
      l[i].second->use = true;
    }

    if (i > 0 and l[i].first == 0 and l[i - 1].first == 0
	and duplicate_clones(l[i - 1].second, l[i].second)) {
      // remove this clone, it's a duplicate
      cerr << "removing duplicate [" << l[i].second->name << "]" << endl;
      l[i].second->use = false;
    }
  }

  filter_evidence(l);

  cerr << "---------- remove_duplicate_evidence end" << endl;
}

void
remove_non_active(vector<pair<int,Clone*> >& l)
{
  cerr << "---------- remove_non_active start" << endl;

  set<Clone*> active;
  Contig* crt_contig = NULL;
  long long int last_pos = -1;
  for (size_t i = 0; i < l.size(); i++) {
    int endpoint = l[i].first;
    Clone* clone = l[i].second;
    if (endpoint == 0) { // start of bp range
      active.insert(clone);
      if (clone->mappedToRepeatSt[0] or clone->mappedToRepeatSt[1]) {
	clone->use = true;
	if (crt_contig != clone->ref) {
	  crt_contig = clone->ref;
	  last_pos = clone->bp[0].pos[1];
	} else {
	  last_pos = max(last_pos, clone->bp[0].pos[1]);
	}
	cerr << "setting crt_contig=[" << crt_contig->name << "] last_pos=[" << last_pos
	     << "] because of [" << clone->name << "]" << endl;
	for_iterable(set<Clone*>, active, it) {
	  if (!(*it)->use) {
	    cerr << "including [" << (*it)->name << "] because of [" << clone->name << "]" << endl;
	    (*it)->use = true;
	  }
	}
      } else {
	// no repeat evidence in this clone, but might be in active region
	if (clone->ref == crt_contig and clone->bp[0].pos[0] <= last_pos) {
	  cerr << "including [" << clone->name << "] because last_pos=[" << last_pos << "]" << endl;
	  clone->use = true;
	} else {
	  clone->use = false;
	}
      }
    } else { // end of bp range
      active.erase(clone);
    }
  }

  filter_evidence(l);

  cerr << "---------- remove_non_active end" << endl;
}

void
remove_singletons(vector<pair<int,Clone*> >& l)
{
  cerr << "---------- remove_singletons start" << endl;

  set<Clone*> active;
 
  for (size_t i = 0; i < l.size(); ++i) {
    if (l[i].first == 0) {
      l[i].second->use = true;
      active.insert(l[i].second);
    } else {
      if (l[i].second == l[i - 1].second and active.size() == 1) {
	cerr << "removing singleton [" << l[i].second->name << "]" << endl;
	l[i].second->use = false;
      }
      active.erase(l[i].second);
    }
  }

  filter_evidence(l);

  cerr << "---------- remove_singletons end" << endl;
}

void
split_map_clone(Clone* c, ostream& out_str, ostream& err_str)
{
  if (c->bp.size() != 1) {
    cerr << "split_map_clone: clone has no breakpoints: " << *c << endl;
    exit(1);
  }

  if (c->bp.size() == 1 and c->bp[0].pos[0] < c->bp[0].pos[1]) {
    vector<Clone*> v(1, c);
    run(*c->ref,
	c->fragPos[0], c->fragPos[1], c->bp[0].pos[0], c->bp[0].pos[1],
	global::repeatListAll,
	v,
	1,
	out_str, err_str);
  } else {
    err_str << "warning: bad clone: " << *c << endl;
    c->read[0].mapping.clear();
    c->read[1].mapping.clear();
  }

  c->computePosition(false, false);

  c->use = (c->ref != NULL and c->fully_mapped() and c->bp.size() == 1
	    and (c->mappedToRepeatSt[0] or c->mappedToRepeatSt[1]));
  if (c->ref == NULL) {
    err_str << "removing clone [" << c->name << "]: not mapped to reference" << endl;
  } else if (!c->fully_mapped()) {
    err_str << "removing clone [" << c->name << "]: not fully mapped" << endl;
  } else if (c->bp.size() != 1) {
    err_str << "removing clone [" << c->name << "]: more than 1 bp" << endl;
  } else if (!(c->mappedToRepeatSt[0] or c->mappedToRepeatSt[1])) {
    err_str << "removing clone [" << c->name << "]: not mapped to repeat" << endl;
  }
}


void
split_mappings(vector<pair<int,Clone*> >& l)
{
  cerr << "---------- split_mappings start" << endl;

  vector<Clone*> v;
  for (size_t i = 0; i < l.size(); ++i) {
    if (l[i].first == 0) {
      v.push_back(l[i].second);
    }
  }

  priority_queue<Chunk,vector<Chunk>,ChunkComparator> h;
  long long next_chunk = 0;

  // do this in parallel:
#pragma omp parallel num_threads(num_threads)
  {
    int tid = omp_get_thread_num();
#pragma omp critical(output)
    {
      cerr << "thread " << tid << ": started" << endl;
    }
#pragma omp for schedule(dynamic)
    for (long long i = 0; i < (long long int)v.size(); ++i) {
      Chunk chunk;
      chunk.chunk_id = i;
      chunk.thread_id = tid;
      chunk.out_str = new stringstream();
      chunk.err_str = new stringstream();
      //err_str << "thread " << tid << ":" << endl;
      split_map_clone(v[i], *chunk.out_str, *chunk.err_str);

      if (not global::full_splitmap_log) {
	string s;
	getline(*chunk.err_str, s);
	chunk.err_str->str("");
	*chunk.err_str << s << endl;
      }
      *chunk.err_str << *v[i] << endl;

#pragma omp critical(output)
      {
	h.push(chunk);
	while (h.size() > 0) {
	  chunk = h.top();
	  assert(chunk.chunk_id >= next_chunk);
	  if (chunk.chunk_id > next_chunk) {
	    break;
	  }
	  cout << chunk.out_str->str();
	  cout.flush();
	  cerr << "chunk=" << chunk.chunk_id << " work_thread=" << chunk.thread_id
	       << " print_thread=" << tid << endl;
	  cerr << chunk.err_str->str();
	  cerr.flush();
	  delete chunk.out_str;
	  delete chunk.err_str;
	  h.pop();
	  ++next_chunk;
	}
      }
    }
  }

  filter_evidence(l);

  cerr << "---------- split_mappings end" << endl;
}

/*
void
restrict_bp_ranges(vector<pair<int,Clone*> >& l)
{
  cerr << "---------- restrict_bp_ranges start" << endl;

  set<Clone*> active;
  for (size_t i = 0; i < l.size(); i++) {
    int endpoint = l[i].first;
    Clone* clone = l[i].second;
    if (endpoint == 0) { // start of bp range
      active.insert(clone);
    } else {
      active.erase(clone);
    }
  }

  cerr << "---------- restrict_bp_ranges start" << endl;
}
*/

void
process_chunk(vector<Clone*>& chunk)
{
  vector<Clone*>::const_iterator it = chunk.begin();
  Clone* c = *it;
  long long int outter_start = c->fragPos[0];
  long long int outter_end = c->fragPos[1];
  long long int inner_start = c->bp[0].pos[0];
  long long int inner_end = c->bp[0].pos[1];
  ++it;

  while (it != chunk.end()) {
    c = *it;
    outter_start = min(outter_start, c->fragPos[0]);
    outter_end = max(outter_end, c->fragPos[1]);
    inner_start = min(inner_start, c->bp[0].pos[0]);
    inner_end = max(inner_end, c->bp[0].pos[1]);
    ++it;
  }

  if (inner_start >= inner_end) {
    cerr << "error: abnormal chunk:" << endl;
    for (it = chunk.begin(); it != chunk.end(); ++it) {
      cerr << *(*it) << endl;
    }
    cerr << "outter_start: " << outter_start << endl;
    cerr << "outter_end: " << outter_end << endl;
    cerr << "inner_start: " << inner_start << endl;
    cerr << "inner_end: " << inner_end << endl;
    exit(1);
  }

  cerr << c->ref->name << '\t'
       << inner_start + 1 << '\t'
       << inner_end + 1 << '\t'
       << outter_start + 1 << '\t'
       << outter_end + 1 << '\t'
       << chunk.size() << '\t';
  for (size_t i = 0; i < chunk.size(); ++i) {
    if (i > 0) {
      cerr << ',';
    }
    cerr << chunk[i]->name;
  }
  cerr << endl;

  run(*c->ref,
      outter_start, outter_end, inner_start, inner_end,
      global::repeatListSt[c->mappedToRepeatSt[0]? 0 : 1],
      chunk,
      0,
      cout, cerr);
}


void
partition_region_into_chunks(vector<pair<int,Clone*> >& v, size_t begin, size_t end)
{
  // greedy idea: pick highest peak, run split mapping aggregator, remove all clones mapped, repeat
  int n_chunks = 0;
  while (true) {
    size_t peak_location = 0;
    size_t peak_height = 0;
    size_t crt_height = 0;
    for (size_t i = begin; i < end; ++i) {
      if (not v[i].second->use) {
	continue;
      }
      assert(crt_height > 0 or v[i].first == 0);
      crt_height += (v[i].first == 0? 1 : -1);
      if (crt_height > peak_height) {
	peak_height = crt_height;
	peak_location = i;
      }
    }

    if (peak_height == 0) {
      break;
    }

    set<Clone*> active;
    for (size_t i = begin; i <= peak_location; ++i) {
      if (not v[i].second->use) {
	continue;
      }
      if (v[i].first == 0) {
	active.insert(v[i].second);
      } else {
	active.erase(v[i].second);
      }
    }
    
    ++n_chunks;
    vector<Clone*> active_vector = vector<Clone*>(active.begin(), active.end());
    process_chunk(active_vector);

    for (set<Clone*>::iterator it = active.begin(); it != active.end(); ++it) {
      if ((*it)->fully_mapped()) {
	(*it)->use = false;
      }
    }
  }

  if (n_chunks > 1) {
    cerr << "warning: abnormal region with " << n_chunks << " chunks:" << endl;
    for (vector<pair<int,Clone*> >::iterator it = v.begin() + begin; it != v.begin() + end; ++it) {
      cerr << it->second->ref->name << " "
	   << it->second->bp[0].pos[it->first] << " " << it->second->bp[0].pos[1 - it->first] << " "
	   << it->second->mappedToRepeatSt[0] << " "
	   << it->second->mappedToRepeatSt[1] << " "
	   << it->second->name << endl;
    }
  }
}


void
final_pass(vector<pair<int,Clone*> >& l)
{
  cerr << "---------- final_pass start" << endl;

  // first, split into 2 lists:
  vector<pair<int,Clone*> > v[2];
  for (size_t i = 0; i < l.size(); i++) {
    int endpoint = l[i].first;
    Clone* clone = l[i].second;
    if (clone->mappedToRepeatSt[0] and not clone->mappedToRepeatSt[1]) {
      v[0].push_back(l[i]);
    } else if (clone->mappedToRepeatSt[1] and not clone->mappedToRepeatSt[0]) {
      v[1].push_back(l[i]);
    } else {
      if (endpoint == 0) {
	cerr << "warning: skipping clone mapped to neither or both repeat strands: "
	     << *clone << endl;
      }
    }
  }

  for (int k = 0; k < 2; k++) {
    set<Clone*> active;
    vector<Clone*> chunk;
    //size_t chunkStartIdx = 0;
    for (size_t i = 0; i < v[k].size(); ++i) {
      int endpoint = v[k][i].first;
      Clone* clone = v[k][i].second;
      if (endpoint == 0) { // start of bp range
	active.insert(clone);
	chunk.push_back(clone);
      } else {
	active.erase(clone);
	if (active.size() == 0) {
	  //partition_region_into_chunks(v[k], chunkStartIdx, i + 1);
	  //chunkStartIdx = i + 1;
	  process_chunk(chunk);
	  chunk.clear();
	}
      }
    }
  }

  cerr << "---------- final_pass end" << endl;
}


int
main(int argc, char* argv[])
{
  string progName(argv[0]);

  vector<pair<int,Clone*> > l;
  vector<Range> range;

  char c;
  while ((c = getopt(argc, argv, "p:N:r:v")) != -1) {
    switch (c) {
    case 'p':
      global::pairing = Pairing(string(optarg));
      cerr << "set pairing: " << global::pairing << endl;
      break;
    case 'N':
      num_threads = atoi(optarg);
      break;
    case 'r':
      range.push_back(Range(optarg));
      break;
    case 'v':
      global::full_splitmap_log = true;
      break;
    }
  }

  if (optind + 4 != argc) {
    cerr << "use: " << argv[0]
	 << " [options] <ref_fa> <ref_mapping.sam> <rep_fa> <rep_mapping.sam>" << endl;
    exit(1);
  }

  cerr << "number of threads: " << num_threads << endl;

  igzstream refFaIn(argv[optind]);
  if (refFaIn.bad()) {
    cerr << "error opening reference fasta file: " << argv[optind] << endl;
    exit(1);
  }
  readFasta(refFaIn, global::refDict);
  refFaIn.close();

  igzstream repFaIn(argv[optind + 2]);
  if (repFaIn.bad()) {
    cerr << "error opening repeat fasta file: " << argv[optind + 2] << endl;
    exit(1);
  }
  readFasta(repFaIn, global::repDict);
  repFaIn.close();
  addRCToDict(global::repDict);

  // create list of repeats to investigate
  for (SQDict::iterator it = global::repDict.begin(); it != global::repDict.end(); ++it) {
    global::repeatListAll.push_back(pair<Contig*,int>(&it->second, 0));
    global::repeatListSt[0].push_back(pair<Contig*,int>(&it->second, 0));
  }
  for (SQDict::iterator it = global::repDict.begin(); it != global::repDict.end(); ++it) {
    global::repeatListAll.push_back(pair<Contig*,int>(&it->second, 1));
    global::repeatListSt[1].push_back(pair<Contig*,int>(&it->second, 1));
  }

  // parse range restrictions
  for_iterable(vector<Range>, range, it) {
    it->parse(global::refDict);
  }

  igzstream refMapIn(argv[optind + 1]);
  if (refMapIn.bad()) {
    cerr << "error opening reference mapping file: " << argv[optind + 1] << endl;
    exit(1);
  }
  igzstream repMapIn(argv[optind + 3]);
  if (repMapIn.bad()) {
    cerr << "error opening repeat mapping file:" << argv[optind + 3] << endl;
    exit(1);
  }

  collect_evidence(&refMapIn, &repMapIn, l, range);

  refMapIn.close();
  repMapIn.close();

  sort_evidence(l);
  print_evidence(cerr, l);

  remove_duplicate_evidence(l);
  print_evidence(cerr, l);

  remove_non_active(l);
  remove_singletons(l);
  print_evidence(cerr, l);

  split_mappings(l);

  sort_evidence(l);
  remove_non_active(l);
  remove_singletons(l);
  print_evidence(cerr, l);

  final_pass(l);

  // deallocate memory
  for (size_t i = 0; i < l.size(); ++i) {
    l[i].second->use = false;
  }
  filter_evidence(l);
  print_evidence(cerr, l);

  return 0;
}
