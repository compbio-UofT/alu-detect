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
#include "Read.hpp"
#include "Mapping.hpp"
#include "SamMapping.hpp"
#include "SamMappingSetGen.hpp"
#include "common.hpp"


int num_threads = 1;

string (*cnp)(const string&);
void (*fnp)(const string&, Clone&, int&);

class Filter
{
public:
  vector<unsigned long> must_have;
  vector<unsigned long> must_not_have;
  bool stop_on_hit;
  string dest_file;
};

map<string,FILE *> file_map;
vector<Filter> filter_vector;

class Chunk
{
public:
  long long int chunk_id;
  int thread_id;
  map<string,stringstream*> out_str;
  stringstream* err_str;
};

class ChunkComparator
{
public:
  bool operator() (const Chunk& lhs, const Chunk& rhs) { return lhs.chunk_id > rhs.chunk_id; }
};


void
addSQToRefDict(const string& line)
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
    Contig* contig = &global::refDict[contig_name];
    if (contig->name.length() == 0) {
      //cerr << "error: missing sequence for contig [" << contig_name << "]";
      //exit(1);
      contig->name = contig_name;
      contig->len = contig_len;
      contig->idx = global::refDict.size() - 1;
      if (global::verbosity > 0)
	clog << "added contig [" << contig->name << "] of length [" << contig->len << "]\n";
    } else {
      cerr << "error: contig [" << contig->name << "] already exists!" << endl;
      exit(1);
    }
  }
  //cout << line << '\n';
}

void
print_sam_mapping(ostream& os, const SamMapping& m)
{
  os << m.name
     << '\t' << m.flags.to_ulong()
     << '\t' << (m.db != NULL? m.db->name : "*")
     << '\t' << m.dbPos
     << '\t' << m.mqv
     << '\t' << m.cigar
     << '\t' << (m.mp_db != NULL?
		 (m.mp_db == m.db? "=" : m.mp_db->name) : "*")
     << '\t' << m.mp_dbPos
     << '\t' << m.tLen
     << '\t' << m.seq
     << '\t' << m.qvString;
  for (size_t j = 0; j < m.rest.size(); ++j) {
    os << '\t' << m.rest[j].key << ':' << m.rest[j].type << ':' << m.rest[j].value;
  }
  os << '\n';
}

void
process_mapping_set(const string& s, vector<SamMapping>& v,
		    map<string,stringstream*>& out_str) //, ostream* err_str)
{
  if ((global::rg_dict.size() == 0 and v.size() != 1)
      or (global::rg_dict.size() > 0 and v.size() != 2)) {
    cerr << "incorrect number of mappings for clone [" << s << "]" << endl;
    exit(1);
  }

  if (global::rg_dict.size() == 0) {
    for (size_t i = 0; i < filter_vector.size(); ++i) {
      Filter& f = filter_vector[i];
      if (((v[0].flags.to_ulong() & f.must_have[0]) == f.must_have[0])
	  && ((v[0].flags.to_ulong() & f.must_not_have[0]) == 0)) {
	if (out_str.count(f.dest_file) == 0) {
	  out_str[f.dest_file] = new stringstream();
	}
	print_sam_mapping(*(out_str[f.dest_file]), v[0]);
	if (f.stop_on_hit) {
	  break;
	}
      }
    }
  } else {
    for (size_t i = 0; i < filter_vector.size(); ++i) {
      Filter& f = filter_vector[i];
      if (((v[0].flags.to_ulong() & f.must_have[0]) == f.must_have[0])
	  && ((v[0].flags.to_ulong() & f.must_not_have[0]) == 0)
	  && ((v[1].flags.to_ulong() & f.must_have[1]) == f.must_have[1])
	  && ((v[1].flags.to_ulong() & f.must_not_have[1]) == 0)) {
	if (out_str.count(f.dest_file) == 0) {
	  out_str[f.dest_file] = new stringstream();
	}
	print_sam_mapping(*(out_str[f.dest_file]), v[0]);
	print_sam_mapping(*(out_str[f.dest_file]), v[1]);
	if (f.stop_on_hit) {
	  break;
	}
      }
    }
  }
}

string
default_cnp(const string& s)
{
  return string(s);
}

void
parse_bitmask(const string& s, unsigned long& must_have, unsigned long& must_not_have)
{
  must_have = 0;
  must_not_have = 0;
  size_t i = s.find('/');
  string one_side = s.substr(0, i);
  stringstream ss(one_side);
  ss >> hex >> must_have;

  if (i < s.npos) {
    one_side = s.substr(i + 1);
    ss.clear();
    ss.str(one_side);
    ss >> hex >> must_not_have;
  }
  if (must_have & must_not_have) {
    cerr << "invalid bitmask: " << s << " (some bits are both set and unset)" << endl;
    exit(1);
  }
}

void
add_filter(const string& s)
{
  Filter f;

  size_t i = s.find(':');
  if (i == s.npos) {
    cerr << "invalid filter: " << s << endl;
    exit(1);
  }
  string conditions = s.substr(0, i);
  f.dest_file = s.substr(i + 1);

  if (global::rg_dict.size() == 0) {
    f.must_have = vector<unsigned long>(1);
    f.must_not_have = vector<unsigned long>(1);

    i = conditions.find(',');
    string one_condition = conditions.substr(0, i);
    parse_bitmask(one_condition, f.must_have[0], f.must_not_have[0]);
    if (i == conditions.npos) {
      f.stop_on_hit = true;
    } else {
      f.stop_on_hit = not(atoi(conditions.substr(i + 1).c_str()) == 0);
    }
  } else {
    f.must_have = vector<unsigned long>(2);
    f.must_not_have = vector<unsigned long>(2);

    i = conditions.find(',');
    if (i == conditions.npos) {
      cerr << "invalid filter: " << s << endl;
      exit(1);
    }
    string one_condition = conditions.substr(0, i);
    parse_bitmask(one_condition, f.must_have[0], f.must_not_have[0]);

    conditions = conditions.substr(i + 1);
    i = conditions.find(',');
    one_condition = conditions.substr(0, i);
    parse_bitmask(one_condition, f.must_have[1], f.must_not_have[1]);

    if (i == conditions.npos) {
      f.stop_on_hit = true;
    } else {
      f.stop_on_hit = not(atoi(conditions.substr(i + 1).c_str()) == 0);      
    }
  }

  if (file_map.count(f.dest_file) == 0) {
    if (f.dest_file[0] == '&') {
      int fd = atoi(f.dest_file.substr(1).c_str());
      if (fd == 1) {
	file_map[f.dest_file] = stdout;
      } else if (fd == 2) {
	file_map[f.dest_file] = stderr;
      } else {
	file_map[f.dest_file] = fdopen(fd, "w");
      }
    } else {
      file_map[f.dest_file] = fopen(f.dest_file.c_str(), "w");
    }
  }

  filter_vector.push_back(f);

  if (global::verbosity > 0) {
    clog << "added filter: " << "0x" << hex << f.must_have[0] << "/" << "0x" << hex << f.must_not_have[0];
    if (global::rg_dict.size() > 0)
      clog << "," << "0x" << hex << f.must_have[1] << "/" << "0x" << hex << f.must_not_have[1];
    clog << "," << f.stop_on_hit << ":" << f.dest_file << '\n';
  }
}


int
main(int argc, char* argv[])
{
  string progName(argv[0]);
  string pairing_file;
  vector<string> filter_list;
  cnp = &default_cnp;

  char c;
  while ((c = getopt(argc, argv, "l:N:Pf:g:v")) != -1) {
    switch (c) {
    case 'l':
      pairing_file = optarg;
      break;
    case 'N':
      num_threads = atoi(optarg);
      break;
    case 'P':
      cnp = cloneNameParser;
      fnp = fullNameParser;
      break;
    case 'f':
      filter_list.push_back(string(optarg));
      //add_filter(optarg);
      break;
    case 'g':
      global::default_rg = optarg;
      break;
    case 'v':
      global::verbosity++;
      break;
    default:
      cerr << "unrecognized option: " << c << endl;
      exit(1);
    }
  }

  if (optind + 1 < argc) {
    cerr << "use: " << argv[0] << " [options] [<mappings_sam>]" << endl;
    exit(1);
  }

  if (global::verbosity > 0) clog << "number of threads: " << num_threads << '\n';

  if (pairing_file.size() > 0) {
    igzstream pairingIn(pairing_file.c_str());
    if (!pairingIn) { cerr << "error opening pairing file: " << pairing_file << endl; exit(1); }
    load_pairing(pairingIn, global::rg_dict, global::num_rg_dict, global::rg_to_num_rg_dict);
    pairingIn.close();
  }

  for (size_t i = 0; i < filter_list.size(); ++i) {
    add_filter(filter_list[i]);
  }

  igzstream mapIn(optind < argc? argv[optind] : "-");
  if (!mapIn) {
    cerr << "error opening mappings file: " << argv[optind] << endl;
    exit(1);
  }

  SamMappingSetGen mapGen(&mapIn, cnp, addSQToRefDict, &global::refDict, true);
  pair<string,vector<SamMapping> >* m = mapGen.get_next();
  if (m != NULL) {
    map<string,stringstream*> out_str;
    //process_mapping_set(m->first, m->second, out_str, &cerr);
    process_mapping_set(m->first, m->second, out_str);
    for (map<string,stringstream*>::iterator it = out_str.begin();
	 it != out_str.end(); ++it) {
      fputs(it->second->str().c_str(), file_map[it->first]);
      delete it->second;
    }

    delete m;

    priority_queue<Chunk,vector<Chunk>,ChunkComparator> h;
    long long next_chunk_in = 0;
    long long next_chunk_out = 0;
    int chunk_size = 1000;

    //omp_lock_t input_lock;
    //omp_init_lock(&input_lock);
#pragma omp parallel num_threads(num_threads)
    {
      int tid = omp_get_thread_num();
      pair<string,vector<SamMapping> >* local_m;
      vector<pair<string,vector<SamMapping> >* > local_m_vector(chunk_size);
      int load;
      while (true) {
	Chunk chunk;
	//chunk.chunk_id = i;
	chunk.thread_id = tid;

#pragma omp critical(input)
	{
	  //omp_set_lock(&input_lock);
	  for (load = 0; load < chunk_size; ++load) {
	    local_m = mapGen.get_next();
	    if (local_m == NULL) {
	      break;
	    }
	    local_m_vector[load] = local_m;
	  }

	  chunk.chunk_id = next_chunk_in;
	  if (load > 0)
	    ++next_chunk_in;
	  //omp_unset_lock(&input_lock);
	}

	if (load == 0)
	  break;

	//chunk.out_str = new stringstream();
	chunk.err_str = new stringstream();
	*chunk.err_str << "tid=" << tid << " chunk_id=" << chunk.chunk_id
		       << " start:" << local_m_vector[0]->first
		       << " end:" << local_m_vector[load - 1]->first
		       << '\n';

	for (int i = 0; i < load; ++i) {
	  process_mapping_set(local_m_vector[i]->first, local_m_vector[i]->second,
			      chunk.out_str);
	  //#ifdef NDEBUG
	  //		    NULL
	  //#else
	  //		    chunk.err_str
	  //#endif
	  //		    );
	  delete local_m_vector[i];
	}

#pragma omp critical(output)
	{
	  h.push(chunk);
	  while (h.size() > 0) {
	    chunk = h.top();
	    assert(chunk.chunk_id >= next_chunk_out);
	    if (chunk.chunk_id > next_chunk_out) {
	      break;
	    }

	    for (map<string,stringstream*>::iterator it = chunk.out_str.begin();
		 it != chunk.out_str.end(); ++it) {
	      fputs(it->second->str().c_str(), file_map[it->first]);
	      delete it->second;
	    }

	    if (global::verbosity > 0) {
	      cerr << "chunk=" << chunk.chunk_id << " work_thread=" << chunk.thread_id
		   << " print_thread=" << tid << '\n';
	      cerr << chunk.err_str->str();
	      cerr.flush();
	    }
	    delete chunk.err_str;
	    h.pop();
	    ++next_chunk_out;
	  }
	}
      }
    }
    //omp_destroy_lock(&input_lock);

    mapIn.close();
    for (map<string,FILE*>::iterator it = file_map.begin();
	 it != file_map.end(); ++it) {
      fclose(it->second);
    }
  }

  return 0;
}
