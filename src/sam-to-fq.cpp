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


bool save_rgid = false;
string (*cnp)(const string&);
void (*fnp)(const string&, Clone&, int&);


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
	clog << "added contig [" << contig->name << "] of length [" << contig->len << "]" << '\n';
    } else {
      cerr << "error: contig [" << contig->name << "] already exists!" << endl;
      exit(1);
    }
  }
  //cout << line << endl;
}


void
print_read_from_mapping(const string& s, const SamMapping& m, ostream* os)
{
  string tail;
  if (((m.flags.to_ulong() & 0x1) != 0) and m.name == s) {
    if ((m.flags.to_ulong() & 0x40) != 0) {
      tail = "/1";
    } else {
      tail = "/2";
    }
  }
  string seq;
  string qual;
  if (m.cigar.find('H') == string::npos and m.seq != "*" and m.seq.size() == m.qvString.size()) {
    if ((m.flags.to_ulong() & 0x10) != 0) {
      seq = reverseComplement(m.seq);
      qual = reverse(m.qvString);
    } else {
      seq = m.seq;
      qual = m.qvString;
    }
  } else {
    seq = "*";
    qual = "*";
  }
  string rgid;
  if (save_rgid) {
    size_t i = 0;
    while (i < m.rest.size() and m.rest[i].key != "RG") ++i;
    if (i < m.rest.size()) {
      if (global::rg_to_num_rg_dict.count(m.rest[i].value) != 1) {
	cerr << "RG " << m.rest[i].value << " found in dictionary "
	     << global::rg_to_num_rg_dict.count(m.rest[i].value) << " times" << endl;
	exit(1);
      }
      rgid = global::rg_to_num_rg_dict[m.rest[i].value];
    } else {
      if (global::rg_to_num_rg_dict.count(global::default_rg) != 1) {
	cerr << "RG " << global::default_rg << " found in dictionary "
	     << global::rg_to_num_rg_dict.count(global::default_rg) << " times" << endl;
	exit(1);
      }
      rgid = global::rg_to_num_rg_dict[global::default_rg];
    }
  }
  *os << '@' << m.name << tail << '\n'
      << seq << '\n'
      << '+' << rgid << '\n'
      << qual << '\n';
}


void
process_mapping_set(const string& s, vector<SamMapping>& v, ostream* out_str)
{
  /*
  if ((global::rg_dict.size() == 0 and v.size() != 1)
      or (global::rg_dict.size() > 0 and v.size() != 2)) {
    cerr << "incorrect number of mappings for clone [" << s << "]" << endl;
    exit(1);
  }
  */

  for (size_t i = 0; i < v.size(); ++i) {
    print_read_from_mapping(s, v[i], out_str);
  }
}


string
default_cnp(const string& s)
{
  return string(s);
}


int
main(int argc, char* argv[])
{
  string progName(argv[0]);
  string pairing_file;
  cnp = &default_cnp;

  char c;
  while ((c = getopt(argc, argv, "l:N:Pg:vs")) != -1) {
    switch (c) {
    case 'l':
      pairing_file = optarg;
      break;
    case 'N':
      global::num_threads = atoi(optarg);
      break;
    case 'P':
      cnp = cloneNameParser;
      fnp = fullNameParser;
      break;
    case 'g':
      global::default_rg = optarg;
      break;
    case 'v':
      global::verbosity++;
      break;
    case 's':
      save_rgid = true;
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

  if (global::verbosity > 0) clog << "number of threads: " << global::num_threads << '\n';

  if (pairing_file.size() > 0) {
    igzstream pairingIn(pairing_file.c_str());
    if (!pairingIn) { cerr << "error opening pairing file: " << pairing_file << endl; exit(1); }
    load_pairing(pairingIn, global::rg_dict, global::num_rg_dict, global::rg_to_num_rg_dict);
    pairingIn.close();
  }

  igzstream mapIn(optind < argc? argv[optind] : "-");
  if (!mapIn) {
    cerr << "error opening mappings file: " << argv[optind] << endl;
    exit(1);
  }

  SamMappingSetGen mapGen(&mapIn, cnp, addSQToRefDict, &global::refDict, true);
  pair<string,vector<SamMapping> >* m = mapGen.get_next();
  if (m != NULL) {
    process_mapping_set(m->first, m->second, &cout);
    delete m;

    priority_queue<Chunk,vector<Chunk>,ChunkComparator> h;
    long long next_chunk_in = 0;
    long long next_chunk_out = 0;
    int chunk_size = 1000;

    //omp_lock_t input_lock;
    //omp_init_lock(&input_lock);
#pragma omp parallel num_threads(global::num_threads)
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

	chunk.out_str = new stringstream();
	chunk.err_str = new stringstream();
	*chunk.err_str << "tid=" << tid << " chunk_id=" << chunk.chunk_id
		       << " start:" << local_m_vector[0]->first
		       << " end:" << local_m_vector[load - 1]->first
		       << '\n';

	for (int i = 0; i < load; ++i) {
	  process_mapping_set(local_m_vector[i]->first, local_m_vector[i]->second,
			      chunk.out_str);
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

	    cout << chunk.out_str->str();

	    if (global::verbosity > 0) {
	      cerr << "chunk=" << chunk.chunk_id << " work_thread=" << chunk.thread_id
		   << " print_thread=" << tid << '\n';
	      cerr << chunk.err_str->str();
	      cerr.flush();
	    }
	    delete chunk.out_str;
	    delete chunk.err_str;
	    h.pop();
	    ++next_chunk_out;
	  }
	}
      }
    }
    //omp_destroy_lock(&input_lock);
  }

  mapIn.close();

  return 0;
}
