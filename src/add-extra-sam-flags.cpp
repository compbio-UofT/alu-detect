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
#include "CloneGen.hpp"
#include "Cigar.hpp"



int num_threads = 1;

int min_read_len = 20;
int min_mqv = 5;
int min_tail_insert_size = 15;
int min_tail_match_len = 5;

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
addSQToRefDict_then_print(const string& line)
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
      cerr << "added contig [" << contig->name << "] of length [" << contig->len << "]" << endl;
    } else {
      cerr << "error: contig [" << contig->name << "] already exists!" << endl;
      exit(1);
    }
  }
  cout << line << endl;
}


void
process_mapping_set(const string& s, vector<SamMapping>& v,
		    ostream& out_str, ostream& err_str)
{
  if ((!global::pairing.paired && v.size() != 1)
      or (global::pairing.paired && v.size() != 2)) {
    cerr << "incorrect number of mappings for clone [" << s << "]" << endl;
    exit(1);
  }

  for (size_t i = 0; i < v.size(); ++i) {
    int len = (v[i].seq.compare("*")? v[i].seq.size() : 0);
    /*
    err_str << "clone s=" << s << " i=" << i << " len=" << len
	    << " v[i].seq=" << v[i].seq << " v[i].seq.size()=" << v[i].seq.size()
	    << " v[i].seq.compare(\"*\")=" << v[i].seq.compare("*") << endl;
    */
    if (len < min_read_len) {
      v[i].flags[16] = 1;
    }

    if (v[i].flags[2] == 0) { // mapped
      if (v[i].mqv >= min_mqv) {
	v[i].flags[12] = 1;
      }
      if (min_tail_insert_size > 0) {
	vector<int> tails(2);
	get_tail_insert_size(v[i].cigar, min_tail_match_len, tails);
	if (tails[0] >= min_tail_insert_size) {
	  v[i].flags[13] = 1;
	}
	if (tails[1] >= min_tail_insert_size) {
	  v[i].flags[14] = 1;
	}
      }
    }
  }

  if (global::pairing.paired and v[0].flags[2] == 0 and v[1].flags[2] == 0) {
    Clone c;
    for (int i = 0; i < 2; ++i) {
      int nip = v[i].flags[7];
      if (fnp != NULL) {
	fnp(v[i].name, c, nip);
      } else {
	c.read[nip].len = (v[i].seq.compare("*")? v[i].seq.size() : 0);
      }
      Mapping m = convert_SamMapping_to_Mapping(v[i]);
      m.qr = &c.read[nip];
      m.is_ref = true;
      c.read[nip].mapping.push_back(m);
    }

    if (global::pairing.pair_concordant(c.read[0].mapping[0], 0, c.read[1].mapping[0], 0)) {
      v[0].flags[15] = 1;
      v[1].flags[15] = 1;
      err_str << "clone s=" << s << ": concordant" << endl;
    } else {
      err_str << "clone s=" << s << ": discordant" << endl;
    }
  }

  for (size_t i = 0; i < v.size(); ++i) {
    out_str << v[i].name
	    << '\t' << v[i].flags.to_ulong()
	    << '\t' << (v[i].db != NULL? v[i].db->name : "*")
	    << '\t' << v[i].dbPos
	    << '\t' << v[i].mqv
	    << '\t' << v[i].cigar
	    << '\t' << (v[i].mp_db != NULL? (v[i].mp_db == v[i].db? "=" : v[i].mp_db->name) : "*")
	    << '\t' << v[i].mp_dbPos
	    << '\t' << v[i].tLen
	    << '\t' << v[i].seq
	    << '\t' << v[i].qvString;
    for (size_t j = 0; j < v[i].rest.size(); ++j) {
      out_str << '\t' << v[i].rest[j].key << ':' << v[i].rest[j].type << ':' << v[i].rest[j].value;
    }
    out_str << endl;
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
  cnp = &default_cnp;

  char c;
  while ((c = getopt(argc, argv, "p:N:Pq:i:")) != -1) {
    switch (c) {
    case 'p':
      global::pairing = Pairing(string(optarg));
      cerr << "set pairing: " << global::pairing << endl;
      break;
    case 'N':
      num_threads = atoi(optarg);
      break;
    case 'P':
      cnp = cloneNameParser;
      fnp = fullNameParser;
      break;
    case 'q':
      min_mqv = atoi(optarg);
      break;
    case 'i':
      min_tail_insert_size = atoi(optarg);
      break;
    default:
      cerr << "unrecognized option: " << c << endl;
      exit(1);
    }
  }

  if (optind + 1 < argc) {
    cerr << "use: " << argv[0]
	 << " [options] [<mappings_sam>]" << endl;
    exit(1);
  }

  cerr << "number of threads: " << num_threads << endl;

  igzstream mapIn(optind < argc? argv[optind] : "-");
  if (mapIn.bad()) {
    cerr << "error opening mappings file: " << argv[optind] << endl;
    exit(1);
  }

  SamMappingSetGen mapGen(&mapIn, cnp, addSQToRefDict_then_print, &global::refDict);
  pair<string,vector<SamMapping> >* m = mapGen.get_next();
  process_mapping_set(m->first, m->second, cout, cerr);
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
    vector<pair<string,vector<SamMapping> >* > local_m_vector;
    int i;
    while (true) {
      Chunk chunk;
      //chunk.chunk_id = i;
      chunk.thread_id = tid;

#pragma omp critical(input)
      {
        //omp_set_lock(&input_lock);
	for (i = 0; i < chunk_size; ++i) {
	  local_m = mapGen.get_next();
	  if (local_m == NULL) {
	    break;
	  }
	  local_m_vector.push_back(local_m);
	}

	chunk.chunk_id = next_chunk_in;
	if (i > 0)
	  ++next_chunk_in;
        //omp_unset_lock(&input_lock);
      }

      if (local_m_vector.size() == 0)
        break;

      chunk.out_str = new stringstream();
      chunk.err_str = new stringstream();
      *chunk.err_str << "tid=" << tid << " chunk_id=" << chunk.chunk_id
		     << " start:" << local_m_vector[0]->first
		     << " end:" << local_m_vector[i-1]->first
		     << endl;

      for (i = 0; i < (int)local_m_vector.size(); ++i) {
	process_mapping_set(local_m_vector[i]->first, local_m_vector[i]->second,
			    *chunk.out_str, *chunk.err_str);
	delete local_m_vector[i];
      }
      local_m_vector.clear();

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
	  cout.flush();
	  cerr << "chunk=" << chunk.chunk_id << " work_thread=" << chunk.thread_id
	       << " print_thread=" << tid << endl;
	  cerr << chunk.err_str->str();
	  cerr.flush();
	  delete chunk.out_str;
	  delete chunk.err_str;
	  h.pop();
	  ++next_chunk_out;
	}
      }
    }
  }
  //omp_destroy_lock(&input_lock);

  mapIn.close();

  return 0;
}
