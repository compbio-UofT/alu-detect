#include <iostream>
#include <cstdlib>
#include <vector>

#include "gzstream/gzstream.h"
#include "globals.hpp"
#include "Clone.hpp"
#include "CloneGen.hpp"
#include "SamMapping.hpp"
#include "SamMappingSetGen.hpp"
#include "Pairing.hpp"
#include "common.hpp"


using namespace std;


class TSD {
public:
  int start;
  int end;
  int count;
  TSD(int start_, int end_) : start(start_), end(end_), count(0) {}
};

class MP {
public:
  Pairing * pairing;
  int t_len;
  MP(Pairing * pairing_, int t_len_) : pairing(pairing_), t_len(t_len_) {}
};

string prog_name;
vector<TSD> tsd;
vector<vector<MP> > cluster;
string (*cnp)(const string&);
void (*fnp)(const string&, Clone&, int&);
int taboo_len = 10;
int reg_start;
int reg_end;
int total_bp;
int total_bp_left;
int total_bp_right;
int total_bp_mid;


void
process_mapping_set(const string& s, vector<SamMapping>& v)
{
  if ((global::rg_dict.size() == 0 and v.size() != 1)
      or (global::rg_dict.size() > 0 and v.size() != 2)) {
    cerr << "incorrect number of mappings for clone [" << s << "]" << endl;
    exit(1);
  }

  vector<Mapping> v2(2);
  for (size_t j = 0; j < v.size(); ++j) {
    if (v[j].mapped) v2[j] = convert_SamMapping_to_Mapping(v[j]);
    total_bp += int(v2[j].dbPos[1] - v2[j].dbPos[0] + 1);
    if (v2[j].dbPos[1] < tsd[0].start - taboo_len) 
      total_bp_left += int(v2[j].dbPos[1] - v2[j].dbPos[0] + 1);
    else if (v2[j].dbPos[0] > tsd[tsd.size() - 1].end + taboo_len)
      total_bp_right += int(v2[j].dbPos[1] - v2[j].dbPos[0] + 1);
    else if (tsd.size() == 2 and v2[j].dbPos[0] > tsd[0].end + taboo_len
	and v2[j].dbPos[1] < tsd[1].start - taboo_len)
      total_bp_mid += int(v2[j].dbPos[1] - v2[j].dbPos[0] + 1);
  }

  for (size_t i = 0; i < tsd.size(); ++i) {
    for (size_t j = 0; j < v.size(); ++j) {
      if (v[j].mapped and
	  v2[j].dbPos[0] <= tsd[i].start - taboo_len and v2[j].dbPos[1] >= tsd[i].end + taboo_len) {
	// captures this TSD!
	tsd[i].count++;
	return;
      }
    }
  }

  if (v.size() != 2 or not v[0].mapped or not v[1].mapped) return;

  Pairing * pairing = NULL;
  if (fnp == NULL) {
    // get read group info from SAM tags
    size_t k = 0;
    while (k < v[0].rest.size() and v[0].rest[k].key != "RG") ++k;
    if (k >= v[0].rest.size()) {
      cerr << "could not determine read group for clone: " << s << "\n";
      exit(EXIT_FAILURE);
    }
    RGDict::iterator it = global::rg_dict.find(v[0].rest[k].value);
    if (it == global::rg_dict.end()) {
      cerr << "error: no pairing info for RG [" << v[0].rest[k].value
	   << "] of clone [" << s << "]\n";
      exit(1);
    }
    pairing = &it->second;
  } else {
    // use full name parser to get read group info
    Clone c;
    int nip;
    fnp(v[0].name, c, nip);
    pairing = c.pairing;
  }

  // read pair, both mapped, neither straddles a TSD
  long long left_end = min(v2[0].dbPos[0], v2[1].dbPos[0]);
  long long right_end = max(v2[0].dbPos[1], v2[1].dbPos[1]);

  if (tsd.size() == 1) {
    if (left_end < tsd[0].start - taboo_len
	and right_end > tsd[0].end + taboo_len) {
      cluster[0].push_back(MP(pairing, pairing->get_t_len(v2[0], 0, v2[1], 0)));
    }
  } else {
    if (left_end < tsd[0].start - taboo_len and
	right_end > tsd[1].end + taboo_len) {
      cluster[0].push_back(MP(pairing, pairing->get_t_len(v2[0], 0, v2[1], 0)));
    } else if (left_end < tsd[0].start - taboo_len and
	       right_end > tsd[0].end + taboo_len and
	       right_end < tsd[1].start - taboo_len) {
      cluster[1].push_back(MP(pairing, pairing->get_t_len(v2[0], 0, v2[1], 0)));
    } else if (left_end > tsd[0].end + taboo_len and
	       left_end < tsd[1].start - taboo_len and
	       right_end > tsd[1].end + taboo_len) {
      cluster[2].push_back(MP(pairing, pairing->get_t_len(v2[0], 0, v2[1], 0)));
    }
  }
}

void
print_cluster_evidence(const vector<MP>& v, ostream& os)
{
  if (v.size() == 0) {
    os << ".";
  } else {
    for (size_t i = 0; i < v.size(); ++i) {
      if (i > 0) os << ";";
      os << v[i].pairing->mean << "," << v[i].pairing->stddev << "," << v[i].t_len;
    }
  }
}

string
default_cnp(const string& s)
{
  return string(s);
}

void
usage(ostream& os)
{
  os << "use: " << prog_name << " [ -l <pairing_file> ] -t <l_start>,<l_end> [ -t <r_start>,<r_end> ] [ <file> ]\n";
}

int
main(int argc, char* argv[])
{
  prog_name = string(argv[0]);
  cnp = default_cnp;

  string pairing_file;

  char c;
  while ((c = getopt(argc, argv, "l:t:s:PN:g:vh")) != -1) {
    switch (c) {
    case 'l':
      pairing_file = optarg;
      break;
    case 't':
      if (optarg[0] != '.') {
	if (tsd.size() >= 2) {
	  cerr << "wrong number of tsds\n";
	  usage(cerr);
	  exit(EXIT_FAILURE);
	}
	string tmp(optarg);
	size_t i = tmp.find(',');
	if (i == string::npos) {
	  cerr << "error parsing tsd specification: " << tmp << "\n";
	  usage(cerr);
	  exit(EXIT_FAILURE);
	}
	int start = atoi(tmp.substr(0, i).c_str());
	int end = atoi(tmp.substr(i+1).c_str());
	if (start < 0 or end < start) {
	  cerr << "error parsing tsd specification: " << tmp << "\n";
	  usage(cerr);
	  exit(EXIT_FAILURE);
	}
	if (tsd.size() == 1 and start < tsd[0].end) {
	  cerr << "tsds in wrong order\n";
	  usage(cerr);
	  exit(EXIT_FAILURE);
	}
	tsd.push_back(TSD(start, end));
      }
      break;
    case 's':
      {
	string tmp(optarg);
	size_t i = tmp.find(',');
	if (i == string::npos) {
	  cerr << "error parsing region start/end specification: " << tmp << "\n";
	  usage(cerr);
	  exit(EXIT_FAILURE);
	}
	reg_start = atoi(tmp.substr(0, i).c_str());
	reg_end = atoi(tmp.substr(i+1).c_str());
      }
      break;
    case 'P':
      cnp = cloneNameParser;
      fnp = fullNameParser;
      break;
    case 'g':
      global::default_rg = optarg;
      break;
    case 'N':
      global::num_threads = atoi(optarg);
      break;
    case 'v':
      global::verbosity++;
      break;
    case 'h':
      usage(cout);
      exit(EXIT_SUCCESS);
    default:
      cerr << "unrecognized option: " << c << endl;
      usage(cerr);
      exit(EXIT_FAILURE);
    }
  }

  if (optind + 1 < argc) {
    usage(cerr);
    exit(EXIT_FAILURE);
  }

  //if (global::verbosity > 0) clog << "number of threads: " << num_threads << '\n';
  if (tsd.size() == 1) {
    cluster = vector<vector<MP> >(1);
  } else if (tsd.size() == 2) {
    cluster = vector<vector<MP> >(3);
  } else {
    cerr << "wrong number of tsds\n";
    usage(cerr);
    exit(EXIT_FAILURE);
  }
  if (global::verbosity > 0) {
    clog << "tsd 1: [" << tsd[0].start << "," << tsd[0].end << ")\n";
    if (tsd.size() > 1) {
      clog << "tsd 2: [" << tsd[1].start << "," << tsd[1].end << ")\n";
    }
  }    

  if (pairing_file.size() > 0) {
    igzstream pairingIn(pairing_file.c_str());
    if (pairingIn.fail()) { cerr << "error opening pairing file: " << pairing_file << endl; exit(EXIT_FAILURE); }
    load_pairing(pairingIn, global::rg_dict, global::num_rg_dict, global::rg_to_num_rg_dict);
    pairingIn.close();
  }

  igzstream mapIn(optind < argc? argv[optind] : "-");
  if (mapIn.fail()) {
    cerr << "error opening mappings file: " << argv[optind] << endl;
    exit(EXIT_FAILURE);
  }

  SamMappingSetGen map_gen(&mapIn, cnp, NULL, &global::refDict, true);
  pair<string,vector<SamMapping> >* m = map_gen.get_next();
  int n_fragments = 0;
  while (m != NULL) {
    ++n_fragments;
    process_mapping_set(m->first, m->second);
    delete m;
    m = map_gen.get_next();
  }
  mapIn.close();

  if (tsd.size() == 1) {
    cout << tsd[0].count << "\t";
    print_cluster_evidence(cluster[0], cout);
  } else {
    cout << tsd[0].count << "\t";
    cout << tsd[1].count << "\t";
    print_cluster_evidence(cluster[0], cout);
    cout << "\t";
    print_cluster_evidence(cluster[1], cout);
    cout << "\t";
    print_cluster_evidence(cluster[2], cout);
  }
  cout << "\t"
       << (reg_start < tsd[0].start - taboo_len ?
	   double(total_bp_left) / double (tsd[0].start - taboo_len - reg_start + 1)
	   : 0)
       << "\t"
       << (reg_end > tsd[tsd.size() - 1].end + taboo_len ?
	   double(total_bp_right)
	   / double (reg_end - tsd[tsd.size() - 1].end - taboo_len + 1)
	   : 0);
  if (tsd.size() == 2) {
    cout << "\t"
	 << (tsd[0].end + taboo_len < tsd[1].end - taboo_len ?
	     double(total_bp_mid)
	     / double (tsd[1].end - taboo_len - tsd[0].end - taboo_len + 1)
	     : 0);
  }
  cout << "\n";

  return EXIT_SUCCESS;
}
