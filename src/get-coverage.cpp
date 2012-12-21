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
/*
#include "Clone.hpp"
#include "Read.hpp"
#include "Mapping.hpp"
#include "SamMappingSetGen.hpp"
#include "CloneGen.hpp"
#include "Fasta.hpp"
#include "common.hpp"
#include "inalign_core.hpp"
#include "deep_size.hpp"
*/
#include "SamMapping.hpp"
#include "DNASequence.hpp"
#include "Range.hpp"
#include "Interval.hpp"
#include "Cigar.hpp"


int num_threads = 1;
SQDict refDict;
vector<vector<int> > reg;
int reg_size = 100;
vector<double> min_covg;
bool output_bed = false;


void
addSQDict(const string& line, SQDict& dict)
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
      //cerr << "error: missing sequence for contig [" << contig_name << "]";
      //exit(1);
      contig->name = contig_name;
      contig->len = contig_len;
      contig->idx = dict.size() - 1;
      cerr << "found contig [" << contig->name << "] with length [" << contig->len << "]" << endl;
    } else if (contig->len != contig_len) {
      cerr << "error: contig [" << contig->name << "] has different lengths: "
	   << contig->len << " or " << contig_len << endl;
      exit(1);
    }
  }
}


void
process_mapping(const string& line)
{
  SamMapping m(line, &refDict, false);
  //cerr << "processing mapping for read: " << m.name << endl;
  if (m.flags[2])
    return;

  Interval<int> qrPos;
  Interval<long long> dbPos;
  parseCigar(m.cigar, m.dbPos, qrPos, dbPos);

  long long db_start = dbPos[0];
  int db_len = dbPos[1] - dbPos[0];
  long long crt_reg = db_start / reg_size;

  while (db_len > 0) {
    if (db_start + db_len <= (crt_reg + 1) * reg_size) {

#pragma omp atomic
      reg[m.db->idx][crt_reg] += db_len;

      break;
    } else {
      int len_in_crt_reg = int((crt_reg + 1) * reg_size - db_start);

#pragma omp atomic
      reg[m.db->idx][crt_reg] += len_in_crt_reg;

      db_start = (crt_reg + 1) * reg_size;
      db_len -= len_in_crt_reg;
      ++crt_reg;
    }
  }
}


void
usage(const string& progName)
{
  cerr << "Use: " << progName << " [options]" << endl;
  cerr << "Options:" << endl
       << "  -N <num_threads>: Use the specifed number of threads" << endl
       << "  -r <range>: Only process mappings from given range (not currently working)" << endl
       << "  -s <reg_size>: Use the given region size (currently: " << reg_size << ")" << endl
       << "  -m <min_covg>: Minimum coverage per region (can be given multiple times)" << endl
       << "  -b: Output bed file of coverage for the first minimum coverage threshold" << endl;
}


int
main(int argc, char* argv[])
{
  string progName(argv[0]);
  vector<Range> range;

  global::min_tail_len = 1;

  char c;
  while ((c = getopt(argc, argv, "N:r:s:m:b")) != -1) {
    switch (c) {
    case 'N':
      num_threads = atoi(optarg);
      break;
    case 'r':
      range.push_back(Range(optarg));
      break;
    case 's':
      reg_size = atoi(optarg);
      break;
    case 'm':
      min_covg.push_back(atof(optarg));
      break;
    case 'b':
      output_bed = true;
      break;
    default:
      cerr << "invalid option [" << c << "]" << endl;
      usage(progName);
      return 1;
    }
  }
  if (min_covg.size() == 0) {
    min_covg.push_back(1.0);
  }

  cerr.precision(2);
  cerr << "using " << num_threads << " threads" << endl;
  cerr << "using regions of size: " << reg_size << endl;
  cerr << fixed << "using minimum coverage:";
  for (size_t i = 0; i < min_covg.size(); ++i) {
    cerr << ' ' << min_covg[i];
  }
  cerr << endl;

  string file_name;
  if (optind == argc || !strcmp(argv[optind], "-")) {
    cerr << "using stdin" << endl;
    file_name = "/dev/fd/0";
  } else {
    cerr << "using input file: " << argv[optind] << endl;
    file_name = argv[optind];
  }
  igzstream f(file_name.c_str());

  string line;
  getline(f, line);
  while (!f.eof() and line[0] == '@') {
    addSQDict(line, refDict);
    getline(f, line);
  }
  if (f.eof()) {
    cerr << "no mappings in input file" << endl;
    cout.precision(2);
    for (size_t k = 0; k < min_covg.size(); ++k) {
      cout << fixed << min_covg[k] << "\t0\t0" << endl;
    }
    return 0;
  }

  // parse range restrictions
  for_iterable(vector<Range>, range, it) {
    it->parse(refDict);
  }

  // set up region vectors
  reg = vector<vector<int> >(refDict.size());
  for_iterable(SQDict, refDict, it) {
    long long n_regs = (it->second.len - 1) / reg_size + 1;
    reg[it->second.idx] = vector<int>(n_regs, 0);
  }
  process_mapping(line);

  omp_lock_t input_lock;
  omp_init_lock(&input_lock);
#pragma omp parallel num_threads(num_threads)
  {
    string local_line;
    while (true) {
      //#pragma omp critical(input)
      {
	omp_set_lock(&input_lock);
	getline(f, local_line);
	omp_unset_lock(&input_lock);
      }

      if (f.eof())
	break;

      process_mapping(local_line);
    }
  }
  omp_destroy_lock(&input_lock);

  f.close();
  cerr << "done reading mappings" << endl;

  if (!output_bed) {
    for (size_t k = 0; k < min_covg.size(); ++k) {
      long long n_regs = 0;
      long long n_regs_min_coverage = 0;
      long long bp_total = 0;
      long long bp_regs_min_coverage = 0;
      for_iterable(SQDict, refDict, it) {
	long long regs = (it->second.len - 1) / reg_size + 1;
	n_regs += regs;
	for (long long i = 0; i < regs; ++i) {
	  bp_total += reg[it->second.idx][i];
	  if (reg[it->second.idx][i] > (long long)(min_covg[k] * reg_size)) {
	    ++n_regs_min_coverage;
	    bp_regs_min_coverage += reg[it->second.idx][i];
	  }
	}
      }

      cerr << "min_covg: " << min_covg[k] << endl;
      cerr << "n_regs: " << n_regs << endl;
      cerr << "n_regs_min_coverage: " << n_regs_min_coverage << endl;
      cerr << "bp_total: " << bp_total << endl;
      cerr << "bp_regs_min_coverage: " << bp_regs_min_coverage << endl;

      cout.precision(2);
      cout << fixed << min_covg[k]
	   << '\t' << (double(n_regs_min_coverage) / n_regs) * 100.0 << '\t'
	   << (n_regs_min_coverage > 0?
	       double(bp_regs_min_coverage) / (n_regs_min_coverage * reg_size)
	       : 0.0)
	   << endl;
    }
  } else { // output bed intervals with minimum coverage
    for_iterable(SQDict, refDict, it) {
      long long regs = (it->second.len - 1) / reg_size + 1;
      long long last_start = -1;
      bool last_active = false;
      for (long long i = 0; i < regs; ++i) {
	if (reg[it->second.idx][i] > (long long)(min_covg[0] * reg_size)) {
	  if (not last_active) {
	    last_start = i;
	    last_active = true;
	  }
	} else {
	  if (last_active) {
	    cout << it->second.name << '\t' << last_start * reg_size << '\t' << i * reg_size << endl;
	    last_active = false;
	  }
	}
      }
      if (last_active) {
	cout << it->second.name << '\t' << last_start * reg_size << '\t' << it->second.len << endl;
      }
    }
  }

  return 0;
}
