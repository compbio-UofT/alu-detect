using namespace std;

#include <cstdlib>
#include <iostream>
#include <string>
#include <map>
#include <vector>
#include <omp.h>

#include "gzstream/gzstream.h"
#include "strtk/strtk.hpp"
/*
#include "globals.hpp"
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
#include "BedLine.hpp"


int num_threads = 1;
SQDict refDict;
vector<vector<int> > reg;
int reg_size = 100;
vector<double> min_covg;
bool output_bed = false;

typedef pair<Contig*, long long> chrpos;

bool
chrpos_comparator(const chrpos& lhs, const chrpos& rhs)
{
  return (lhs.first->idx < rhs.first->idx
	  || (lhs.first->idx == rhs.first->idx && lhs.second <= rhs.second));
}

map <chrpos, BedLine, bool(*)(const chrpos&, const chrpos&)> bedLine_map(chrpos_comparator);


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
  SamMapping m(line, &refDict);
  //cerr << "processing mapping for read: " << m.name << endl;
  if (m.flags[2])
    return;

  Interval<int> qrPos;
  Interval<long long> dbPos;
  parseCigar(m.cigar, m.dbPos, qrPos, dbPos);

  //cerr << "found mapping at " << m.db->name << ":" << dbPos[0] << "-" << dbPos[1] << endl;

  chrpos tmp_chrpos(m.db, dbPos[1]);
  map<chrpos, BedLine>::iterator it;
  it = bedLine_map.lower_bound(tmp_chrpos);

  while (true) {
    if (it == bedLine_map.begin())
      break;
    --it;

    //cerr << "candidate bed interval: " << it->second.line << endl;
    if (it->second.db->idx != m.db->idx or dbPos[0] > it->second.pos[1])
      break;

    // mapping intersects the bed interval
    int bp = min(dbPos[1], it->second.pos[1]) - max(dbPos[0], it->second.pos[0]) + 1;
    //cerr << "intersecting for " << bp << "bp" << endl;
    assert(bp >= 1);
#pragma omp atomic
    it->second.support += bp;
  }
}


void
usage(const string& progName)
{
  cerr << "Use: " << progName << " [options] <sam_file> <bed_file>" << endl;
  cerr << "Options:" << endl
       << "  -N <num_threads>: Use the specifed number of threads" << endl;
  //<< "  -r <range>: Only process mappings from given range (not currently working)" << endl
  //<< "  -s <reg_size>: Use the given region size (currently: " << reg_size << ")" << endl
  //<< "  -m <min_covg>: Minimum coverage per region (can be given multiple times)" << endl
  //<< "  -b: Output bed file of coverage for the first minimum coverage threshold" << endl;
}


int
main(int argc, char* argv[])
{
  string progName(argv[0]);
  //vector<Range> range;

  char c;
  while ((c = getopt(argc, argv, "N:")) != -1) {
    switch (c) {
    case 'N':
      num_threads = atoi(optarg);
      break;
      /*
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
      */
    default:
      cerr << "invalid option [" << c << "]" << endl;
      usage(progName);
      return 1;
    }
  }
  if (optind + 2 != argc) {
    usage(progName);
    return 1;
  }

  cerr << "using " << num_threads << " threads" << endl;

  string sam_file_name;
  if (!strcmp(argv[optind], "-")) {
    cerr << "sam mappings: stdin" << endl;
    sam_file_name = "/dev/fd/0";
  } else {
    cerr << "sam mappings: " << argv[optind] << endl;
    sam_file_name = argv[optind];
  }
  igzstream sam_file(sam_file_name.c_str());

  string bed_file_name;
  if (!strcmp(argv[optind + 1], "-")) {
    cerr << "bed file: stdin" << endl;
    bed_file_name = "/dev/fd/0";
  } else {
    cerr << "bed file: " << argv[optind + 1] << endl;
    bed_file_name = argv[optind + 1];
  }
  igzstream bed_file(bed_file_name.c_str());

  string sam_line;
  getline(sam_file, sam_line);
  while (!sam_file.eof() and sam_line[0] == '@') {
    addSQDict(sam_line, refDict);
    getline(sam_file, sam_line);
  }
  if (sam_file.eof()) {
    cerr << "no mappings in input file" << endl;
  } else {

    // load bed intervals
    while (true) {
      string bed_line;
      getline(bed_file, bed_line);

      if (bed_file.eof())
	break;

      BedLine tmp_bedLine(bed_line, &refDict);
      chrpos tmp_chrpos(tmp_bedLine.db, tmp_bedLine.pos[0]);
      bedLine_map[tmp_chrpos] = tmp_bedLine;
      //cerr << "added interval with key: " << tmp_chrpos.first->name << ":" << tmp_chrpos.second << endl;
    }

    process_mapping(sam_line);

    omp_lock_t input_lock;
    omp_init_lock(&input_lock);
#pragma omp parallel num_threads(num_threads)
    {
      string local_sam_line;
      while (true) {
	//#pragma omp critical(input)
	{
	  omp_set_lock(&input_lock);
	  getline(sam_file, local_sam_line);
	  omp_unset_lock(&input_lock);
	}

	if (sam_file.eof())
	  break;

	process_mapping(local_sam_line);
      }
    }
    omp_destroy_lock(&input_lock);

    cerr << "done reading mappings" << endl;
  }
  sam_file.close();
  bed_file.close();

  map<chrpos, BedLine>::iterator it;
  for (it = bedLine_map.begin(); it != bedLine_map.end(); ++it) {
    cout << it->second.line << "\t" << it->second.support << endl;
  }

  return 0;
}
