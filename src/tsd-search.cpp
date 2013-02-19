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
#include "common.hpp"
#include "inalign_core.hpp"
#include "deep_size.hpp"
*/
#include "DNASequence.hpp"
#include "Range.hpp"
#include "Interval.hpp"
#include "BedLine.hpp"
#include "Fasta.hpp"

#include "common/sw-vector.h"
#include "common/sw-full-common.h"
#include "common/sw-full-ls.h"


int num_threads = 1;
SQDict refDict;
vector<BedLine> bedLine;


int const match_score = 10;
int const mismatch_score = -15;
int const gap_open_score = -33;
int const gap_extend_score = -7;


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
find_tsd(BedLine * b)
{
  long long left_window_start = max(0ll, b->pos[0] - 100);
  long long left_window_len = b->pos[0] - left_window_start;

  long long right_window_start = b->pos[1] + 1;
  long long right_window_len = min(100ll, b->db->len - right_window_start);

  /*
  cerr << "processing line: " << b->line << endl
       << " chr len=" << b->db->len << endl
       << " left_window_start=" << left_window_start << endl
       << " left_window_len=" << left_window_len << endl
       << " right_window_start=" << right_window_start << endl
       << " right_window_len=" << right_window_len << endl;
  */

  if (left_window_len < 8 or right_window_len < 8)
    return;

  int thresh = 8 * match_score;
  int score = sw_vector(&b->db->seq[0].c_str()[left_window_start], 0, left_window_len,
			&b->db->seq[0].c_str()[right_window_start], right_window_len,
			NULL, -1, false);

  if (score < thresh)
    return;

  struct sw_full_results sfr;
  memset(&sfr, 0, sizeof(sfr));
  sw_full_ls(&b->db->seq[0].c_str()[left_window_start], 0, left_window_len,
	     &b->db->seq[0].c_str()[right_window_start], right_window_len,
	     thresh, score, &sfr,
	     false, NULL, 0, true);

  b->skip_left = (left_window_len - sfr.genome_start) - sfr.gmapped;
  b->len_left = sfr.gmapped;
  b->skip_right = sfr.read_start;
  b->len_right = sfr.rmapped;

  free(sfr.dbalign);
  free(sfr.qralign);
}


void
usage(const string& progName)
{
  cerr << "Use: " << progName << " [options] <fasta_file> <bed_file>" << endl;
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

  {
    string fasta_file_name;
    if (!strcmp(argv[optind], "-")) {
      fasta_file_name = "/dev/fd/0";
    } else {
      fasta_file_name = argv[optind];
    }
    cerr << "using fasta file: " << fasta_file_name << endl;
    igzstream fasta_file(fasta_file_name.c_str());
    if (!fasta_file) {
      cerr << "error opening fasta file" << endl;
      return 1;
    }
    readFasta(fasta_file, refDict);
    fasta_file.close();
    cerr << "done reading fasta file" << endl;
  }

  {
    string bed_file_name;
    if (!strcmp(argv[optind + 1], "-")) {
      bed_file_name = "/dev/fd/0";
    } else {
      bed_file_name = argv[optind + 1];
    }
    cerr << "using bed file: " << bed_file_name << endl;
    igzstream bed_file(bed_file_name.c_str());
    if (!bed_file) {
      cerr << "error opening bed file" << endl;
      return 1;
    }
    while (true) {
      string bed_line;
      getline(bed_file, bed_line);
      if (bed_file.eof())
	break;
      bedLine.push_back(BedLine(bed_line, &refDict));
    }
    bed_file.close();
  }

#pragma omp parallel num_threads(num_threads)
  {
    sw_vector_setup(200, 200,
		    gap_open_score, gap_extend_score,
		    gap_open_score, gap_extend_score,
		    match_score, mismatch_score,
		    false, true);
    sw_full_ls_setup(200, 200,
		     gap_open_score, gap_extend_score,
		     gap_open_score, gap_extend_score,
		     match_score, mismatch_score,
		     true, 8);

#pragma omp for
    for (int i = 0; i < (int)bedLine.size(); ++i) {
      find_tsd(&bedLine[i]);
    }

    sw_vector_cleanup();
    sw_full_ls_cleanup();
  }

  for (vector<BedLine>::iterator it = bedLine.begin(); it != bedLine.end(); ++it) {
    cout << it->line << "\t" << it->skip_left
	 << "\t" << it->len_left
	 << "\t" << it->skip_right
	 << "\t" << it->len_right
	 << endl;
  }

  return 0;
}
