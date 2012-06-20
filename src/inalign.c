#include <cstdlib>
#include <getopt.h>
#include <iostream>

#include "gzstream/gzstream.h"
#include "inalign_core.hpp"
#include "Clone.hpp"
#include "DNASequence.hpp"
#include "Fasta.hpp"
#include "common.hpp"

#define MIN_OVERLAP_START 10
#define MIN_OVERLAP_END 10


char const short_options[] = "rs:e:";
struct option const long_options[] =
  {
    {"report-mappings", 0, 0, 'r'},
    {"overlap-start", 1, 0, 's'},
    {"overlap-end", 1, 0, 'e'}
  };


int
main(int argc, char* argv[])
{
  char const * const prog_name = argv[0];

  SQDict refDict;
  long long int overlapStart = -1;
  long long int overlapEnd = -1;
  int relOverlapStart;
  int relOverlapEnd;
  SQDict repDict;
  vector<pair<Contig*,int> > repeat;
  vector<Clone> clone;
  int update_mappings = 0;

  char c;
  while ((c = getopt_long(argc, argv, short_options, long_options, NULL)) != -1) {
    switch (c) {
    case 'r':
      update_mappings = 1;
      cerr << "reporting mappings only" << endl;
      break;
    case 's':
      overlapStart = atoll(optarg);
      break;
    case 'e':
      overlapEnd = atoll(optarg);
      break;
    default:
      cerr << "usage: " << prog_name << " [options] <reference.fa> <repeats.fa> <reads.fq>" << endl;
      exit(1);
      break;
    }
  }
  argc -= optind - 1;
  argv += optind - 1;

  if (argc != 4) {
    cerr << "usage: " << prog_name << " [options] <reference.fa> <repeats.fa> <reads.fq>" << endl;
    exit(1);
  }

  // reference file
  igzstream refFaIn(argv[1]);
  if (refFaIn.bad()) {
    cerr << "error opening reference fasta file: " << argv[1] << endl;
    exit(1);
  }
  readFasta(refFaIn, refDict, true);
  refFaIn.close();
  if (refDict.size() == 0) {
    cerr << "error: did not find reference sequence" << endl;
    exit(1);
  } else if (refDict.size() >= 2) {
    cerr << "error: found more than 1 reference sequence" << endl;
    exit(1);
  }
  Contig& refContig = refDict.begin()->second;

  // adjust overlap start and end
  if (overlapStart < 0) {
    relOverlapStart = MIN_OVERLAP_START;
  } else {
    relOverlapStart = int(overlapStart - 1 - refContig.seqOffset[0] + 1);
    if (relOverlapStart < MIN_OVERLAP_START) {
      cerr << "overlap start too small [" << relOverlapStart
	   << "]; adjusting to minimum [" << MIN_OVERLAP_START << "]" << endl;
      relOverlapStart = MIN_OVERLAP_START;
    }
  }
  overlapStart = refContig.seqOffset[0] + relOverlapStart - 1;
  if (overlapEnd < 0) {
    relOverlapEnd = int(refContig.len - MIN_OVERLAP_END + 1);
  } else {
    relOverlapEnd = int(overlapEnd - 1 - refContig.seqOffset[0] + 1);
    if (relOverlapEnd > int(refContig.len - MIN_OVERLAP_END + 1)) {
      cerr << "overlap end too large [" << relOverlapEnd
	   << "]; adjusting to maximum [" << refContig.len - MIN_OVERLAP_END + 1 << "]" << endl;
      relOverlapStart = int(refContig.len - MIN_OVERLAP_END + 1);
    }
  }
  overlapEnd = refContig.seqOffset[0] + relOverlapEnd - 1;
  if (overlapEnd <= overlapStart) {
    cerr << "overlap start and end in wrong order" << endl;
    exit(1);
  }

  // repeats file
  igzstream repFaIn(argv[2]);
  if (repFaIn.bad()) {
    cerr << "error opening repeat fasta file: " << argv[2] << endl;
    exit(1);
  }
  readFasta(repFaIn, repDict);
  repFaIn.close();
  addRCToDict(repDict);

  // create list of repeats to investigate
  for (SQDict::iterator it = repDict.begin(); it != repDict.end(); ++it) {
    repeat.push_back(pair<Contig*,int>(&it->second, 0));
  }
  for (SQDict::iterator it = repDict.begin(); it != repDict.end(); ++it) {
    repeat.push_back(pair<Contig*,int>(&it->second, 1));
  }

  // clones file
  igzstream cloneFqIn(argv[3]);
  if (cloneFqIn.bad()) {
    cerr << "error opening clone fastq file: " << argv[3] << endl;
    exit(1);
  }
  clone = readAllFastq(cloneFqIn, cloneNameParser, fullNameParser);
  cloneFqIn.close();

  // make a vector of pointers
  vector<Clone*> clone_p;
  for (unsigned int i = 0; i < clone.size(); i++) {
    clone_p.push_back(&clone[i]);
  }

  run(refContig,
      refContig.seqOffset[0], refContig.seqOffset[0] + refContig.len - 1,
      overlapStart, overlapEnd,
      repeat,
      clone_p,
      update_mappings,
      cout, cerr);

  if (update_mappings) {
    for (vector<Clone>::iterator it = clone.begin(); it != clone.end(); ++it) {
      for (int nip = 0; nip < 2; ++nip) {
	cout << it->read[nip].name;
	if (it->read[nip].mapping.size() == 0) {
	  cout << "\t*";
	} else {
	  for (vector<Mapping>::iterator it2 = it->read[nip].mapping.begin();
	       it2 != it->read[nip].mapping.end(); ++it2) {
	    cout << "\t" << it2->db->name;
	    if (!it2->is_ref) {
	      cout << (it2->st == 0? "+" : "-");
	    }
	    cout << "\t" << it2->dbPos[0] + 1 << "\t" << it2->dbPos[1] + 1;
	  }
	}
	cout << endl;
      }
    }
  }

  return 0;
}
