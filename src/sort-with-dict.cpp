#include <iostream>
#include <fstream>
#include <cstdlib>
#include <map>
#include <queue>
#include <string>

using namespace std;


string prog_name;
int verbosity;
char delimiter = '\t';
int field_number = 1;

void
usage(ostream& os)
{
  os << "use: " << prog_name
     << " [ -t <delim> ] [ -k <field_num> ] <dict_file> <main_file>\n";
}

int
main(int argc, char * argv[])
{
  prog_name = string(argv[0]);

  char c;
  while ((c = getopt(argc, argv, "t:k:vh")) != -1) {
    switch (c) {
    case 't':
      delimiter = optarg[0];
      break;
    case 'k':
      field_number = atoi(optarg);
      if (field_number <= 0) {
	cerr << "field_number must be positive: " << optarg << "\n";
	exit(EXIT_FAILURE);
      }
      break;
    case 'v':
      verbosity++;
      break;
    case 'h':
      usage(cout);
      exit(EXIT_SUCCESS);
      break;
    default:
      cerr << "unrecognized option: " << c << "\n";
      usage(cerr);
      exit(EXIT_FAILURE);
    }
  }

  if (argc - optind != 2) {
    cerr << "error: expecting exactly 2 arguments\n";
    usage(cerr);
    exit(EXIT_FAILURE);
  }

  ifstream dict_file(argv[optind]);
  if (not dict_file.is_open()) {
    cerr << "error: could not open file: " << argv[optind] << "\n";
    exit(EXIT_FAILURE);
  }

  ifstream main_file(argv[optind + 1]);
  if (not main_file.is_open()) {
    cerr << "error: could not open file: " << argv[optind + 1] << "\n";
    exit(EXIT_FAILURE);
  }

  if (verbosity > 0) {
    clog << "using delimiter: [" << delimiter << "]\n";
    clog << "using field number: [" << field_number << "]\n";
    clog << "using dict_file: [" << argv[optind] << "]\n";
    clog << "using main_file: [" << argv[optind + 1] << "]\n";
  }

  map<string, queue<string> > buffer;
  map<string, queue<string> >::iterator it;

  while (true) {
    string dict_key;
    getline(dict_file, dict_key);
    if (dict_file.eof())
      break;
    if (verbosity >= 2) {
      clog << "looking for key: [" << dict_key << "]\n";
    }
    it = buffer.find(dict_key);
    if (it == buffer.end()) {
      while (true) {
	string line;
	getline(main_file, line);
	if (main_file.eof()) {
	  cerr << "EOF reached while looking for dict_key: " << dict_key << "\n";
	  exit(EXIT_FAILURE);
	}
	size_t k = 0;
	int i = 0;
	while (k != string::npos and i < field_number - 1) {
	  k = line.find_first_of(delimiter, k);
	  //cerr << "find_first_of: " << k << "\n";
	  if (k != string::npos) k++;
	  i++;
	}
	if (k == string::npos) {
	  cerr << "line does not contain " << field_number << " fields separated by [" << delimiter << "]\n";
	  cerr << line << "\n";
	  exit(EXIT_FAILURE);
	}
	size_t l = line.find_first_of(delimiter, k);
	string crt_key = line.substr(k, l != string::npos? l-k : string::npos);
	if (verbosity >= 2) {
	  clog << "found key: [" << crt_key << "]\n";
	}
	buffer[crt_key].push(line);
	if (crt_key == dict_key) {
	  it = buffer.find(dict_key);
	  break;
	}
      }
    }
    cout << it->second.front() << "\n";
    it->second.pop();
    if (it->second.size() == 0)
      buffer.erase(it);
  }

  return EXIT_SUCCESS;
}
