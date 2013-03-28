#include <iostream>
#include <cstdlib>
#include <string>

using namespace std;

string delim("\t");
bool print_newline = true;

int
main(int argc, char * argv[])
{
  char c;
  while ((c = getopt(argc, argv, "t:n")) != -1) {
    switch (c) {
    case 't':
      delim = optarg;
      break;
    case 'n':
      print_newline = false;
      break;
    default:
      cerr << "unrecognized option: " << c << endl;
      exit(EXIT_FAILURE);
    }
  }
  
  if (argc > optind) {
    cout << argv[optind];
    for (int i = optind + 1; i < argc; ++i)
      cout << delim << argv[i];
  }
  if (print_newline)
    cout << '\n';
  return EXIT_SUCCESS;
}
