#include <iostream>
#include <cstdlib>
#include <string>

std::string delim("\t");
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
      std::cerr << "unrecognized option: " << c << std::endl;
      exit(EXIT_FAILURE);
    }
  }
  
  if (argc > optind) {
    std::cout << argv[optind];

    for (int i = optind + 1; i < argc; ++i)
      std::cout << delim << argv[i];

    if (print_newline)
      std::cout << std::endl;
  }
  return EXIT_SUCCESS;
}
