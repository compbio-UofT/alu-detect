#include <iostream>
#include <fstream>
#include <cstdlib>
#include <csignal>
#include <cstring>
#include <cerrno>
#include <vector>

using namespace std;


#define DEF_BLOCK_SIZE (1u << 20)

size_t const block_size = DEF_BLOCK_SIZE;
char buff[block_size];

int
main(int argc, char * argv[])
{
  ios_base::openmode mode = ios_base::out;
  vector<ostream*> out;

  char c;
  while ((c = getopt(argc, argv, "a")) != -1) {
    switch (c) {
    case 'a':
      mode |= ios_base::ate;
      break;
    default:
      cerr << "unrecognized option: " << c << endl;
      exit(EXIT_FAILURE);
    }
  }

  // ignore SIGPIPE
  struct sigaction act;
  act.sa_handler=SIG_IGN;
  sigemptyset(&act.sa_mask);
  act.sa_flags=0;
  sigaction(SIGPIPE, &act, NULL);

  out.push_back(&cout);
  for (int i = optind; i < argc; ++i) {
    ofstream * tmp = new ofstream(argv[i], mode);
    if (tmp->fail()) {
      cerr << "error opening file [" << argv[i] << "]: " << strerror(errno) << endl;
      exit(EXIT_FAILURE);
    }
    out.push_back(tmp);
  }

  bool done;
  do {
    cin.read(buff, block_size);
    streamsize n = cin.gcount();
    done = true;
    for (size_t i = 0; i < out.size(); ++i) {
      if (!out[i]->fail()) {
	done = false;
	out[i]->write(buff, n);
      }
    }
  } while (not cin.eof() and not done);

  for (size_t i = 1; i < out.size(); ++i) {
    ((ofstream *)out[i])->close();
    delete out[i];
  }

  return EXIT_SUCCESS;
}
