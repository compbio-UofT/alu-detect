#include <iostream>
#include <cstdlib>
#include <boost/thread.hpp>

#include "concurrent_queue.hpp"


using namespace std;

char const * prog_name;
int verbosity;

class Block
{
public:
  static const int size = (1 << 20);
  char * block;
  int load;
  int bid;

  Block() : block(NULL), load(0), bid(0) {}
};

concurrent_queue<Block> q;


void
out_flow()
{
  stringstream err_str;
  if (verbosity > 0)
    cerr << "out_flow: started\n";
  Block b;
  while (q.wait_and_pop(b)) {
    cout.write(b.block, b.load);

    if (verbosity > 0) {
      err_str << "out_flow: wrote block " << b.bid << "\n";
      cerr << err_str.str();
      err_str.str("");
    }
    delete [] b.block;
  }
  if (verbosity > 0)
    cerr << "out_flow: ended\n";
}

void
in_flow()
{
  stringstream err_str;
  if (verbosity > 0)
    cerr << "in_flow: started\n";
  int i = 0;
  while (!cin.eof()) {
    Block b;
    b.block = new char[Block::size];
    if (b.block == NULL) {
      cerr << "prog_name: out of memory\n";
      exit(1);
    }	

    cin.read(b.block, Block::size);
    b.load = cin.gcount();

    b.bid = i++;
    if (verbosity > 0) {
      err_str << "in_flow: read block " << b.bid << " (" << b.load << " bytes)\n";
      cerr << err_str.str();
      err_str.str("");
    }
    q.timed_wait_and_push(b);
  }
  q.set_done();
  if (verbosity > 0)
    cerr << "in_flow: ended\n";
}

void
usage(ostream& os)
{
  os << prog_name << ": copy stdin to stdout; keep reading stdin if stdout stuck" << endl;
}

int
main(int argc, char * argv[])
{
  prog_name = argv[0];
  ios_base::sync_with_stdio(false);

  char c;
  while ((c = getopt(argc, argv, "vh")) != -1) {
    switch (c) {
    case 'v':
      verbosity++;
      break;
    case 'h':
      usage(cout);
      exit(EXIT_SUCCESS);
      break;
    default:
      cerr << "unrecognized option: " << c << endl;
      usage(cerr);
      exit(EXIT_FAILURE);
    }
  }
  if (optind < argc) {
    cerr << "unrecognized parameter: " << argv[optind] << endl;
    usage(cerr);
    exit(EXIT_FAILURE);
  }

  boost::thread out_flow_thread(out_flow);
  in_flow();
  out_flow_thread.join();

  return EXIT_SUCCESS;
}
