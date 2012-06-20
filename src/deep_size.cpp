#include "deep_size.hpp"


size_t size_below(bool) { return 0; }
size_t size_below(char) { return 0; }
size_t size_below(int) { return 0; }
size_t size_below(long long int) { return 0; }
size_t size_below(void*) { return 0; }


size_t
size_below(const string& s)
{
  return s.capacity() * sizeof(char);
}
