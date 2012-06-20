#ifndef deep_size_hpp
#define deep_size_hpp

using namespace std;

#include <string>
#include <vector>


size_t size_below(bool);
size_t size_below(char);
size_t size_below(int);
size_t size_below(long long int);
size_t size_below(void*);
size_t size_below(const string&);


template<class T, class U>
size_t
size_below(const pair<T,U>& p)
{
  return size_below(p.first) + size_below(p.second);
}


template<class T>
size_t
size_below(const vector<T>& v)
{
  size_t res = v.capacity() * sizeof(T);
  for (size_t i = 0; i < v.size(); ++i) {
    res += size_below(v[i]);
  }
  return res;
}


#endif
