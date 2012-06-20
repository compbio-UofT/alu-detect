#ifndef Interval_hpp_
#define Interval_hpp_

using namespace std;

#include <ostream>
#include <algorithm>


template<class T>
class Interval
{
public:
  T endpoint[2];

  T& operator [](int i) { return endpoint[i]; }
  T operator[](int i) const { return endpoint[i]; }
  void sort() { if (endpoint[0] > endpoint[1]) swap(endpoint[0], endpoint[1]); }

  bool intersects(const Interval<T>& other) {
    return (endpoint[0] <= other.endpoint[0] and other.endpoint[0] <= endpoint[1])
      or (other.endpoint[0] <= endpoint[0] and endpoint[0] <= other.endpoint[1]);
  }
};

template <class T>
ostream&
operator <<(ostream& ostr, const Interval<T>& i)
{
  ostr << "[" << i[0] << "," << i[1] << "]";
  return ostr;
}


#endif
