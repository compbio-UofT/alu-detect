#ifndef ScoreMatrix_hpp_
#define ScoreMatrix_hpp_

using namespace std;

#include <vector>


class ScoreMatrix
{
public:
  vector<int> score_;

  ScoreMatrix(int, int, int = 0);

  int& score(const char x, const char y) {
    return score_[256 * int(x) + int(y)];
  }
};


#endif
