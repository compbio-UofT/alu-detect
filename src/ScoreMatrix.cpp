#include "ScoreMatrix.hpp"

#include <cctype>


ScoreMatrix::ScoreMatrix(int match_score, int mismatch_score, int )
  : score_(256 * 256, mismatch_score)
{
  // match
  score('a', 'a') = match_score;
  score('c', 'c') = match_score;
  score('g', 'g') = match_score;
  score('t', 't') = match_score;

  // half-match, if any

  // symmetrize with respect to case
  for (int i = 0; i < 256; i++) {
    if (!isalpha(i))
      continue;
    for (int j = 0; j < 256; j++) {
      if (!isalpha(j))
	continue;
      score(i, j) = score(tolower(i), tolower(j));
    }
  }
}
