#include "inalign_core.hpp"

#include <cstdlib>
#include <cassert>
#include <climits>
#include <cstring>
#include <iostream>
#include <iomanip>

#include "globals.hpp"
#include "Clone.hpp"
#include "DNASequence.hpp"
#include "Fasta.hpp"
#include "common.hpp"
#include "ScoreMatrix.hpp"


int const cell_nw = 0;
int const cell_n = 1;
int const cell_w = 2;
int const cell_max_nw = 3;
int const cell_ptr_w = 4;
int const cell_se = 5;
int const cell_s = 6;
int const cell_e = 7;
int const cell_max_se = 8;
int const cell_ptr_e = 9;
int const cell_break = 10;
int const cell_ptr_overlap_w = 11;
int const cell_ptr_overlap_e = 12;
int const cell_max_nw_nobreak = 13;
int const cell_ptr_w_nobreak = 14;

int match_score_ref	= 10;
int match_score_rep	= 9;
int mismatch_score	= -15;
int db_gap_open_score	= -33;
int db_gap_extend_score	= -7;
int qr_gap_open_score	= -33;
int qr_gap_extend_score	= -3;
int db_breakpoint_score	= -100;
int db_breakpoint_gap_score	= -2;

#define MIN_READ_LEN 20

#define GOOD_REFERENCE_ALIGNMENT_THRESHOLD 0.7
#define GOOD_REPEAT_ALIGNMENT_THRESHOLD 0.7
#define MIN_MATCHES_REFERENCE 10
#define MIN_MATCHES_REPEAT 10
#define MAX_PCT_FLANKING_REGION 0.2
#define BREAKPOINT_SAFETY_LEN 8
#define REPEAT_CALL_SAFETY_DIFF 50
#define NULL_HYPOTHESIS_SAFETY_THRESHOLD 500
#define REPEAT_TAIL_LENGTH 40
#define MAX_LEN_FLANKING_DUPLICATION 100

//double const bp_consensus_threshold = .6; // min fraction of clones supporting a bp (should be >.5)
int const min_reads_to_bypass_evidence_against_bp = 3;

ScoreMatrix ref_scoreMatrix(match_score_ref, mismatch_score);
ScoreMatrix rep_scoreMatrix(match_score_rep, mismatch_score);


class read_entry
{
 public:
  const string*	name;
  const char*	seq;

  vector<vector<string> >	db_align;
  vector<vector<string> >	qr_align;
  vector<vector<int> >	db_match_start;
  vector<vector<int> >	db_match_end;
  vector<vector<int> >	qr_match_start;
  vector<vector<int> >	qr_match_end;
  vector<vector<int> >	score;
  vector<vector<int> >	mqv;
  vector<int>		start_in_base;
  vector<int>		num_pieces;

  int		mapped_to_base;
  int		mapped_to_repeat;
  int		evidence_against_bp_start;
  int		evidence_against_bp_end;
  int           read_len;

  read_entry()
    : name(NULL)
    , seq(NULL)
    , mapped_to_base(0)
    , mapped_to_repeat(0) {}
};

class pair_entry
{
 public:
  read_entry	re[2];
  const string*	clone_name;

  vector<int>	repeat_support;

  int		left_read;

  pair_entry()
    : clone_name(NULL) {}
};

typedef struct {
  int score[15];
  int back[15];
} cell;


inline void
max4_with_index(int score_a, int back_a,
		int score_b, int back_b,
		int score_c, int back_c,
		int score_d, int back_d,
		int * score_max, int * back_max) {
  int _score_max;
  int _back_max;
  if (score_max == NULL) score_max = &_score_max;
  if (back_max == NULL) back_max = &_back_max;

  *score_max = score_a;
  *back_max = back_a;
  if (score_b > *score_max) {
    *score_max = score_b;
    *back_max = back_b;
  }
  if (score_c > *score_max) {
    *score_max = score_c;
    *back_max = back_c;
  }
  if (score_d > *score_max) {
    *score_max = score_d;
    *back_max = back_d;
  }
}

inline void
max3_with_index(int score_a, int back_a,
		int score_b, int back_b,
		int score_c, int back_c,
		int * score_max, int * back_max) {
  int _score_max;
  int _back_max;
  if (score_max == NULL) score_max = &_score_max;
  if (back_max == NULL) back_max = &_back_max;

  *score_max = score_a;
  *back_max = back_a;
  if (score_b > *score_max) {
    *score_max = score_b;
    *back_max = back_b;
  }
  if (score_c > *score_max) {
    *score_max = score_c;
    *back_max = back_c;
  }
}

inline void
max2_with_index(int score_a, int back_a,
		int score_b, int back_b,
		int * score_max, int * back_max) {
  int _score_max;
  int _back_max;
  if (score_max == NULL) score_max = &_score_max;
  if (back_max == NULL) back_max = &_back_max;

  *score_max = score_a;
  *back_max = back_a;
  if (score_b > *score_max) {
    *score_max = score_b;
    *back_max = back_b;
  }
}


int
get_bp_gap_score(const char * qr, int start, int end)
{
  int result = db_breakpoint_gap_score;
  for (int i = start + 1; i <= end; i++)
    if (qr[i] != qr[i - 1])
      result += db_breakpoint_gap_score;
  return result;
}


void
dp_from_nw(const char * db, const char * qr, int db_len, int qr_len,
	   vector<vector<cell> >& m,
	   int overlap_start = 0)
{
  int i, j;

  if (overlap_start == 0)
    overlap_start = db_len;

  for (i = 0; i <= db_len; i++) {
    m[0][i].score[cell_nw] = 0;
    m[0][i].score[cell_n] = INT_MIN/2;
    m[0][i].score[cell_w] = INT_MIN/2;
    m[0][i].score[cell_max_nw] = 0;
  }
  for (j = 1; j <= qr_len; j++) {
    m[j][0].score[cell_nw] = INT_MIN/2;
    m[j][0].score[cell_n] = INT_MIN/2;
    m[j][0].score[cell_w] = INT_MIN/2;
    m[j][0].score[cell_max_nw] = m[j][0].score[cell_nw];
    m[j][0].back[cell_ptr_w] = 0;
    m[j][0].back[cell_ptr_overlap_w] = 0;
  }

  for (j = 1; j <= qr_len; j++) {
    for (i = 1; i <= db_len; i++) {

      /* nw */
      max3_with_index(m[j-1][i-1].score[cell_nw], cell_nw,
		      m[j-1][i-1].score[cell_n], cell_n,
		      m[j-1][i-1].score[cell_w], cell_w,
		      &m[j][i].score[cell_nw], &m[j][i].back[cell_nw]);
      m[j][i].score[cell_nw] += ref_scoreMatrix.score(db[i-1], qr[j-1]);
      //(db[i-1] == qr[j-1]? match_score_ref : mismatch_score);

      /* n */
      max3_with_index(m[j-1][i].score[cell_nw] + qr_gap_open_score, cell_nw,
		      m[j-1][i].score[cell_n], cell_n,
		      m[j-1][i].score[cell_w] + qr_gap_open_score, cell_w,
		      &m[j][i].score[cell_n], &m[j][i].back[cell_n]);
      m[j][i].score[cell_n] += qr_gap_extend_score;

      /* w */
      max3_with_index(m[j][i-1].score[cell_nw] + db_gap_open_score, cell_nw,
		      m[j][i-1].score[cell_n] + db_gap_open_score, cell_n,
		      m[j][i-1].score[cell_w], cell_w,
		      &m[j][i].score[cell_w], &m[j][i].back[cell_w]);
      m[j][i].score[cell_w] += db_gap_extend_score;

      /* max_nw */
      max3_with_index(m[j][i].score[cell_nw], cell_nw,
		      m[j][i].score[cell_n], cell_n,
		      m[j][i].score[cell_w], cell_w,
		      &m[j][i].score[cell_max_nw], &m[j][i].back[cell_max_nw]);

      /* ptr_w */
      max2_with_index(m[j][i].score[cell_max_nw], i,
		      m[j][m[j][i-1].back[cell_ptr_w]].score[cell_max_nw], m[j][i-1].back[cell_ptr_w],
		      NULL, &m[j][i].back[cell_ptr_w]);

      /* ptr_overlap_w */
      if (i <= overlap_start)
	m[j][i].back[cell_ptr_overlap_w] = i;
      else
	max2_with_index(m[j][i].score[cell_max_nw], i,
			m[j][m[j][i-1].back[cell_ptr_overlap_w]].score[cell_max_nw], m[j][i-1].back[cell_ptr_overlap_w],
			NULL, &m[j][i].back[cell_ptr_overlap_w]);

    }
  }
}


void
dp_from_se(const char * db, const char * qr, int db_len, int qr_len,
	   vector<vector<cell> >& m,
	   int overlap_end = 0)
{
  int i, j;

  if (overlap_end == 0)
    overlap_end = 1;

  for (i = db_len + 1; i >= 1; i--) {
    m[qr_len + 1][i].score[cell_se] = 0;
    m[qr_len + 1][i].score[cell_s] = INT_MIN/2;
    m[qr_len + 1][i].score[cell_e] = INT_MIN/2;
    m[qr_len + 1][i].score[cell_max_se] = 0;
  }
  for (j = qr_len; j >= 1; j--) {
    m[j][db_len + 1].score[cell_se] = INT_MIN/2;
    m[j][db_len + 1].score[cell_s] = INT_MIN/2;
    m[j][db_len + 1].score[cell_e] = INT_MIN/2;
    m[j][db_len + 1].score[cell_max_se] = INT_MIN/2;
    m[j][db_len + 1].back[cell_ptr_e] = db_len + 1;
    m[j][db_len + 1].back[cell_ptr_overlap_e] = db_len + 1;
  }

  for (j = qr_len; j >= 1; j--) {
    for (i = db_len; i >= 1; i--) {

      /* se */
      max3_with_index(m[j+1][i+1].score[cell_se], cell_se,
		      m[j+1][i+1].score[cell_s], cell_s,
		      m[j+1][i+1].score[cell_e], cell_e,
		      &m[j][i].score[cell_se], &m[j][i].back[cell_se]);
      m[j][i].score[cell_se] += ref_scoreMatrix.score(db[i-1], qr[j-1]);
      //(db[i-1] == qr[j-1]? match_score_ref : mismatch_score);

      /* s */
      max3_with_index(m[j+1][i].score[cell_se] + qr_gap_open_score, cell_se,
		      m[j+1][i].score[cell_s], cell_s,
		      m[j+1][i].score[cell_e] + qr_gap_open_score, cell_e,
		      &m[j][i].score[cell_s], &m[j][i].back[cell_s]);
      m[j][i].score[cell_s] += qr_gap_extend_score;

      /* e */
      max3_with_index(m[j][i+1].score[cell_se] + db_gap_open_score, cell_se,
		      m[j][i+1].score[cell_s] + db_gap_open_score, cell_s,
		      m[j][i+1].score[cell_e], cell_e,
		      &m[j][i].score[cell_e], &m[j][i].back[cell_e]);
      m[j][i].score[cell_e] += db_gap_extend_score;

      /* max_se */
      max3_with_index(m[j][i].score[cell_se], cell_se,
		      m[j][i].score[cell_s], cell_s,
		      m[j][i].score[cell_e], cell_e,
		      &m[j][i].score[cell_max_se], &m[j][i].back[cell_max_se]);

      /* ptr_e */
      max2_with_index(m[j][i].score[cell_max_se], i,
		      m[j][m[j][i+1].back[cell_ptr_e]].score[cell_max_se], m[j][i+1].back[cell_ptr_e],
		      NULL, &m[j][i].back[cell_ptr_e]);

      /* ptr_overlap_e */
      if (i >= overlap_end)
	m[j][i].back[cell_ptr_overlap_e] = i;
      else
	max2_with_index(m[j][i].score[cell_max_se], i,
			m[j][m[j][i+1].back[cell_ptr_overlap_e]].score[cell_max_se], m[j][i+1].back[cell_ptr_overlap_e],
			NULL, &m[j][i].back[cell_ptr_overlap_e]);

    }
  }
}


void
dp_from_nw_with_base(const char * db, const char * qr, int db_len, int qr_len,
		     vector<vector<cell> >& m,
		     vector<vector<cell> >& m_base,
		     int overlap_start, int overlap_end)
{
  int i, j;

  for (i = 0; i <= db_len; i++) {
    m[0][i].score[cell_nw] = 0;
    m[0][i].score[cell_n] = INT_MIN/2;
    m[0][i].score[cell_w] = INT_MIN/2;
    m[0][i].score[cell_max_nw] = 0;
    m[0][i].score[cell_break] = INT_MIN/2;
  }

  for (j = 1; j <= qr_len; j++) {
    // choose point to split between reference and repeat mapping
    // j+1 first letter in qr mapped to repeat
    m[j][0].score[cell_break] = m_base[j][m_base[j][overlap_end - 1].back[cell_ptr_overlap_w]].score[cell_max_nw];
    m[j][0].back[cell_break] = j;
    for (int j_ref = j - 1; j_ref > 1; j_ref--) {
      max2_with_index(m[j][0].score[cell_break], m[j][0].back[cell_break],
		      m_base[j_ref][m_base[j_ref][overlap_end - 1].back[cell_ptr_overlap_w]].score[cell_max_nw] + get_bp_gap_score(qr, j_ref + 1, j), j_ref,
		      &m[j][0].score[cell_break], &m[j][0].back[cell_break]);
    }
    m[j][0].score[cell_break] += db_breakpoint_score;

    m[j][0].score[cell_nw] = INT_MIN/2;
    m[j][0].score[cell_n] = INT_MIN/2;
    m[j][0].score[cell_w] = INT_MIN/2;
    m[j][0].score[cell_max_nw] = m[j][0].score[cell_break];
    m[j][0].back[cell_max_nw] = cell_break;
    m[j][0].back[cell_ptr_w] = 0;
    m[j][0].score[cell_max_nw_nobreak] = m[j][0].score[cell_nw];
    m[j][0].back[cell_max_nw_nobreak] = cell_nw;
    m[j][0].back[cell_ptr_w_nobreak] = 0;
  }

  for (j = 1; j <= qr_len; j++) {
    for (i = 1; i <= db_len; i++) {

      /* nw */
      max4_with_index(m[j-1][i-1].score[cell_nw], cell_nw,
		      m[j-1][i-1].score[cell_n], cell_n,
		      m[j-1][i-1].score[cell_w], cell_w,
		      m[j-1][0].score[cell_break], cell_break,
		      &m[j][i].score[cell_nw], &m[j][i].back[cell_nw]);
      m[j][i].score[cell_nw] += rep_scoreMatrix.score(db[i-1], qr[j-1]);
      //(db[i-1] == qr[j-1]? match_score_rep : mismatch_score);

      /* n */
      max4_with_index(m[j-1][i].score[cell_nw] + qr_gap_open_score, cell_nw,
		      m[j-1][i].score[cell_n], cell_n,
		      m[j-1][i].score[cell_w] + qr_gap_open_score, cell_w,
		      m[j-1][0].score[cell_break] + qr_gap_open_score, cell_break, // this should never be used
		      &m[j][i].score[cell_n], &m[j][i].back[cell_n]);
      m[j][i].score[cell_n] += qr_gap_extend_score;

      /* w */
      max4_with_index(m[j][i-1].score[cell_nw] + db_gap_open_score, cell_nw,
		      m[j][i-1].score[cell_n] + db_gap_open_score, cell_n,
		      m[j][i-1].score[cell_w], cell_w,
		      m[j][0].score[cell_break] + db_gap_open_score, cell_break, // this should never be used
		      &m[j][i].score[cell_w], &m[j][i].back[cell_w]);
      m[j][i].score[cell_w] += db_gap_extend_score;

      /* max_nw */
      max4_with_index(m[j][i].score[cell_nw], cell_nw,
		      m[j][i].score[cell_n], cell_n,
		      m[j][i].score[cell_w], cell_w,
		      m[j][0].score[cell_break], cell_break,
		      &m[j][i].score[cell_max_nw], &m[j][i].back[cell_max_nw]);

      /* max_nw_nobreak */
      max3_with_index(m[j][i].score[cell_nw], cell_nw,
                      m[j][i].score[cell_n], cell_n,
                      m[j][i].score[cell_w], cell_w,
                      &m[j][i].score[cell_max_nw_nobreak], &m[j][i].back[cell_max_nw_nobreak]);

      /* ptr_w */
      max2_with_index(m[j][i].score[cell_max_nw], i,
		      m[j][m[j][i-1].back[cell_ptr_w]].score[cell_max_nw], m[j][i-1].back[cell_ptr_w],
		      NULL, &m[j][i].back[cell_ptr_w]);

      /* ptr_w_nobreak */
      max2_with_index(m[j][i].score[cell_max_nw_nobreak], i,
                      m[j][m[j][i-1].back[cell_ptr_w_nobreak]].score[cell_max_nw_nobreak], m[j][i-1].back[cell_ptr_w_nobreak],
                      NULL, &m[j][i].back[cell_ptr_w_nobreak]);

    }
  }

  // connect end of repeat mapping with reference mapping
  for (j = qr_len - 1; j >= 1; j--) {
    // compute the score of matching the j-th letter in qr to repeat,
    // followed by a breakpoint sequence, then by reference matches
    m[j][db_len + 1].score[cell_break] = m_base[j + 1][m_base[j + 1][overlap_start + 1].back[cell_ptr_overlap_e]].score[cell_max_se];
    m[j][db_len + 1].back[cell_break] = j + 1;
    for (int j_ref = j + 2; j_ref <= qr_len; j_ref++) {
      max2_with_index(m[j][db_len + 1].score[cell_break], m[j][db_len + 1].back[cell_break],
		      m_base[j_ref][m_base[j_ref][overlap_start + 1].back[cell_ptr_overlap_e]].score[cell_max_se] + get_bp_gap_score(qr, j + 1, j_ref - 1), j_ref,
		      &m[j][db_len + 1].score[cell_break], &m[j][db_len + 1].back[cell_break]);
    }
    m[j][db_len + 1].score[cell_break] += m[j][m[j][db_len].back[cell_ptr_w_nobreak]].score[cell_max_nw_nobreak] + db_breakpoint_score;
  }

}


void
backtrace_from_se(const char * db, const char * qr, int db_len, int qr_len,
		  vector<vector<cell> >& m,
		  int db_start, int qr_start,
		  int min_db_start,
		  int use_break,
		  string& db_align, string& qr_align,
		  int& db_match_start, int& db_match_end, int& qr_match_start, int& qr_match_end,
		  int& score, int& mqv)
{
  int i, j, k, l, l_end;
  int j_prev = -1, j_next, j_exit_loop;

  assert(0 < min_db_start && min_db_start <= db_start);

  //*dbalign_p = (char *)calloc(db_len + qr_len + 1 + qr_len, sizeof(char));
  //*qralign_p = (char *)calloc(db_len + qr_len + 1 + qr_len, sizeof(char));
  //db_align = string(db_len + qr_len + 1 + qr_len, char(' '));
  //qr_align = string(db_len + qr_len + 1 + qr_len, char(' '));
  vector<char> tmp_db_align(db_len + qr_len + 1 + qr_len);
  vector<char> tmp_qr_align(db_len + qr_len + 1 + qr_len);

  j = qr_start;
  /*
  i_max = db_start;
  for (i = i_max - 1; i >= min_db_start; i--)
    if ((!use_break && m[j][i].score[cell_max_nw] > m[j][i_max].score[cell_max_nw])
	|| (use_break && m[j][i].score[cell_max_nw_nobreak] > m[j][i_max].score[cell_max_nw_nobreak]))
      i_max = i;
  i = i_max;
  */
  /*
  i = m[j][db_start].back[!use_break? cell_ptr_w : cell_ptr_w_nobreak];
  assert(i > 0);
  //k = m[j][i].back[!use_break? cell_max_nw : cell_max_nw_nobreak];
  k = m[j][i].back[!use_break? cell_max_nw : cell_max_nw_nobreak];
  assert(k != cell_break);
  */
  if (!use_break && qr_start == qr_len) {
    // in reference, no jump into repeat after this
    i = m[j][db_start].back[cell_ptr_w];
    k = m[j][i].back[cell_max_nw];
    // check if the jump is unique:
    int max_cnt = 0;
    for (int ii = db_start; ii >= min_db_start; --ii) {
      if (m[j][ii].score[cell_max_nw] == m[j][i].score[cell_max_nw]) {
	++max_cnt;
      }
    }
    if (max_cnt < 1) {
      cerr << "warning: max_cnt<1: this is bad" << endl;
      exit(1);
    }
    if (max_cnt == 1) {
      mqv = 10;
    } else {
      mqv = 0;
    }
  } else if (!use_break && qr_start < qr_len) {
    // in reference, but followed by jump into repeat
    i = m[j][db_start].back[cell_ptr_overlap_w];
    k = m[j][i].back[cell_max_nw];
    // check if the jump is unique:
    int max_cnt = 0;
    for (int ii = db_start; ii >= min_db_start; --ii) {
      if (m[j][ii].score[cell_max_nw] == m[j][i].score[cell_max_nw]) {
	++max_cnt;
      }
    }
    if (max_cnt < 1) {
      cerr << "warning: max_cnt<1: this is bad" << endl;
      exit(1);
   }
    if (max_cnt == 1) {
      mqv = 10;
    } else {
      mqv = 0;
    }
  } else {
    // in repeat
    assert(use_break);
    i = m[j][db_start].back[cell_ptr_w_nobreak];
    k = m[j][i].back[cell_max_nw_nobreak];
    // check if the jump is unique:
    int max_cnt = 0;
    for (int ii = db_start; ii >= min_db_start; --ii) {
      if (m[j][ii].score[cell_max_nw_nobreak] == m[j][i].score[cell_max_nw_nobreak]) {
	++max_cnt;
      }
    }
    if (max_cnt < 1) {
      cerr << "warning: max_cnt<1: this is bad" << endl;
      exit(1);
    }
    if (max_cnt == 1) {
      mqv = 10;
    } else {
      mqv = 0;
    }
  }
  assert(i > 0);
  assert(k != cell_break);
  l = db_len + qr_len - 1;
  l_end = l + 1;
  score = m[j][i].score[cell_max_nw];
  db_match_end = i;
  qr_match_end = j;
  do {
    db_match_start = i;
    qr_match_start = j;
    switch (k) {
    case cell_nw:
      tmp_db_align[l] = db[i-1];
      tmp_qr_align[l] = qr[j-1];
      l--;
      k = m[j][i].back[k];
      j--;
      i--;
      break;
    case cell_n:
      tmp_db_align[l] = '-';
      tmp_qr_align[l] = qr[j-1];
      l--;
      k = m[j][i].back[k];
      j--;
      break;
    case cell_w:
      tmp_db_align[l] = db[i-1];
      tmp_qr_align[l] = '-';
      l--;
      k = m[j][i].back[k];
      i--;
      break;
    }
    if (k == cell_break) {
      i = 0;
    }
  } while (i > 0 && j > 0);
  j_exit_loop = j;

  if (use_break) {
    // extend left
    if (k == cell_break) {
      j_prev = m[j][0].back[cell_break];
      while (j > j_prev) {
	tmp_db_align[l] = '=';
	tmp_qr_align[l] = qr[j-1];
	l--;
	j--;
	qr_match_start--;
      }
    }

    // extend right
    if (qr_start < qr_len) {
      j = qr_start;
      j_next = m[j][db_len + 1].back[cell_break];
      while (j_next > j + 1) {
	tmp_db_align[l_end] = '=';
	tmp_qr_align[l_end] = qr[j];
	l_end++;
	j++;
	qr_match_end++;
      }
    }

    // compute score of this part only
    if (qr_start < qr_len) {
      assert(db_start == db_len);
      score += db_breakpoint_score + get_bp_gap_score(qr, qr_start + 1, j_next - 1);
	//(j_next - qr_start - 1) * db_breakpoint_gap_score;
    }
    if (k == cell_break) {
      score += ( - m[j_exit_loop][0].score[cell_break] ) + db_breakpoint_score + get_bp_gap_score(qr, j_prev + 1, j_exit_loop);
	//(j_exit_loop - j_prev) * db_breakpoint_gap_score;
    }
  }

  //memmove(*dbalign_p, &(*dbalign_p)[l + 1], (strlen(&(*dbalign_p)[l + 1]) + 1) * sizeof(char));
  //memmove(*qralign_p, &(*qralign_p)[l + 1], (strlen(&(*qralign_p)[l + 1]) + 1) * sizeof(char));
  //*dbalign_p = (char *)realloc(*dbalign_p, (strlen(*dbalign_p) + 1) * sizeof(char));
  //*qralign_p = (char *)realloc(*qralign_p, (strlen(*qralign_p) + 1) * sizeof(char));
  //db_align.substr(l + 1, l_end - (l + 1));
  //qr_align = qr_align.substr(l + 1, l_end - (l + 1));
  db_align = string(tmp_db_align.begin() + l + 1, tmp_db_align.begin() + l_end);
  qr_align = string(tmp_qr_align.begin() + l + 1, tmp_qr_align.begin() + l_end);
}


void
backtrace_from_nw(const char * db, const char * qr, int db_len, int qr_len,
		  vector<vector<cell> >& m,
		  int db_start, int qr_start, int max_db_start,
		  string& db_align, string& qr_align,
		  int& db_match_start, int& db_match_end, int& qr_match_start, int& qr_match_end,
		  int& score, int& mqv)
{
  int i, j, k, l;

  //*dbalign_p = (char *)calloc(db_len + qr_len + 1, sizeof(char));
  //*qralign_p = (char *)calloc(db_len + qr_len + 1, sizeof(char));
  //db_align = string(db_len + qr_len + 1, char(' '));
  //qr_align = string(db_len + qr_len + 1, char(' '));
  vector<char> tmp_db_align(db_len + qr_len + 1);
  vector<char> tmp_qr_align(db_len + qr_len + 1);

  j = qr_start;
  i = m[j][db_start].back[cell_ptr_overlap_e];
  /*
  i_max = db_start;
  for (i = i_max + 1; i <= max_db_start; i++)
    if (m[j][i].score[cell_max_se] > m[j][i_max].score[cell_se])
      i_max = i;
  i = i_max;
  */
  assert(i > 0);
  k = m[j][i].back[cell_max_se];
  assert(k != cell_break);

  // check if jump is unique
  int max_cnt = 0;
  for (int ii = db_start; ii <= max_db_start; ++ii) {
    if (m[j][ii].score[cell_max_se] == m[j][i].score[cell_max_se]) {
      ++max_cnt;
    }
  }
  if (max_cnt < 1) {
    cerr << "warning: max_cnt<1: this is bad" << endl;
    exit(1);
  }
  if (max_cnt == 1) {
    mqv = 10;
  } else {
    mqv = 0;
  }

  l = 0;
  score = m[j][i].score[cell_max_se];
  db_match_start = i;
  qr_match_start = j;
  do {
    db_match_end = i;
    qr_match_end = j;
    switch (k) {
    case cell_se:
      tmp_db_align[l] = db[i-1];
      tmp_qr_align[l] = qr[j-1];
      l++;
      k = m[j][i].back[k];
      j++;
      i++;
      break;
    case cell_s:
      tmp_db_align[l] = '-';
      tmp_qr_align[l] = qr[j-1];
      l++;
      k = m[j][i].back[k];
      j++;
      break;
    case cell_e:
      tmp_db_align[l] = db[i-1];
      tmp_qr_align[l] = '-';
      l++;
      k = m[j][i].back[k];
      i++;
      break;
    }
    if (k == cell_break) {
      i = db_len;
    }
  } while (i <= db_len && j <= qr_len);

  //*dbalign_p = (char *)realloc(*dbalign_p, (l + 1) * sizeof(char));
  //*qralign_p = (char *)realloc(*qralign_p, (l + 1) * sizeof(char));
  //db_align = db_align.substr(0, l + 1);
  //qr_align = qr_align.substr(0, l + 1);
  db_align = string(tmp_db_align.begin(), tmp_db_align.begin() + l);
  qr_align = string(tmp_qr_align.begin(), tmp_qr_align.begin() + l);
}


void
backtrace_with_base(const char * db_base, const char * db, const char * qr,
		    int db_base_len, int db_len, int qr_len,
		    vector<vector<cell> >& m_base,
		    vector<vector<cell> >& m,
		    int overlap_start, int overlap_end,
		    vector<string>& db_align, vector<string>& qr_align,
		    vector<int>& db_match_start, vector<int>& db_match_end,
		    vector<int>& qr_match_start, vector<int>& qr_match_end,
		    int& start_in_base, vector<int>& score, vector<int>& mqv)
{
  int j;
  int max_score, tmp_score, j_max_score = -1, align_type;

  // first, try alignment to db_base only
  max_score = m_base[qr_len][m_base[qr_len][db_base_len].back[cell_ptr_w]].score[cell_max_nw];
  align_type = 0;

  // next, try alignment that ends in db (might or might not cross from db_base into db)
  tmp_score = m[qr_len][m[qr_len][db_len].back[cell_ptr_w]].score[cell_max_nw];
  if (tmp_score > max_score) {
    max_score = tmp_score;
    align_type = 1;
    //assert(m[qr_len][db_len].back[cell_ptr_w] > 0);
  }

  // finally, try alignments that cross from db into db_base (might or might not cross from db_base to db)
  for (j = qr_len - 1; j >= 1; j--) {
    tmp_score = m[j][db_len + 1].score[cell_break];
    if (tmp_score > max_score) {
      max_score = tmp_score;
      align_type = 2;
      j_max_score = j;
    }
  }

  if (align_type == 0)
    {
      // align to db_base only
      start_in_base = 1;
      backtrace_from_se(db_base, qr, db_base_len, qr_len,
			m_base,
			db_base_len, qr_len, 1,
			0,
			db_align[0], qr_align[0],
			db_match_start[0], db_match_end[0], qr_match_start[0], qr_match_end[0],
			score[0], mqv[0]);
      assert(qr_match_start[0] == 1 && qr_match_end[0] == qr_len);
    }
  else if (align_type == 1)
    {
      // alignment ends in db
      backtrace_from_se(db, qr, db_len, qr_len,
			m,
			db_len, qr_len, 1,
			1,
			db_align[0], qr_align[0],
			db_match_start[0], db_match_end[0], qr_match_start[0], qr_match_end[0],
			score[0], mqv[0]);
      assert(qr_match_end[0] == qr_len);

      if (qr_match_start[0] == 1) {
	// alignment to db only
	start_in_base = 0;
      } else {
	// alignment starts in db_base
	start_in_base = 1;
	db_align[1] = db_align[0];
	qr_align[1] = qr_align[0];
	db_match_start[1] = db_match_start[0];
	db_match_end[1] = db_match_end[0];
	qr_match_start[1] = qr_match_start[0];
	qr_match_end[1] = qr_match_end[0];
	score[1] = score[0] - db_breakpoint_score;
	mqv[1] = mqv[0];

	backtrace_from_se(db_base, qr, db_base_len, qr_len,
			  m_base,
			  overlap_end - 1, qr_match_start[1] - 1, overlap_start,
			  0,
			  db_align[0], qr_align[0],
			  db_match_start[0], db_match_end[0], qr_match_start[0], qr_match_end[0],
			  score[0], mqv[0]);
	assert(qr_match_start[0] == 1 && qr_match_end[0] == qr_match_start[1] - 1);
	assert(db_match_end[0] >= overlap_start);
	assert(db_match_end[0] < overlap_end);
      }
    }
  else // align_type == 2
    {
      // alignment crosses from db into db_base
      backtrace_from_se(db, qr, db_len, qr_len,
                        m,
                        db_len, j_max_score, 1,
			1,
                        db_align[0], qr_align[0],
                        db_match_start[0], db_match_end[0], qr_match_start[0], qr_match_end[0],
                        score[0], mqv[0]);
      score[0] -= db_breakpoint_score;
      assert(qr_match_end[0] < qr_len);

      if (qr_match_start[0] == 1) {
	// alignment starts in db
	start_in_base = 0;
      } else {
        // algnment starts in db_base
	start_in_base = 1;
        db_align[1] = db_align[0];
        qr_align[1] = qr_align[0];
        db_match_start[1] = db_match_start[0];
        db_match_end[1] = db_match_end[0];
        qr_match_start[1] = qr_match_start[0];
        qr_match_end[1] = qr_match_end[0];
	score[1] = score[0] - db_breakpoint_score;
	mqv[1] = mqv[0];

        backtrace_from_se(db_base, qr, db_base_len, qr_len,
                          m_base,
                          overlap_end - 1, qr_match_start[1] - 1, overlap_start,
			  0,
                          db_align[0], qr_align[0],
                          db_match_start[0], db_match_end[0], qr_match_start[0], qr_match_end[0],
                          score[0], mqv[0]);
        assert(qr_match_start[0] == 1 && qr_match_end[0] == qr_match_start[1] - 1);
	assert(db_match_end[0] >= overlap_start);
	assert(db_match_end[0] < overlap_end);
      }
      int idx = (start_in_base ? 2 : 1);
      backtrace_from_nw(db_base, qr, db_base_len, qr_len,
			m_base,
			overlap_start + 1, qr_match_end[idx - 1] + 1, overlap_end,
			db_align[idx], qr_align[idx],
			db_match_start[idx], db_match_end[idx], qr_match_start[idx], qr_match_end[idx],
			score[idx], mqv[idx]);
      assert(qr_match_start[idx] == qr_match_end[idx - 1] + 1 && qr_match_end[idx] == qr_len);
      assert(db_match_start[idx] <= overlap_end);
      assert(db_match_start[idx] > overlap_start);
    }
}


void
fprint_reference_mapping(ostream& ostr,
			 const string* read_name,
			 const string& reference_name, long long int base_offset,
			 const string& db_align, const string& qr_align,
			 int db_match_start, int db_match_end,
			 int qr_match_start, int qr_match_end,
			 int score)
{
  ostr << "[" << *read_name << "] vs [" << reference_name << "]:" << '\n';
  ostr << db_align
       << " [" << base_offset + db_match_start - 1
       << "-" << base_offset + db_match_end - 1
       << "] (" << score
       << "(" << (int)((double)score /
		       (double)((qr_match_end - qr_match_start + 1) * match_score_ref) * 100.0)
       << "%))" << '\n';
  ostr << qr_align << '\n';
}


void
fprint_split_mapping(ostream& ostr,
		     const string* read_name,
		     const string& reference_name, long long int base_offset,
		     const string& repeat_name, int repeat_rc,
		     const vector<string>& db_align, const vector<string>& qr_align,
		     const vector<int> db_match_start, const vector<int> db_match_end,
		     const vector<int> qr_match_start, const vector<int> qr_match_end,
		     const vector<int> score,
		     int num_pieces, int start_in_base)
{
  ostr << "[" << *read_name << "] vs [" << reference_name
       << "+" << repeat_name << (repeat_rc? "-" : "+") << "]" << '\n';
  if (!start_in_base) 
    {
      if (num_pieces == 1) {
	ostr << "|" << db_align[0] << "|"
	     << " [|" << db_match_start[0] << "-" << db_match_end[0] << "|]"
	     << " (|" << score[0]
	     << "(" << (int)((double)score[0] /
			     (double)((qr_match_end[0] - qr_match_start[0] + 1) * match_score_rep) * 100.0)
	     << "%)|)" << '\n';
	ostr << "|" << qr_align[0] << "|" << '\n';
      } else {
	ostr << "|" << db_align[0] << "|" << db_align[1]
	     << " [|" << db_match_start[0] << "-" << db_match_end[0] << "|"
	     << base_offset + db_match_start[1] - 1 << "-" << base_offset + db_match_end[1] - 1
	     << "] (|" << score[0]
	     << "(" << (int)((double)score[0] /
			     (double)((qr_match_end[0] - qr_match_start[0] + 1) * match_score_rep) * 100.0)
	     << "%)|" << score[1]
	     << "(" << (int)((double)score[1]/(double)((qr_match_end[1] - qr_match_start[1] + 1) * match_score_ref) * 100.0)
	     << "%))" << '\n';
	ostr << "|" << qr_align[0] << "|" << qr_align[1] << '\n';
      }
    } 
  else // starts in base
    {
      if (num_pieces == 1) {
	ostr << db_align[0] << "|| [" << base_offset + db_match_start[0] - 1
	     << "-" << base_offset + db_match_end[0] - 1
	     << "||] (" << score[0]
	     << "(" << (int)((double)score[0] /
			     (double)((qr_match_end[0] - qr_match_start[0] + 1) * match_score_ref) * 100.0)
	     << "%)||)" << '\n';
	ostr << qr_align[0] << "||" << '\n';
      } else if (num_pieces == 2) {
	ostr << db_align[0] << "|" << db_align[1] << "|"
	     << " [" << base_offset + db_match_start[0] - 1
	     << "-" << base_offset + db_match_end[0] - 1
	     << "|" << db_match_start[1] << "-" << db_match_end[1] << "|]"
	     << " (" << score[0]
	     << "(" << (int)((double)score[0] /
			     (double)((qr_match_end[0] - qr_match_start[0] + 1) * match_score_ref) * 100.0)
	     << "%)|" << score[1]
	     << "(" << (int)((double)score[1] /
			     (double)((qr_match_end[1] - qr_match_start[1] + 1) * match_score_rep) * 100.0)
	     << "%)|)" << '\n';
	ostr << qr_align[0] << "|" << qr_align[1] << "|" << '\n';
      } else {
	ostr << db_align[0] << "|" << db_align[1] << "|" << db_align[2]
	     << " [" << base_offset + db_match_start[0] - 1
	     << "-" << base_offset + db_match_end[0] - 1
	     << "|" << db_match_start[1] << "-" << db_match_end[1]
	     << "|" << base_offset + db_match_start[2] - 1
	     << "-" << base_offset + db_match_end[2] - 1 << "]"
	     << " (" << score[0] << "|" << score[1] << "|" << score[2] << ")" << '\n';
	ostr << qr_align[0] << "|" << qr_align[1] << "|" << qr_align[2] << '\n';
      }
    }
}


void
run(const Contig& refContig,
    long long int refStart, long long int refEnd,
    long long int overlapStart, long long int overlapEnd,
    const vector<pair<Contig*,int> >& repeat,
    vector<Clone*>& clone,
    int only_update_mappings,
    ostream& out_str, ostream& err_str)
{
  const char* reference_seq = &refContig.seq[0].c_str()[refStart - refContig.seqOffset[0]];
  long long int base_offset = refStart + 1;
  int reference_len = int(refEnd - refStart + 1);
  int overlap_start = int(overlapStart - refStart + 1);
  int overlap_end = int(overlapEnd - refStart + 1);
  int n_repeats = repeat.size();
  int n_pairs = clone.size();
  int pair_id, nip, rep_id;

  err_str << "region=" << refContig.name << ":" << refStart + 1 << "-" << refEnd + 1
	  << ":" << overlapStart + 1 << "-" << overlapEnd + 1
	  << " repeats=" << n_repeats << " clones=" << n_pairs << '\n';

  err_str << "using reference_len [" << reference_len << "]"  << '\n';
  err_str << "using overlap_start [" << overlap_start << "]"  << '\n';
  err_str << "using overlap_end [" << overlap_end << "]"  << '\n';

  if (overlap_start >= overlap_end) {
    cerr << "warning: bad call to run(): " << *clone[0] << endl;
    return;
  }

  int max_repeat_len = -1;
  for (rep_id = 0; rep_id < n_repeats; ++rep_id) {
    max_repeat_len = max(max_repeat_len, (int)repeat[rep_id].first->len);
  }

  int max_read_len = -1;
  for (pair_id = 0; pair_id < n_pairs; pair_id++) {
    max_read_len = max(max_read_len, clone[pair_id]->read[0].len);
    max_read_len = max(max_read_len, clone[pair_id]->read[1].len);
  }

  // allocate global matrices
  cell zero_cell;
  memset((void*)&zero_cell, 0, sizeof(cell));
  vector<vector<cell> > m_reference(max_read_len + 2,
				    vector<cell>(reference_len + 2, zero_cell));
  vector<vector<cell> > m_repeat(max_read_len + 2,
				 vector<cell>(max_repeat_len + 2, zero_cell));

  // allocate structures
  vector<pair_entry> pe(n_pairs);
  for (pair_id = 0; pair_id < n_pairs; pair_id++) {
    pe[pair_id].repeat_support = vector<int>(1 + n_repeats);
    for (nip = 0; nip < 2; nip++) {
      pe[pair_id].re[nip].db_align = vector<vector<string> >(1 + n_repeats, vector<string>(3));
      pe[pair_id].re[nip].qr_align = vector<vector<string> >(1 + n_repeats, vector<string>(3));
      pe[pair_id].re[nip].db_match_start = vector<vector<int> >(1 + n_repeats, vector<int>(3));
      pe[pair_id].re[nip].db_match_end = vector<vector<int> >(1 + n_repeats, vector<int>(3));
      pe[pair_id].re[nip].qr_match_start = vector<vector<int> >(1 + n_repeats, vector<int>(3));
      pe[pair_id].re[nip].qr_match_end = vector<vector<int> >(1 + n_repeats, vector<int>(3));
      pe[pair_id].re[nip].score = vector<vector<int> >(1 + n_repeats, vector<int>(3));
      pe[pair_id].re[nip].mqv = vector<vector<int> >(1 + n_repeats, vector<int>(3));
      pe[pair_id].re[nip].start_in_base = vector<int>(1 + n_repeats);
      pe[pair_id].re[nip].num_pieces = vector<int>(1 + n_repeats);
    }
  }

  // copy reads into old structures
  for (pair_id = 0; pair_id < n_pairs; pair_id++) {
    Clone& c = *clone[pair_id];
    err_str << c << '\n';
    pe[pair_id].clone_name = &c.name;
    pe[pair_id].left_read = (c.read[0].st == 0? 0 : 1); // XXX: ASSUMES "opp-in"!!!
    for (nip = 0; nip < 2; nip++) {
      read_entry* rep = &pe[pair_id].re[nip];
      rep->name = &c.read[nip].name;
      rep->seq = c.read[nip].seq.c_str();
      rep->read_len = c.read[nip].len;
    }
  }

  vector<int> evidence_against_bp(overlap_end - overlap_start, int(0));
  vector<int> total_repeat_score(n_repeats, int(0));

  // align each read to the reference, then to all repeats
  for (pair_id = 0; pair_id < n_pairs; pair_id++) {
    for (nip = 0; nip < 2; nip++) {
      read_entry * rep = &pe[pair_id].re[nip];
      if (rep->read_len < MIN_READ_LEN) {
	err_str << "skipping read [" << *(rep->name) << "]: too short" << '\n';
	continue;
      }

      dp_from_nw(reference_seq, rep->seq,
		 reference_len, rep->read_len,
		 m_reference, overlap_start);
      dp_from_se(reference_seq, rep->seq,
		 reference_len, rep->read_len,
		 m_reference, overlap_end);
      backtrace_from_se(reference_seq, rep->seq,
			reference_len, rep->read_len,
			m_reference,
			reference_len, rep->read_len, 1, 0,
			rep->db_align[0][0], rep->qr_align[0][0],
			rep->db_match_start[0][0], rep->db_match_end[0][0],
			rep->qr_match_start[0][0], rep->qr_match_end[0][0],
			rep->score[0][0], rep->mqv[0][0]);

      // if read maps well to the reference only
      if (rep->score[0][0] >= (int)(GOOD_REFERENCE_ALIGNMENT_THRESHOLD
				    * (double)(rep->read_len * match_score_ref)))
	{
	  rep->mapped_to_base = 1;

	  fprint_reference_mapping(err_str,
				   rep->name,
				   refContig.name, base_offset,
				   rep->db_align[0][0], rep->qr_align[0][0],
				   rep->db_match_start[0][0], rep->db_match_end[0][0],
				   rep->qr_match_start[0][0], rep->qr_match_end[0][0],
				   rep->score[0][0]);

	  // mark breakpoints which are ruled out by this read
	  rep->evidence_against_bp_start = max(overlap_start,
					       rep->db_match_start[0][0] + BREAKPOINT_SAFETY_LEN - 1);
	  rep->evidence_against_bp_end = min(overlap_end - 1,
					     rep->db_match_end[0][0] - BREAKPOINT_SAFETY_LEN);

	  for (int i = rep->evidence_against_bp_start; i <= rep->evidence_against_bp_end; i++) {
	    evidence_against_bp[i - overlap_start]++;
	  }
	}
      else // map it to repeats
	{
	  for (rep_id = 0; rep_id < n_repeats; rep_id++) {
	    vector<string>& db_align = rep->db_align[1 + rep_id];
	    vector<string>& qr_align = rep->qr_align[1 + rep_id];
	    vector<int>& db_match_start = rep->db_match_start[1 + rep_id];
	    vector<int>& db_match_end = rep->db_match_end[1 + rep_id];
	    vector<int>& qr_match_start = rep->qr_match_start[1 + rep_id];
	    vector<int>& qr_match_end = rep->qr_match_end[1 + rep_id];
	    vector<int>& score = rep->score[1 + rep_id];
	    vector<int>& mqv = rep->mqv[1 + rep_id];
	    int& start_in_base = rep->start_in_base[1 + rep_id];

	    dp_from_nw_with_base(repeat[rep_id].first->seq[repeat[rep_id].second].c_str(), rep->seq,
				 repeat[rep_id].first->len, rep->read_len,
				 m_repeat, m_reference,
				 overlap_start, overlap_end);
	    backtrace_with_base(reference_seq, repeat[rep_id].first->seq[repeat[rep_id].second].c_str(), rep->seq,
				reference_len, repeat[rep_id].first->len, rep->read_len,
				m_reference, m_repeat,
				overlap_start, overlap_end,
				db_align, qr_align,
				db_match_start, db_match_end,
				qr_match_start, qr_match_end,
				start_in_base, score, mqv);
	    // count number of alignment pieces
	    rep->num_pieces[1 + rep_id] = 1;
	    if (qr_match_end[0] < rep->read_len) {
	      rep->num_pieces[1 + rep_id]++;
	      if (qr_match_end[1] < rep->read_len)
		rep->num_pieces[1 + rep_id]++;
	    }

	    fprint_split_mapping(err_str,
				 rep->name,
				 refContig.name, base_offset,
				 repeat[rep_id].first->name, repeat[rep_id].second,
				 db_align, qr_align,
				 db_match_start, db_match_end,
				 qr_match_start, qr_match_end,
				 score,
				 rep->num_pieces[1 + rep_id], start_in_base);

	    // accounting stage: check mapping quality
	    switch (rep->num_pieces[1 + rep_id]) {

	    case 1:
	      if (!start_in_base
		  && score[0] >= (int)(GOOD_REPEAT_ALIGNMENT_THRESHOLD
				       * (double)(rep->read_len * match_score_rep))) {
		// no breakpoint evidence, but count score
		total_repeat_score[rep_id] += score[0];
		rep->mapped_to_repeat = 1;
	      } else {
		rep->num_pieces[1 + rep_id] = 0;
	      }
	      break;

	    case 2:
	      int qr_len_in_reference;
	      int qr_len_in_repeat;
	      int score_reference;
	      int score_repeat;
	      int db_len_in_repeat;

	      if (start_in_base) {
		qr_len_in_reference = qr_match_end[0] - qr_match_start[0] + 1;
		qr_len_in_repeat = qr_match_end[1] - qr_match_start[1] + 1;
		score_reference = score[0];
		score_repeat = score[1];
		db_len_in_repeat = db_match_end[1] - db_match_start[1] + 1;
	      } else {
		qr_len_in_reference = qr_match_end[1] - qr_match_start[1] + 1;
		qr_len_in_repeat = qr_match_end[0] - qr_match_start[0] + 1;
		score_reference = score[1];
		score_repeat = score[0];
		db_len_in_repeat = db_match_end[0] - db_match_start[0] + 1;
	      }
	      if (qr_len_in_reference >= MIN_MATCHES_REFERENCE
		  && qr_len_in_repeat >= MIN_MATCHES_REPEAT
		  && (qr_len_in_repeat - db_len_in_repeat) <= (int)(MAX_PCT_FLANKING_REGION * (double)db_len_in_repeat)
		  && score_reference >= (int)(GOOD_REFERENCE_ALIGNMENT_THRESHOLD * (double)(qr_len_in_reference * match_score_ref))
		  && score_repeat >= (int)(GOOD_REPEAT_ALIGNMENT_THRESHOLD * (double)(qr_len_in_repeat * match_score_rep))) {
		total_repeat_score[rep_id] += score[start_in_base? 1 : 0];
		rep->mapped_to_repeat = 1;
	      } else {
		rep->num_pieces[1 + rep_id] = 0;
	      }
	      break;

	    default:
	      assert(rep->num_pieces[1 + rep_id] == 3);
	      if (false // don't allow these mappings
		  && qr_match_end[0] - qr_match_start[0] + 1 >= MIN_MATCHES_REFERENCE
		  && qr_match_end[1] - qr_match_start[1] + 1 >= MIN_MATCHES_REPEAT
		  && qr_match_end[2] - qr_match_start[2] + 1 >= MIN_MATCHES_REFERENCE
		  && ((qr_match_end[1] - qr_match_start[1]) - (db_match_end[1] - db_match_start[1])) <= (int)(MAX_PCT_FLANKING_REGION * (double)(db_match_end[1] - db_match_start[1]))
		  && score[0] >= (int)(GOOD_REFERENCE_ALIGNMENT_THRESHOLD * (double)((qr_match_end[0] - qr_match_start[0] + 1) * match_score_ref))
		  && score[1] >= (int)GOOD_REPEAT_ALIGNMENT_THRESHOLD * (double)((qr_match_end[1] - qr_match_start[1] + 1) * match_score_rep)
		  && score[2] >= (int)(GOOD_REFERENCE_ALIGNMENT_THRESHOLD * (double)((qr_match_end[2] - qr_match_start[2] + 1) * match_score_ref))) {
		total_repeat_score[rep_id] += score[1];
		rep->mapped_to_repeat = 1;
	      } else {
		rep->num_pieces[1 + rep_id] = 0;
	      }
	      break;
	    }

	  }
	  if (!rep->mapped_to_repeat) {
	    err_str << "read [" << *(rep->name) << "] not mapped to this region" << '\n';
	  }
	}
    }
  }

  //int total_clones = n_pairs;
  vector<int> total_clones_supporting_insertion(1 + n_repeats, int(0));
  vector<int> total_reads_spanning_insert_start(1 + n_repeats, int(0));
  vector<int> total_reads_spanning_insert_end(1 + n_repeats, int(0));
  vector<int> min_repeat_position(1 + n_repeats, int(0));
  vector<int> max_repeat_position(1 + n_repeats, int(0));
  vector<int> bp_mapped_to_repeat(1 + n_repeats, int(0));
  vector<int> sum_scores_mapped_to_repeat(1 + n_repeats, int(0));
  vector<int> repeat_score(1 + n_repeats, int(0));
  vector<int> bp_null_hypothesis(1 + n_repeats, int(0));
  vector<int> sum_scores_null_hypothesis(1 + n_repeats, int(0));
  vector<int> null_hypothesis_score(1 + n_repeats, int(0));
  vector<int> null_hypothesis_possible(1 + n_repeats, int(0));
  vector<vector<int> > evidence_for_insert_start(1 + n_repeats,
						 vector<int>(overlap_end - overlap_start, int(0)));
  vector<vector<int> > evidence_for_insert_end(1 + n_repeats,
					       vector<int>(overlap_end - overlap_start, int(0)));

  // done mapping reads, time to aggregate results
  // compute which clones support which repeats
  for (rep_id = 0; rep_id < n_repeats; rep_id++) {
    for (pair_id = 0; pair_id < n_pairs; pair_id++) {
      // does this clone support this repeat?
      read_entry * left_re = &pe[pair_id].re[pe[pair_id].left_read];
      read_entry * right_re = &pe[pair_id].re[(pe[pair_id].left_read + 1) % 2];
      if (// at least one side touches the reference
	  (left_re->mapped_to_base || left_re->num_pieces[1 + rep_id] == 2
	   || right_re->mapped_to_base || right_re->num_pieces[1 + rep_id] == 2)
	  && // at least one side touches the repeat
	  (left_re->num_pieces[1 + rep_id] > 0
	   || right_re->num_pieces[1 + rep_id] > 0)
	  && // both sides are either mapped or too short
	  ((left_re->mapped_to_base || left_re->num_pieces[1 + rep_id] > 0 || left_re->read_len < MIN_READ_LEN)
	   && (right_re->mapped_to_base || right_re->num_pieces[1 + rep_id] > 0 || right_re->read_len < MIN_READ_LEN)))
	{
	  pe[pair_id].repeat_support[1 + rep_id] = 1;
	  total_clones_supporting_insertion[1 + rep_id]++;

	  for (int k = 0; k < 2; k++) {
	    nip = (pe[pair_id].left_read + k) % 2;
	    read_entry * rep = &pe[pair_id].re[nip];
	    vector<int>& db_match_start = rep->db_match_start[1 + rep_id];
	    vector<int>& db_match_end = rep->db_match_end[1 + rep_id];
	    vector<int>& qr_match_start = rep->qr_match_start[1 + rep_id];
	    vector<int>& qr_match_end = rep->qr_match_end[1 + rep_id];
	    vector<int>& score = rep->score[1 + rep_id];
	    vector<int>& mqv = rep->mqv[1 + rep_id];
	    int& start_in_base = rep->start_in_base[1 + rep_id];
	    int num_pieces = rep->num_pieces[1 + rep_id];

	    if (rep->mapped_to_base || rep->read_len < MIN_READ_LEN)
	      continue;

	    // update min, max, and base counts
	    int repeat_chunk = 0;
	    if (num_pieces == 2 && start_in_base)
	      repeat_chunk = 1;

	    if (min_repeat_position[1 + rep_id] == 0
		|| db_match_start[repeat_chunk] < min_repeat_position[1 + rep_id]) {
	      min_repeat_position[1 + rep_id] = db_match_start[repeat_chunk];
	    }
	    if (max_repeat_position[1 + rep_id] == 0
		|| db_match_end[repeat_chunk] > max_repeat_position[1 + rep_id]) {
	      max_repeat_position[1 + rep_id] = db_match_end[repeat_chunk];
	    }
	    bp_mapped_to_repeat[1 + rep_id] +=
	      qr_match_end[repeat_chunk] - qr_match_start[repeat_chunk] + 1;

	    sum_scores_mapped_to_repeat[1 + rep_id] += score[repeat_chunk];
	    bp_null_hypothesis[1 + rep_id] += rep->read_len;
	    sum_scores_null_hypothesis[1 + rep_id] += rep->score[0][0];

	    if (num_pieces == 2) {
	      // evidence of one breakpoint
	      if (start_in_base and mqv[0] >= 4) {
		// spanning insert start
		total_reads_spanning_insert_start[1 + rep_id]++;
		evidence_for_insert_start[1 + rep_id][db_match_end[0] - overlap_start]++;
	      } else if (!start_in_base and mqv[1] >= 4) {
		// spanning insert end
		total_reads_spanning_insert_end[1 + rep_id]++;
		evidence_for_insert_end[1 + rep_id][db_match_start[1] - 1 - overlap_start]++;
	      }
	    }
	  }
	}
    } // end for loop over pair_id
    if (sum_scores_mapped_to_repeat[1 + rep_id] < 0)
      sum_scores_mapped_to_repeat[1 + rep_id] = 0;
    repeat_score[1 + rep_id] = (int)(1000.0 * (double)sum_scores_mapped_to_repeat[1 + rep_id]
				     / (double)(bp_mapped_to_repeat[1 + rep_id] * match_score_rep));
    if (sum_scores_null_hypothesis[1 + rep_id] < 0)
      sum_scores_null_hypothesis[1 + rep_id] = 0;
    null_hypothesis_score[1 + rep_id] = (int)(1000.0 * (double)sum_scores_null_hypothesis[1 + rep_id]
					      / (double)(bp_null_hypothesis[1 + rep_id] * match_score_ref));
    if (null_hypothesis_score[1 + rep_id] > NULL_HYPOTHESIS_SAFETY_THRESHOLD)
      null_hypothesis_possible[1 + rep_id] = 1;
  }

  // compute the winner
  int max_rep_id = -1;
  for (rep_id = 0; rep_id < n_repeats; rep_id++)
    if ((max_rep_id < 0 && total_repeat_score[rep_id] > 0)
	|| (max_rep_id >= 0 && total_repeat_score[rep_id] > total_repeat_score[max_rep_id]))
      max_rep_id = rep_id;

  // also the runners up
  //int n_possible_repeats = 0;
  //int * possible_repeat = (int *)malloc(n_repeats * sizeof(int));
  vector<int> possible_repeat;

  if (max_rep_id >= 0) {
    //n_possible_repeats = 1;
    possible_repeat.push_back(max_rep_id);
    for (rep_id = 0; rep_id < n_repeats; rep_id++)
      if (rep_id != max_rep_id
	  && total_clones_supporting_insertion[1 + rep_id] == total_clones_supporting_insertion[1 + max_rep_id]
	  && sum_scores_mapped_to_repeat[1 + rep_id] > sum_scores_mapped_to_repeat[1 + max_rep_id] - REPEAT_CALL_SAFETY_DIFF)
	possible_repeat.push_back(rep_id);
  }

  if (max_rep_id < 0) {
    err_str << "no repeats supported by the reads" << '\n';
  } else {
    err_str << "repeats supported by the reads (with score):" << '\n';
    for (size_t i = 0; i < possible_repeat.size(); i++) {
      rep_id = possible_repeat[i];
      err_str << repeat[rep_id].first->name << "\t" << total_repeat_score[rep_id] << '\n';
    }
  }

  if (possible_repeat.size() > 0)
    rep_id = possible_repeat[0];
  else
    rep_id = -1;
  for (pair_id = 0; pair_id < n_pairs; pair_id++) {
    for (int k = 0; k < 2; k++) {
      nip = (pe[pair_id].left_read + k) % 2;
      read_entry * rep = &pe[pair_id].re[nip];
      vector<int>& db_match_start = rep->db_match_start[1 + rep_id];
      vector<int>& db_match_end = rep->db_match_end[1 + rep_id];
      vector<int>& qr_match_start = rep->qr_match_start[1 + rep_id];
      vector<int>& qr_match_end = rep->qr_match_end[1 + rep_id];
      vector<int>& mqv = rep->mqv[1 + rep_id];
      int& start_in_base = rep->start_in_base[1 + rep_id];
      int num_pieces = rep->num_pieces[1 + rep_id];

      // update mappings:
      Read& r = clone[pair_id]->read[nip];
      //int repeat_st = (clone[pair_id]->read[nip].st + repeat[rep_id].second) % 2;
      int repeat_st = (rep_id >= 0? repeat[rep_id].second : 0);
      Mapping m;
      m.qr = &r;
      r.mapping.clear();
      r.mappedToRepeatSt[0] = false;
      r.mappedToRepeatSt[1] = false;

      if (rep->read_len < MIN_READ_LEN
	  || (!rep->mapped_to_base && (rep_id < 0 || num_pieces == 0 || num_pieces > 2)))
	{
	  err_str << *(rep->name) << "\t*" << '\n';
	}
      else if (rep->mapped_to_base)
	{
	  err_str << *(rep->name) << '\t'
		  << refContig.name << '\t'
		  << base_offset - 1 + rep->db_match_start[0][0] << '\t'
		  << base_offset - 1 + rep->db_match_end[0][0] << '\n';
	  m.db = (Contig*)&refContig;
	  m.dbPos[0] = base_offset - 1 + rep->db_match_start[0][0] - 1;
	  m.dbPos[1] = base_offset - 1 + rep->db_match_end[0][0] - 1;
	  m.qrPos[0] = rep->qr_match_start[0][0] - 1;
	  m.qrPos[1] = rep->qr_match_end[0][0] - 1;
	  m.st = 0;
	  m.mqv = rep->mqv[0][0];
	  m.is_ref = true;
	  r.mapping.push_back(m);
	}
      else if (num_pieces == 1)
	{
	  err_str << *(rep->name) << '\t'
		  << repeat[rep_id].first->name << (repeat[rep_id].second ? "-" : "+") << '\t'
		  << db_match_start[0] << '\t'
		  << db_match_end[0] << '\n';
	  m.db = repeat[rep_id].first;
	  m.dbPos[0] = db_match_start[0] - 1;
	  m.dbPos[1] = db_match_end[0] - 1;
	  m.qrPos[0] = qr_match_start[0] - 1;
	  m.qrPos[1] = qr_match_end[0] - 1;
	  m.st = repeat[rep_id].second;
	  m.mqv = mqv[0];
	  m.is_ref = false;
	  r.mapping.push_back(m);

	  r.mappedToRepeatSt[repeat_st] = true;
	}
      else if (num_pieces == 2 && start_in_base)
	{
	  err_str << *(rep->name) << '\t'
		  << refContig.name << '\t'
		  << base_offset - 1 + db_match_start[0] << '\t'
		  << base_offset - 1 + db_match_end[0] << '\t'
		  << repeat[rep_id].first->name << (repeat[rep_id].second ? "-" : "+") << '\t'
		  << db_match_start[1] << '\t'
		  << db_match_end[1] << '\n';
	  m.db = (Contig*)&refContig;
	  m.dbPos[0] = base_offset - 1 + db_match_start[0] - 1;
	  m.dbPos[1] = base_offset - 1 + db_match_end[0] - 1;
	  m.qrPos[0] = qr_match_start[0] - 1;
	  m.qrPos[1] = qr_match_end[0] - 1;
	  m.st = 0;
	  m.mqv = mqv[0];
	  m.is_ref = true;
	  r.mapping.push_back(m);

	  m.db = repeat[rep_id].first;
	  m.dbPos[0] = db_match_start[1] - 1;
	  m.dbPos[1] = db_match_end[1] - 1;
	  m.qrPos[0] = qr_match_start[1] - 1;
	  m.qrPos[1] = qr_match_end[1] - 1;
	  m.st = repeat[rep_id].second;
	  m.mqv = mqv[1];
	  m.is_ref = false;
	  r.mapping.push_back(m);

	  r.mappedToRepeatSt[repeat_st] = true;
	}
      else if (num_pieces == 2 && !start_in_base)
	{
	  err_str << *(rep->name) << '\t'
		  << repeat[rep_id].first->name << (repeat[rep_id].second ? "-" : "+") << '\t'
		  << db_match_start[0] << '\t'
		  << db_match_end[0] << '\t'
		  << refContig.name << '\t'
		  << base_offset - 1 + db_match_start[1] << '\t'
		  << base_offset - 1 + db_match_end[1] << '\n';
	  m.db = repeat[rep_id].first;
	  m.dbPos[0] = db_match_start[0] - 1;
	  m.dbPos[1] = db_match_end[0] - 1;
	  m.qrPos[0] = qr_match_start[0] - 1;
	  m.qrPos[1] = qr_match_end[0] - 1;
	  m.st = repeat[rep_id].second;
	  m.mqv = mqv[0];
	  m.is_ref = false;
	  r.mapping.push_back(m);

	  m.db = (Contig*)&refContig;
	  m.dbPos[0] = base_offset - 1 + db_match_start[1] - 1;
	  m.dbPos[1] = base_offset - 1 + db_match_end[1] - 1;
	  m.qrPos[0] = qr_match_start[1] - 1;
	  m.qrPos[1] = qr_match_end[1] - 1;
	  m.st = 0;
	  m.mqv = mqv[1];
	  m.is_ref = true;
	  r.mapping.push_back(m);

	  r.mappedToRepeatSt[repeat_st] = true;
	}
    }
  }

  if (not only_update_mappings)
    // not updating mappings; this is a "real" run, print mapping summary
    {
      for (pair_id = 0; pair_id < n_pairs; pair_id++) {
	for (int k = 0; k < 2; k++) {
	  nip = (pe[pair_id].left_read + k) % 2;
	  read_entry * rep = &pe[pair_id].re[nip];
	  int& start_in_base = rep->start_in_base[1 + rep_id];
	  int num_pieces = rep->num_pieces[1 + rep_id];

	  if (k == 0)
	    err_str << *pe[pair_id].clone_name << ':';

	  if (rep->read_len < MIN_READ_LEN
	      || (!rep->mapped_to_base && (rep_id < 0 || num_pieces == 0 || num_pieces > 2)))
	    err_str << setw(12) << " -";
	  else if (rep->mapped_to_base)
	    err_str << setw(12) << " ref";
	  else if (num_pieces == 1)
	    err_str << setw(12) << " rep";
	  else if (num_pieces == 2 && start_in_base)
	    err_str << setw(12) << " ref+rep";
	  else if (num_pieces == 2 && !start_in_base)
	    err_str << setw(12) << " rep+ref";

	  if (k == 1) {
	    if (rep_id < 0 || !pe[pair_id].repeat_support[1 + rep_id])
	      err_str << " !";
	    err_str << '\n';
	  }
	}
      }
    }

  // if this was a "mock" run, we're done
  if (only_update_mappings)
    return;

  // do we have a candidate repeat
  if (rep_id < 0)
    return;

  // is it supported by at least 2 clones
  if (total_clones_supporting_insertion[1 + rep_id] <= 1)
    return;

  // if we mapped only to poly-a tail, disregard mapping
  if ((repeat[rep_id].second == 0
       and min_repeat_position[1 + rep_id] > repeat[rep_id].first->len - REPEAT_TAIL_LENGTH)
      or (repeat[rep_id].second == 1
	  and max_repeat_position[1 + rep_id] < REPEAT_TAIL_LENGTH)) {
    err_str << "only matched polyA tail; disregarding" << '\n';
    return;
  }

  // this is a "real" run, and there is a candidate repeat
  for (pair_id = 0; pair_id < n_pairs; pair_id++) {
    if (!pe[pair_id].repeat_support[1 + rep_id])
      continue;

    for(int k = 0; k < 2; k++) {
      nip = (pe[pair_id].left_read + k) % 2;
      read_entry * rep = &pe[pair_id].re[nip];
      vector<string>& db_align = rep->db_align[1 + rep_id];
      vector<string>& qr_align = rep->qr_align[1 + rep_id];
      vector<int>& db_match_start = rep->db_match_start[1 + rep_id];
      vector<int>& db_match_end = rep->db_match_end[1 + rep_id];
      vector<int>& qr_match_start = rep->qr_match_start[1 + rep_id];
      vector<int>& qr_match_end = rep->qr_match_end[1 + rep_id];
      vector<int>& score = rep->score[1 + rep_id];
      int& start_in_base = rep->start_in_base[1 + rep_id];
      //int num_pieces = rep->num_pieces[1 + rep_id];

      if (rep->read_len < MIN_READ_LEN)
	continue;
      if (rep->mapped_to_base)
	{
	  fprint_reference_mapping(err_str,
				   rep->name,
				   refContig.name, base_offset,
				   rep->db_align[0][0], rep->qr_align[0][0],
				   rep->db_match_start[0][0], rep->db_match_end[0][0],
				   rep->qr_match_start[0][0], rep->qr_match_end[0][0],
				   rep->score[0][0]);
	} 
      else
	{
	  fprint_split_mapping(err_str,
			       rep->name,
			       refContig.name, base_offset,
			       repeat[rep_id].first->name, repeat[rep_id].second,
			       db_align, qr_align,
			       db_match_start, db_match_end,
			       qr_match_start, qr_match_end,
			       score,
			       rep->num_pieces[1 + rep_id], start_in_base);
	}
    }
  }

  // process breakpoint evidence:
  // take the max insert start and min insert end
  //   with >= 2 support (unless evidence=1, in which case take lone evidence)
  int ins_start_evidence = 0;
  for (int i = 0; i < overlap_end - overlap_start; i++) {
    ins_start_evidence += evidence_for_insert_start[1 + rep_id][i];
  }
  int ins_start = -1;
  if (ins_start_evidence > 0) {
    for (int i = overlap_end - overlap_start - 1; i >= 0; --i) {
      //if (evidence_for_insert_start[1 + rep_id][i] > int(bp_consensus_threshold * double(ins_start_evidence))) {
      if ((ins_start_evidence == 1 && evidence_for_insert_start[1 + rep_id][i] >= 1)
	  || evidence_for_insert_start[1 + rep_id][i] >= 2) {
	ins_start = i;
	break;
      }
    }
    if (ins_start < 0) {
      ins_start = -2;
    }
  }

  int ins_end_evidence = 0;
  for (int i = 0; i < overlap_end - overlap_start; i++) {
    ins_end_evidence += evidence_for_insert_end[1 + rep_id][i];
  }
  int ins_end = -1;
  if (ins_end_evidence > 0) {
    for (int i = 0; i < overlap_end - overlap_start; i++) {
      //if (evidence_for_insert_end[1 + rep_id][i] > int(bp_consensus_threshold * double(ins_end_evidence))) {
      if ((ins_end_evidence == 1 && evidence_for_insert_end[1 + rep_id][i] >= 1)
	  || evidence_for_insert_end[1 + rep_id][i] >= 2) {
	ins_end = i;
	break;
      }
    }
    if (ins_end < 0) {
      ins_end = -2;
    }
  }

  if (ins_start == -2) { // ambiguous ins_start
    err_str << "warning: ambiguous evidence for insert start" << '\n';
  } else if (ins_start == -1) {
    err_str << "no evidence for insert start" << '\n';
  } else if (evidence_for_insert_start[1 + rep_id][ins_start] < min_reads_to_bypass_evidence_against_bp
	     && evidence_against_bp[ins_start]) {
    err_str << "insert start at [" << base_offset + overlap_start - 1 + ins_start
	    << "] supported by [" << evidence_for_insert_start[1 + rep_id][ins_start]
	    << "] reads but disproved by [" << evidence_against_bp[ins_start]
	    << "] reads" << '\n';
    ins_start = -2;
  }
  if (ins_start < 0) {
    total_reads_spanning_insert_start[1 + rep_id] = 0;
  }

  if (ins_end == -2) { // ambiguous ins_end
    err_str << "warning: ambiguous evidence for insert end" << '\n';
  } else if (ins_end == -1) {
    err_str << "no evidence for insert end" << '\n';
  } else if (evidence_for_insert_end[1 + rep_id][ins_end] < min_reads_to_bypass_evidence_against_bp
	     && evidence_against_bp[ins_end]) {
    err_str << "insert end at [" << base_offset + overlap_start - 1 + ins_end + 1
	    << "] supported by [" << evidence_for_insert_end[1 + rep_id][ins_end]
	    << "] reads but disproved by [" << evidence_against_bp[ins_end]
	    << "] reads" << '\n';
    ins_end = -2;
  }
  if (ins_end < 0) {
    total_reads_spanning_insert_end[1 + rep_id] = 0;
  }

  int tsd_len = 0;
  if (ins_start >= 0 && ins_end >= 0) {
    tsd_len = ins_start - ins_end;
    if (ins_start < ins_end) { // TS loss, not duplication; tsd_len < 0
      swap(ins_start, ins_end);
    }

    if (abs(tsd_len) > MAX_LEN_FLANKING_DUPLICATION) {
      err_str << "warning: tsd_len [" << tsd_len
	      << "], possible evidence of 2 distinct insertions" << '\n';
    }
  }

  // trim insert range when only one endpoint is known
  if (ins_start >= 0 && ins_end < 0) {
    ins_end = ins_start - 50;
  }
  if (ins_start < 0 && ins_end >= 0) {
    ins_start = ins_end + 50;
  }
  assert(ins_start < 0 or ins_end < 0 or ins_end <= ins_start);

  // make alternative score nonnegative
  if (sum_scores_null_hypothesis[1 + rep_id] <= 0)
    sum_scores_null_hypothesis[1 + rep_id] = 0;

  // notes to self:
  //   base_offset: 1-based, absolute
  //   overlap_start/end: 1-based, off base_offset
  //   ins_start/end: 0-based, off overlap_start
  out_str << refContig.name << '\t'
	  << (ins_end >= 0?
	      (base_offset - 1) + (overlap_start - 1) + ins_end
	      : (base_offset - 1) + overlap_start - 1) << '\t' // 0-based, closed
	  << (ins_start >= 0?
	      (base_offset - 1) + (overlap_start - 1) + ins_start + 1
	      : (base_offset - 1) + (overlap_end - 1)) + 1 << '\t' // 0-based, open
	  << repeat[rep_id].first->name << '\t'
	  << repeat_score[1 + rep_id] << '\t'
	  << (repeat[rep_id].second? "-" : "+") << '\t'
	  << (repeat[rep_id].second?
	      int(repeat[rep_id].first->len) + 1 - min_repeat_position[1 + rep_id]
	      : min_repeat_position[1 + rep_id]) << '\t'
	  << (repeat[rep_id].second?
	      int(repeat[rep_id].first->len) + 1 - max_repeat_position[1 + rep_id]
	      : max_repeat_position[1 + rep_id]) << '\t'
	  << total_clones_supporting_insertion[1 + rep_id] << '\t'
	  << total_reads_spanning_insert_start[1 + rep_id] << '\t'
	  << total_reads_spanning_insert_end[1 + rep_id] << '\t'
	  << tsd_len << '\t';
  if (possible_repeat.size() > 1 || null_hypothesis_possible[1 + rep_id]) {
    if (null_hypothesis_possible[1 + rep_id]) {
      out_str << "Null:" << null_hypothesis_score[1 + rep_id];
    }
    for (size_t i = 1; i < possible_repeat.size(); i++) {
      if (i > 1 || null_hypothesis_possible[1 + rep_id])
	out_str << ',';
      out_str << repeat[possible_repeat[i]].first->name << ':' << repeat_score[1 + rep_id];
    }
  } else {
    out_str << '.';
  }
  out_str << '\n';
}
