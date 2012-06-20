#ifndef _UTIL_H
#define _UTIL_H

#include <assert.h>
#include <ctype.h>
#include <dirent.h>
#include <unistd.h>
#include <zlib.h>
#include <sys/types.h>
#include <sys/stat.h>

#include "stats.h"

#define MAX(_a, _b) ((_a) > (_b) ? (_a) : (_b))
#define MIN(_a, _b) ((_a) < (_b) ? (_a) : (_b))

struct _strbuf_t {
	char   *string;
	u_int	string_length;
	u_int	string_alloced;
};
typedef struct _strbuf_t * strbuf_t;

typedef long long int llint;

uint64_t	gettimeinusecs(void);
uint64_t	rdtsc(void);
double		cpuhz(void);
u_int		strchrcnt(const char *, const char);
bool		is_number(const char *);
bool		is_whitespace(const char *);
void		xstat(const char *, struct stat *);
void *		xmalloc(size_t);
void *		xmalloc_m(size_t, char const *);
void *		xmalloc_c(size_t, count_t *);
void *		xcalloc(size_t);
void *		xcalloc_m(size_t, char const *);
void *		xcalloc_c(size_t, count_t *);
void	       *xrealloc(void *, size_t);
void	       *xrealloc_c(void *, size_t, size_t, count_t *);
char	       *xstrdup(const char *);
double		ls_factorial(u_int);
double		ls_choose(int64_t, int64_t);
char	       *trim_brackets(char *);
void		progress_bar(FILE *, uint64_t, uint64_t, uint);
uint64_t	file_iterator(char *, void (*)(char *, struct stat *, void *), void *);
uint64_t	file_iterator_n(char **, int, void (*)(char *, struct stat *, void *), void *);
char const *	get_compiler(void);
char *		strrev(char *);
char	*	strtrim(char *);
strbuf_t	strbuf_create(void);
char *		strbuf_string(strbuf_t, int *);
void		strbuf_append(strbuf_t, char const *, ...);
void		strbuf_destroy(strbuf_t);
char *		comma_integer(uint64_t);
void		xgzwrite(gzFile, voidp, unsigned);
void		xgzread(gzFile, voidp, size_t);
size_t		removedups(void *, size_t, size_t, int (*)(void const *, void const *));
void		crash(int, int, char const *, ...);
void		logit(int, char const *, ...);
long long	nchoosek(int, int);
double		log_nchoosek(int, int);
void		cat(FILE *, FILE *);


static inline uint
ceil_div(uint a, uint b) {
  assert(b > 0);

  if (a == 0) return 0;
  return ((a - 1) / b) + 1;
}

/* compute base^power */
static inline size_t
power(size_t base, size_t exp)
{
  size_t result = 1;

  while (exp > 0) {
    if ((exp % 2) == 1)
      result *= base;
    base *= base;
    exp /= 2;
  }

  return result;
}

/* compute 4^power */
static inline llint
power4(int exp) {
  return (llint)1 << (2 * exp);
}


/* standard quality value from probability of error */
static inline int
qv_from_pr_err(double pr_err)
{
  if (pr_err > .99999999)
    return 0;
  else if (pr_err < 1E-25)
    return 250;
  else
    return (int)(-10.0 * log(pr_err) / log(10.0));
}

static inline int
qv_from_pr_corr(double pr_corr)
{
  return qv_from_pr_err(1 - pr_corr);
}

static inline double
pr_err_from_qv(int qv)
{
  if (qv <= 0)
    return .99999999;
  else if (qv >= 250)
    return 1E-25;
  else
    return pow(10.0, -(double)qv/10.0);
}


static inline int
double_to_neglog(double x, int shift = 1000)
{
  return (int)((double)shift * -log(x));
}


static inline double
neglog_to_double(int y, int shift = 1000)
{
  return exp(-(double)y / (double)shift);
}


static inline double
normal_cdf(double x, double mean, double stddev)
{
  double y = (x - mean) / stddev;
  if (y < 0) y = -y;
  double b0 = 0.2316419;
  double b1 = 0.319381530;
  double b2 = -0.356563782;
  double b3 = 1.781477937;
  double b4 = -1.821255978;
  double b5 = 1.330274429;
  double pi = 3.141592653589;
  double t = 1.0 / (1.0 + b0 * y);
  double res = (exp(- y * y / 2) / sqrt(2.0 * pi)) * ((((b5 * t + b4) * t + b3) * t + b2) * t + b1) * t;
  if (x > mean) res = 1 - res;
  return res;
}


#endif
