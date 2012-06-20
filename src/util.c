#define __STDC_FORMAT_MACROS
#include <assert.h>
#include <inttypes.h>
#include <errno.h>
#include <math.h>
#include <stdarg.h>
#include <stdbool.h>
#include <stdio.h>
#include <stdlib.h>
#include <stdint.h>
#include <string.h>
#include <unistd.h>
#include <zlib.h>
#include <limits.h>
#include <sys/time.h>
#include <sys/types.h>

#include "util.h"


uint64_t
gettimeinusecs()
{
	struct timeval tv;
	gettimeofday(&tv, NULL);
	return ((uint64_t)tv.tv_sec * 1000000 + tv.tv_usec);	
}

uint64_t
rdtsc()
{
	uint32_t lo, hi;

#ifdef __GNUC__
	__asm__ __volatile__ ("rdtsc" : "=a" (lo), "=d" (hi));
#else
	asm("rdtsc" : "=a" (lo), "=d" (hi));
#endif

	return (((uint64_t)hi << 32) | lo);
}

double
cpuhz()
{
	uint64_t before;
	struct timeval tv1, tv2;
	int diff;

	/* XXX - abusive, poor man's calc; needs good (2ms) clock granularity */
	gettimeofday(&tv1, NULL);
	before = rdtsc();
	do {
		gettimeofday(&tv2, NULL);

		diff = tv2.tv_usec - tv1.tv_usec;
		if (diff < 0)
			diff = 1000000 - tv1.tv_usec + tv2.tv_usec;
	} while (diff < 2000);

	return (((double)(rdtsc() - before) / diff) * 1.0e6);
}

u_int
strchrcnt(const char *str, const char c)
{
	int i;

	i = 0;
	while (*str != '\0') {
		if (*str++ == c)
			i++;
	}

	return (i);
}

bool
is_number(const char *str)
{

	while (*str != '\0')
		if (!isdigit((int)*str++))
			return (false);

	return (true);
}

bool
is_whitespace(const char *str)
{

	while (*str != '\0')
		if (!isspace((int)*str++))
			return (false);

	return (true);
}

void
xstat(const char *path, struct stat *sbp)
{
	
	if (stat(path, sbp) != 0) {
		fprintf(stderr, "error: failed to stat [%s]: %s\n", path,
		    strerror(errno));
		exit(1);
	}
}

void *
xmalloc(size_t size)
{
	void *ptr;

	ptr = malloc(size);
	if (ptr == NULL) {
		fprintf(stderr, "error: malloc failed: %s\n", strerror(errno));
		exit(1);
	}

	return (ptr);
}

void * xmalloc_m(size_t size, char const * error_msg)
{
  void * p;

  p = malloc(size);
  if (p == NULL) {
    fprintf(stderr, "error: malloc failed: %s [%s]\n", strerror(errno), error_msg);
    exit(1);
  }
  return p;
}


void *
xmalloc_c(size_t size, count_t * c)
{
  if (c != NULL)
    count_add(c, size);
  return xmalloc(size);
}

void *
xcalloc(size_t size)
{
  void *ptr;

  ptr = calloc(size, 1);
  if (ptr == NULL) {
    fprintf(stderr, "error: calloc failed: %s\n", strerror(errno));
    exit(1);
  }

  return ptr;
}

void * xcalloc_m(size_t size, char const * error_msg)
{
  void * p;

  p = calloc(size, 1);
  if (p == NULL) {
    fprintf(stderr, "error: calloc failed: %s [%s]\n", strerror(errno), error_msg);
    exit(1);
  }
  return p;
}

void *
xcalloc_c(size_t size, count_t * c)
{
  if (c != NULL)
    count_add(c, size);
  return xcalloc(size);
}

void *
xrealloc(void *ptr, size_t size)
{
  ptr = realloc(ptr, size);
  if (ptr == NULL) {
    fprintf(stderr, "error: realloc failed: %s\n", strerror(errno));
    exit(1);
  }

  return (ptr);
}

void *
xrealloc_c(void *ptr, size_t size, size_t old_size, count_t * c)
{
  if (c != NULL)
    count_add(c, (int64_t)size - (int64_t)old_size);
  return xrealloc(ptr, size);
}


char *
xstrdup(const char *str)
{
	char *dup;

	assert(str != NULL);

	dup = strdup(str);
	if (dup == NULL) {
		fprintf(stderr, "error: strdup failed: %s\n", strerror(errno));
		exit(1);
	}

	return (dup);
}

/* factorial using stirling's approximation after 20 */
double
ls_factorial(u_int n)
{
	const double fact[21] = {
		1.0,
		1.0,
		2.0,
		6.0,
		24.0,
		120.0,
		720.0,
		5040.0,
		40320.0,
		362880.0,
		3628800.0,
		39916800.0,
		479001600.0,
		6227020800.0,
		87178291200.0,
		1307674368000.0,
		20922789888000.0,
		355687428096000.0,
		6402373705728000.0,
		121645100408832000.0,
		2432902008176640000.0
	};
	double a, b;

	if (n <= 20)
		return log (fact[n]);

	a = log(sqrt(2 * M_PI * n));
	b = n * log(n / M_E);

	return (a + b);
}

/* choose in log space */
double
ls_choose(int64_t n, int64_t k)
{
	double a, b, c;

	if (k < 0 || k > n)
		return (0);

	a = ls_factorial(n);
	b = ls_factorial(k);
	c = ls_factorial(n - k);

	return (a - (b + c));
}

char *
trim_brackets(char *str)
{

	if (str[0] == '[')
		str++;
	if (str[strlen(str) - 1] == ']')
		str[strlen(str) - 1] = '\0';

	return (str);
}


void
progress_bar(FILE *out, uint64_t at, uint64_t of, uint incr)
{
	static int lastperc, beenhere;
	static char whirly = '\\';

	char progbuf[52 + 16];
	int perc, i, j, dec;

	if (at == 0 && of == 0) {
		beenhere = lastperc = 0;
		whirly = '\\';
		return;
	}

	perc = (at * 100 * incr) / of;

	if (beenhere && perc == lastperc)
		return;

	beenhere = 1;
	lastperc = perc;

	dec = perc % incr;
	perc /= incr;

	/* any excuse to have a whirly gig */
	switch (whirly) {
	case '|':
		whirly = '/';
		break;
	case '/':
		whirly = '-';
		break;
	case '-':
		whirly = '\\';
		break;
	case '\\':
		whirly = '|';
		break;
	}
	if (at >= of)
		whirly = '|';

	progbuf[25] = whirly;
		
	for (i = j = 0; i <= 100; i += 2) {
		if (j != 25) {
			if (i <= perc)
				progbuf[j++] = '=';
			else
				progbuf[j++] = ' ';
		} else {
			j++;
		}
	}
	memset(&progbuf[51], 0, 16);	/* XXX - valgrind */

	if (incr == 100)
		fprintf(out, "\rProgress: [%s] %3d.%02d%%", progbuf, perc, dec);
	else if (incr == 10)
		fprintf(out, "\rProgress: [%s] %3d.%d%%", progbuf, perc, dec);
	else
		fprintf(out, "\rProgress: [%s] %3d%%", progbuf, perc);
	
	fflush(out);
}

/*
 * Given a path, if it's a regular file (or symlink to one), call fh on it.
 * If it's a directory, call fh on all regular files within it (or symlinks to
 * regular files).
 *
 * Returns the number of files fh was called on.
 */
uint64_t
file_iterator(char *path, void (*fh)(char *, struct stat *, void *),
    void *arg)
{
	char fpath[2048];
	struct stat sb;
	DIR *dp;
	struct dirent *de;
	uint64_t files;

	/* is a regular file... */
	xstat(path, &sb);
	if (S_ISREG(sb.st_mode)) {
		fh(path, &sb, arg);
		return (1);
	}

	/* is (hopefully) a directory... */
	dp = opendir(path);
	if (dp == NULL) {
		fprintf(stderr, "error: failed to open directory [%s]: %s\n",
		    path, strerror(errno));
		exit(1);
	}

	files = 0;
	while (1) {
		de = readdir(dp);
		if (de == NULL)
			break;

		strcpy(fpath, path);
		if (fpath[strlen(path) - 1] != '/')
			strcat(fpath, "/");
		strcat(fpath, de->d_name);
		xstat(fpath, &sb);

#if defined(DT_REG) && defined(DT_LNK)
		if (de->d_type != DT_REG && de->d_type != DT_LNK)
			continue;
#else
		if (!S_ISREG(sb.st_mode) && !S_ISLNK(sb.st_mode))
			continue;
#endif

		/* ensure it's a regular file or link to one */
		if (S_ISREG(sb.st_mode)) {
			fh(fpath, &sb, arg);
			files++;
		} else {
			fprintf(stderr, "warning: [%s] is neither a regular "
			    "file, nor a link to one; skipping...", fpath);
			continue;
		}
	}

	closedir(dp);

	return (files);
}

uint64_t
file_iterator_n(char **paths, int npaths,
    void (*fh)(char *, struct stat *, void *), void *arg)
{
	uint64_t files;
	int i;

	for (i = files = 0; i < npaths; i++)
		files += file_iterator(paths[i], fh, arg);

	return (files);
}

char const *
get_compiler()
{

#if defined(__GNUC__)
	if (strstr(__VERSION__, "Intel(R)"))
		return ("ICC " __VERSION__);
	else
		return ("GCC " __VERSION__);
#elif defined(__SUNPRO_C)
	return ("Sun Pro C");
#elif defined(__SUNPRO_CC)
	return ("Sun Pro C++");
#elif defined(__cplusplus)
	return ("unknown C++");
#else
	return ("unknown C");
#endif
}

/* reverse the string `str' in place */
char *
strrev(char *str)
{
	char c;
	int i, j;

	j = strlen(str) - 1;
	for (i = 0; i < j; i++, j--) {
		c = str[j];
		str[j] = str[i];
		str[i] = c;
	}

	return (str);
}

/* trim whitespace in `str' in place at beginning and end */
char *
strtrim(char *str)
{
	char *ret;

        assert(str != NULL);

        while (isspace((int)*str) && *str != '\0')
                str++;

        ret = str;

	if (*str != '\0') {
		while (*str != '\0')
			str++;

		str--;
		
		while (isspace((int)*str))
			str--;
		str++;
	}

        *str = '\0';

        return (ret);
}

strbuf_t
strbuf_create()
{
	strbuf_t sbp;

	sbp = (strbuf_t)xmalloc(sizeof(*sbp));
	memset(sbp, 0, sizeof(*sbp));
	sbp->string_alloced = 4096;
	sbp->string = (char *)xmalloc(4096);

	return (sbp);
}

char *
strbuf_string(strbuf_t sbp, int *length)
{

	assert(sbp->string_length == strlen(sbp->string));
	if (length != NULL)
		*length = sbp->string_length;
	return (xstrdup(sbp->string));
}

void
strbuf_append(strbuf_t sbp, char const *fmt, ...)
{
	va_list ap;
	int bytes;

	assert(sbp->string_length < sbp->string_alloced);
	if (sbp->string_alloced == sbp->string_length) {
		sbp->string_alloced += 4096;
		sbp->string = (char *)xrealloc(sbp->string, sbp->string_alloced);
	}

	va_start(ap, fmt);
	do {
		bytes = vsnprintf(&sbp->string[sbp->string_length],
		    sbp->string_alloced - sbp->string_length, fmt, ap);

		/* wasn't enough space. resize and try again */
		if (sbp->string_length + bytes >= sbp->string_alloced) {
			sbp->string_alloced += 4096;
			sbp->string = (char *)xrealloc(sbp->string, sbp->string_alloced);
		}
	} while (sbp->string_length + bytes >= sbp->string_alloced);
	va_end(ap);

	sbp->string_length += bytes;
	assert(sbp->string_length < sbp->string_alloced);
}

void
strbuf_destroy(strbuf_t sbp)
{

	free(sbp->string);
	free(sbp);
}

void xgzread(gzFile fp, voidp buf, size_t len) {
        size_t total = 0;
        while (total < len) {
                int res;
                int to_read = (len - total <= (size_t)INT_MAX? (int)(len - total) : INT_MAX);
                res = gzread(fp, (char *)buf+total, to_read);
                if (res <= 0) {
                        fprintf(stderr, "error: gzread returned %u\n", res);
                        exit(0);
                }
                total += (size_t)res;
        }

}

void xgzwrite(gzFile fp, voidp buf, unsigned len){
	uint total = 0;
	while (total < len){
		uint res;
		res = gzwrite(fp,(char *)buf+total,len-total);
		if (res <= 0){
			fprintf(stderr,"error: gzwrite returned %u\n",res);
			exit(0);
		}
		total += res;
	}
}

/*
 * Return a string on the stack corresponding to an unsigned integer that also
 * features commas. E.g.: int 1000 yields "1,000".
 *
 * There is a pool of sequentially allocated buffers returned, so this should
 * be safe to use multiple times in function arguments.
 */
char *
comma_integer(uint64_t val)
{
	static char rets[50][32];	// no malloc, allow uses in fn args, etc
	static int col = 0;

	char *ret = rets[(col++ % (sizeof(rets) / sizeof(rets[0])))];
	char str[sizeof(rets[0])];
	int skip, i, j;

	memset(str, 0, sizeof(str));	// XXX - shut up, valgrind
	snprintf(str, sizeof(str), "%" PRIu64, val);

	skip = 3 - (strlen(str) % 3);
	for (i = j = 0; str[i] != '\0'; i++) {
		if ((i + skip) % 3 == 0 && i != 0)
			ret[j++] = ',';
		ret[j++] = str[i];
	}
	ret[j] = '\0';

	return (ret);
}


// pre: array is (q)sorted
size_t removedups(void * a, size_t n, size_t sz, int (*cmp)(void const *, void const *))
{
  char * p = (char *)a;
  size_t m = 0;
  size_t i, j;
  i = 0;
  while (i < n) {
    j = i+1;
    while (j < n && !cmp(p+i*sz, p+j*sz)) j++;
    if (m < i)
      memcpy(p+m*sz, p+i*sz, sz);
    m++;
    i = j;
  }
  return m;
}


void
crash(int exit_code, int display_errno, char const * msg, ...)
{
  va_list fmtargs;
  char new_msg[strlen(msg) + 1000];

  if (display_errno)
    sprintf(new_msg, "error (%s): %s\n", strerror(errno), msg);
  else
    sprintf(new_msg, "error: %s\n", msg);

  va_start(fmtargs, msg);
  vfprintf(stderr, new_msg, fmtargs);
  va_end(fmtargs);

  exit(exit_code);
}


void
logit(int display_errno, char const * msg, ...)
{
  va_list fmtargs;
  char new_msg[strlen(msg) + 1000];

  if (display_errno)
    sprintf(new_msg, "note (%s): %s\n", strerror(errno), msg);
  else
    sprintf(new_msg, "note: %s\n", msg);

  va_start(fmtargs, msg);
  vfprintf(stderr, new_msg, fmtargs);
  va_end(fmtargs);
}


long long
nchoosek(int n, int k)
{
  long long res = 1;
  int i;
  for (i = 0; i < k; i++)
    res *= (n - i);
  for (i = 0; i < k; i++)
    res /= (i + 1);
  return res;
}


double
log_nchoosek(int n, int k)
{
  double res = 0.0;
  int i;
  for (i = 0; i < k; i++)
    res += log(n - i) - log(i + 1);
  return res;
}

void
cat(FILE * src, FILE * dest)
{
  size_t buffer_size = 2046;
  char buffer[buffer_size];
  size_t read;
  bool ends_in_newline = true;
  while ((read = fread(buffer, 1, buffer_size - 1, src))) {
    buffer[read] = '\0';
    fprintf(dest, "%s", buffer);
    if (buffer[read - 1] == '\n') {
      ends_in_newline = true;
    } else {
      ends_in_newline = false;
    }
  }
  if (!ends_in_newline) {
    fprintf(dest,"\n");
  }
}
