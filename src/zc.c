#include <stdio.h>
#include <stdlib.h>
#include <zlib.h>
#include <string.h>
#include <errno.h>

#define BLOCK_SIZE (1 << 20)

char buffer[BLOCK_SIZE];

void
process_file(gzFile f)
{
  int bytes_in;

  while ((bytes_in = gzread(f, buffer, BLOCK_SIZE)) > 0) {
    fwrite(buffer, bytes_in, sizeof(char), stdout);
  }
}


int
main(int argc, char * argv[])
{
  int i;
  gzFile f;
  int seen_stdin = 0;

  if (argc < 2) {
    f = gzdopen(fileno(stdin), "r");
    process_file(f);
    gzclose(f);
  } else {
    for (i = 1; i < argc; i++) {
      if (strncmp(argv[i], "-", 2) == 0) {
	if (seen_stdin) {
	  // ignore it
	  //fprintf(stderr, "cannot cat stdin twice\n");
	  //exit(1);
	} else {
	  seen_stdin = 1;
	  f = gzdopen(fileno(stdin), "r");
	  process_file(f);
	  gzclose(f);
	}
      } else {
	f = gzopen(argv[i], "r");
	if (f == NULL) {
	  fprintf(stderr, "error opening file [%s]: %s\n", argv[i], strerror(errno));
	}
	process_file(f);
	gzclose(f);
      }
    }
  }

  fclose(stdout);
  return 0;
}
