#!/usr/bin/awk -f

BEGIN {
  if (minlen == 0) minlen = 35;
}

{
  line[0] = $0;
  for (k = 1; k < 4; k++) {
    getline;
    line[k] = $0;
  }

  len = length(line[1]);

  if (length(line[3]) != len) {
    printf("invalid read entry:\n%s\n%s\n%s\n%s\n",
      line[0], line[1], line[2], line[3]);
    exit(1);
  }

  if (len < minlen)
    next;

  for (k = 0; k < 4; k++) {
    print line[k];
  }
}
