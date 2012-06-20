#!/usr/bin/awk -f

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

  line[0] = sprintf("@%09d:%d:%s", n, nip + 1, substr(line[0], 2));

  for (k = 0; k < 4; k++) {
    print line[k];
  }

  if (nip == 1) {
    n++;
    nip = 0;
  } else {
    nip++;
  }
}
