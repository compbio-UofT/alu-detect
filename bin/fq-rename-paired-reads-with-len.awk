#!/usr/bin/awk -f

{
  line[0] = $0;
  for (k = 1; k < 8; k++) {
    getline;
    line[k] = $0;
  }

  if (length(line[1]) != length(line[3])) {
    printf("invalid read entry:\n%s\n%s\n%s\n%s\n",
      line[0], line[1], line[2], line[3]);
    exit(1);
  }
  if (length(line[5]) != length(line[7])) {
    printf("invalid read entry:\n%s\n%s\n%s\n%s\n",
      line[4], line[5], line[6], line[7]);
    exit(1);
  }

  line[0] = sprintf("@%09d:%d:%d:%d:0:%s", n, 1, length(line[1]), length(line[5]), substr(line[0], 2));
  line[4] = sprintf("@%09d:%d:%d:%d:0:%s", n, 2, length(line[1]), length(line[5]), substr(line[4], 2));

  for (k = 0; k < 8; k++) {
    print line[k];
  }

  n++;
}
