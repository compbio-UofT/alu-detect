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

  trim_start = match(line[3], "B*$");
  if (trim_start > 0) {
    if (trim_start == 1) trim_start = 2;
    line[1] = substr(line[1], 1, trim_start - 1);
    line[3] = substr(line[3], 1, trim_start - 1);
  }

  for (k = 0; k < 4; k++) {
    print line[k];
  }
}
