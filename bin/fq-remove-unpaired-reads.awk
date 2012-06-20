#!/usr/bin/awk -f

BEGIN {
  if (tail_regexp == "") tail_regexp = "/[12]$";
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

  tail_start = match(line[0], tail_regexp);
  if (tail_start == 0) {
    printf("read name [%s] does not match tail_regexp [%s]\n", line[0], tail_regexp);
    exit(1);
  }
  name = substr(line[0], 1, tail_start - 1);

  if (name == last_name) {
    for (k = 0; k < 4; k++) {
      print last_line[k];
    }
    for (k = 0; k < 4; k++) {
      print line[k];
    }
    last_name = "";
  } else {
    for (k = 0; k < 4; k++) {
      last_line[k] = line[k];
    }
    last_name = name;
  }
}
