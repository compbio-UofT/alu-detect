#!/usr/bin/awk -f

BEGIN {
    FS="\t"
    OFS="\t"
}

function detect_input(s)
{
    if (substr(s,1,1) == "@" && !match(s, "\t"))
	ifq = "fq";
    else if (match(s, "\t")) {
	split(s, p, "\t");
	if (substr(p[2],1,2) == "zs")
	    ifq = "sfq";
	else
	    ifq = "tfq";
    } else {
	printf("could not detect input format! [%s]\n", s) > "/dev/stderr";
	exit(1);
    }
    if (verbose)
      printf("detected input format [%s]\n", ifq) > "/dev/stderr";
    if (ofq == "") {
	ofq = ifq;
        if (verbose)
	  printf("set output format to [%s]\n", ifq) > "/dev/stderr";
    }
}

function get_read(s)
{
    if (ifq == "") detect_input($0);
    if (ifq == "fq") {
	r[1] = s;
	for (k = 2; k <= 4; k++) {
	    getline;
	    r[k] = $0;
	}
	if (substr(r[1], 1, 1) != "@") {
	    printf("read name does not start with '@' [%s]\n", r[1]) > "/dev/stderr";
	    exit(1);
	}
	if (match(r[1], "\t")) {
	    printf("read name contains tabs! [%s]\n", r[1]) > "/dev/stderr";
	    exit(1);
	}
	if (substr(r[3], 1, 1) != "+") {
	    printf("comment does not start with '+' [%s:%s]\n", r[1], r[3]) > "/dev/stderr";
	    exit(1);
	}
	if (match(r[3], "\t")) {
	    printf("read comment contains tabs! [%s:%s]\n", r[1], r[3]) > "/dev/stderr";
	    exit(1);
	}
	r[1] = substr(r[1], 2);
	r[3] = substr(r[3], 2);
    } else if (ifq == "tfq" || ifq == "sfq") {
	cnt = split(s, r, "\t");
	if (cnt != 4) {
	    printf("read entry [%s] contains [%d] fields instead of 4!\n", s, cnt) > "/dev/stderr";
	    exit(1);
	}
	if (ifq == "sfq") {
	    for (k = 2; k <= 4; k++) {
		r[k] = substr(r[k], 6);
	    }
	}
    }
    return r[1] "\t" r[2] "\t" r[3] "\t" r[4];
}

function put_read(s)
{
    split(s, r, "\t");
    if (ofq == "fq") {
	print "@" r[1];
	print r[2];
	print "+" r[3];
	print r[4];
    } else if (ofq == "tfq") {
	print r[1], r[2], r[3], r[4];
    } else if (ofq == "sfq") {
	print r[1], "zs:Z:" r[2], "zc:Z:" r[3], "zq:Z:" r[4];
    } else {
	ofq="fq"
	put_read(s)
    }
}

function get_readpair(s)
{
    s1 = get_read(s)
    getline
    s2 = get_read($0)
    return s1 "\t" s2
}

function put_readpair(s)
{
    split(s, rp, "\t")
    put_read(rp[1] "\t" rp[2] "\t" rp[3] "\t" rp[4])
    put_read(rp[5] "\t" rp[6] "\t" rp[7] "\t" rp[8])
}
