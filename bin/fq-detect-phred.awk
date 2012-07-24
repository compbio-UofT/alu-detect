#!/usr/bin/awk -f

function detect_phred(s) {
    # must autodetect input_phred
    # this is done on the first qual string
    # if autodetection fails, crash
    h[1] = 33;
    h[2] = 64;
    n_h = 2;
    min_mqv = -5;
    max_mqv = 45;
    for (i = 1; i <= n_h; i++) {
	works[i] = 1;
	for (j = 1; j <= length(s); j++) {
	    c = substr(s, j, 1)
	    if (ord[c] - h[i] < min_mqv || ord[c] - h[i] > max_mqv) {
		print "character [" c "] with ord [" ord[c] "] contradicts phred [" h[i] "]" >"/dev/stderr";
		works[i] = 0;
		break;
	    }
	}
    }
    k = 0;
    for (i = 1; i <= n_h; i++) {
	if (works[i] > 0) {
	    if (k == 0) {
		k = i;
	    } else {
		print "could not detect phred: " h[k] " and " h[i] " are both possible" >"/dev/stderr";
		exit(1);
	    }
	}
    }
    print "autodetected phred [" h[k] "]" >"/dev/stderr";
    return h[k]
}
