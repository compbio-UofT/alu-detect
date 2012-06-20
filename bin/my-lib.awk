#!/usr/bin/awk -f

function n_index(s, t, n) {
    start_idx = 1
    for (i = 0; i < n; i++) {
	k = index(substr(s, start_idx), t);
	#print "i= " i " k=" k > "/dev/fd/2"
	if (k == 0) return 0;
	start_idx += k
    }
    return start_idx - 1;
}
