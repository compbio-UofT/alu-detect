#!/usr/bin/env awk

BEGIN{i=0}
{
	if ($0 ~ "^>") {
		if (!(SEQ=="")) {
			print SEQ
		}
		SEQ=substr($0,2)"\t"
	} else {
		SEQ=SEQ$0
	}
}
END{print SEQ}
