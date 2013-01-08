#!/bin/bash extend-hashbang
#make -f

prefix := ${WORK_DIR}table.${NGS_NAME}
suffix := csv

files := ${foreach len,${LEN_VALS},${foreach supp,${SUPP_VALS},${prefix}.len${len}.supp${supp}.${suffix}}}
get_len = ${shell echo $@ | grep -o 'len[0-9]*' | cut -c 4-}
get_supp = ${shell echo $@ | grep -o 'supp[0-9]*' | cut -c 5-}

.PHONY : all test guard-%

all : guard-NGS_NAME \
	guard-LEN_VALS guard-SUPP_VALS guard-NULL_VALS guard-CI_LEN_VALS \
	guard-OUTPUT_FD guard-WORK_DIR guard-CALLS_TRUTH_PAIRS \
	${files}
	cat ${prefix}.len*.supp*.${suffix} >&${OUTPUT_FD}

test : ${files}
	@ echo files: ${files}
	@ echo CALLS_TRUTH_PAIRS: "${CALLS_TRUTH_PAIRS}"

guard-% :
	@ if [ -z "${${*}}" ]; then echo "$* not set"; exit 1; fi

${prefix}.%.${suffix} :
	@ echo $@: len=${get_len}, supp=${get_supp}
	LEN_VALS=${get_len} SUPP_VALS=${get_supp} run-filters >>$@
