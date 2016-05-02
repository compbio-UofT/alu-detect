#!/usr/bin/make -f

inputs := ${FAKE_CALLS} ${FAKE_CALLS_REF_ALU_BP} ${FAKE_CALLS_NEAR_ALUS} ${FAKE_CALLS_NEAR_ALUS_REF_ALU_BP} ${REAL_CALLS} ${REAL_CALLS_REF_ALU_BP}
prefix := ${WORK_DIR}table.${NGS_NAME}
suffix := csv

files := ${foreach len,${LEN_VALS},${foreach supp,${SUPP_VALS},${prefix}.len${len}.supp${supp}.${suffix}}}
get_len = ${shell echo $@ | grep -o 'len[0-9]*' | cut -c 4-}
get_supp = ${shell echo $@ | grep -o 'supp[0-9]*' | cut -c 5-}

.PHONY : all test guard-%

all : guard-NGS_NAME \
	guard-FAKE_CALLS guard-FAKE_CALLS_REF_ALU_BP \
	guard-FAKE_CALLS_NEAR_ALUS guard-FAKE_CALLS_NEAR_ALUS_REF_ALU_BP \
	guard-REAL_CALLS guard-REAL_CALLS_REF_ALU_BP \
	guard-LEN_VALS guard-SUPP_VALS guard-NULL_VALS guard-CI_LEN_VALS \
	guard-TARGETS guard-TARGETS_IN_ALUS guard-STANDARD_ALUS guard-KNOWN_ALUS \
	guard-OUTPUT_FD guard-WORK_DIR \
	${files}
	cat ${prefix}.len*.supp*.${suffix} >&${OUTPUT_FD}

test : ${files}
	@ echo files: ${files}

guard-% :
	@ if [ -z "${${*}}" ]; then echo "$* not set"; exit 1; fi

${prefix}.%.${suffix} :
	@ echo $@: len=${get_len}, supp=${get_supp}
#	LEN_VALS=${get_len} SUPP_VALS=${get_supp} run-filters ${FAKE_CALLS} ${FAKE_CALLS_NEAR_ALUS} ${REAL_CALLS} >$@
	LEN_VALS=${get_len} SUPP_VALS=${get_supp} run-filters ${FAKE_CALLS_REF_ALU_BP} ${FAKE_CALLS_NEAR_ALUS_REF_ALU_BP} ${REAL_CALLS_REF_ALU_BP} | tawk '{$$5="ref-alus"; print}' >>$@

