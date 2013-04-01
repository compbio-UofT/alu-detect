#!/bin/bash
source lib.common.sh
set_explicit_errtrap

NCPU=${NCPU:-4}
if cmd_exists pigz; then
    ZIP=${ZIP:-"pigz -9 -p $NCPU"}
    UNZIP=${UNZIP:-"pigz -d -p $NCPU"}
else
    ZIP=${ZIP:-"gzip -9"}
    UNZIP=${UNZIP:-gunzip}
fi

if cmd_exists pv; then
    PV=${PV:-"pv -f -i 10"}
else
    PV=${PV:-cat}
fi

if [ ! "${BASH_XTRACEFD:-}" ]; then
    XTRACE=${XTRACE:-/dev/null}
    exec {BASH_XTRACEFD}>"$XTRACE"
fi
export BASH_XTRACEFD=$BASH_XTRACEFD

setup_extra_paths () {
    [ -r "$extra_paths_file" ] || return 0
    local oldifs=$IFS
    IFS=$'\n'
    local line
    while read -r line; do
	add_to_path "$line" PATH
    done <"$extra_paths_file"
    IFS=$oldifs
}

# once BASE_DIR is defined, set up variable names
set_global_var_names () {
    extra_paths_file=$BASE_DIR/data/extra-paths
    alus_fa=$BASE_DIR/data/alus.fa
    alus_no_polya_fa=$BASE_DIR/data/alus.hidden-polya.fa
}
if [ "${BASE_DIR:-}" ]; then
    set_global_var_names
    setup_extra_paths
fi

set_ref_var_names () {
    ref_fa=$BASE_DIR/data/ref.$1.fa
    ref_fai=$BASE_DIR/data/ref.$1.fa.fai
    ref_bt2_idx=$BASE_DIR/data/ref.$1
    ref_alus_bed=$BASE_DIR/data/ref.$1.alus.bed.gz
}

set_real_and_fake_ref_var_names () {
    set_ref_var_names $1
    real_ref_name=$1
    real_ref_fa=$ref_fa
    real_ref_fai=$ref_fai
    real_ref_bt2_idx=$ref_bt2_idx
    real_ref_alus_bed=$ref_alus_bed

    set_ref_var_names ${2:-fake_$1}
    fake_ref_name=${2:-fake_$1}
    fake_ref_fa=$ref_fa
    fake_ref_fai=$ref_fai
    fake_ref_bt2_idx=$ref_bt2_idx
    fake_ref_alus_bed=$ref_alus_bed

    real_ref_alus_tsd_bed=$BASE_DIR/data/ref.$real_ref_name.alus.tsd.bed.gz
    real_deletions_bed=$BASE_DIR/data/deletions.$real_ref_name.to.$fake_ref_name.bed
    real_alus_to_remove_bed=$BASE_DIR/data/alus.to-remove.$real_ref_name.to.$fake_ref_name.bed
    fake_targets_bed=$BASE_DIR/data/targets.$fake_ref_name.bed

    real_known_novel_alus_bed=$BASE_DIR/data/known-novel-alus.$real_ref_name.bed
    fake_known_novel_alus_bed=$BASE_DIR/data/known-novel-alus.$fake_ref_name.bed
}

set -x
