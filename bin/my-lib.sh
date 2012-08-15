#!/bin/bash
#set -x
set -o nounset
#set -o errexit
BASH_XTRACEFD="${BASH_XTRACEFD:-}"

MY_LOG_LEVEL=0
make_note() {
    suspend_xtrace
    MSG_LEVEL=0
    OPTIND=1
    while getopts "l:" OPT "$@"; do
	#echo "got option: $OPT" 1>&2
	case $OPT in
	    l)
		MSG_LEVEL="$OPTARG"
		;;
	esac
    done
    shift $(($OPTIND - 1))
    if [ "$MSG_LEVEL" = 0 ] ; then
	echo "note: $@" 1>&2
    elif [ "$MSG_LEVEL" -le "$MY_LOG_LEVEL" ] ; then
	echo "note [$MSG_LEVEL]: $@" 1>&2
    fi
    restore_xtrace
}

crash() {
    suspend_xtrace
    echo "error: $@" 1>&2
    exit 1
}

find_my_name_and_dir() {
    suspend_xtrace
    if [ $# -ne 1 ] ; then crash "use: $0 \${BASH_SOURCE[0]}" ; fi
    SOURCE="$1"
    while [ -h "$SOURCE" ] ; do
	SOURCE="$(readlink "$SOURCE")"
    done
    MY_NAME="$(basename "$SOURCE")"
    MY_DIR="$( cd -P "$( dirname "$SOURCE" )" && pwd )"
    restore_xtrace
}

make_temp_dir() {
    suspend_xtrace
    TDIR="$(mktemp -qd -t "${1:-tmp}"-XXXXXXXXXX)"
    if [ -z "$TDIR" ] ; then crash "could not make temporary dir!" ; fi
    if [ ! -d "$TDIR" ] ; then crash "temporary dir not a dir!" ; fi
    if [ ! -w "$TDIR" ] ; then crash "temporary dir not writable!" ; fi
    restore_xtrace
}

print_with_delim() {
    suspend_xtrace
    DELIM="\t"
    OPTIND=1
    while getopts "d:" OPT "$@" ; do
	case $OPT in
	    d)
		DELIM="$OPTARG"
		;;
	esac
    done
    shift $(($OPTIND - 1))

    FIRST=1
    while [ $# -gt 0 ] ; do
	if [ $FIRST = 0 ] ; then
	    printf "$DELIM"
	fi
	printf '%s' "$1"
	FIRST=0
	shift
    done
    printf '\n'
    restore_xtrace
}

suspend_xtrace() {
    OPTS="$SHELLOPTS"
    set +x
    if [[ "$OPTS" =~ xtrace ]] ; then
	XTRACE=1
    else
	XTRACE=0
    fi
}

restore_xtrace() {
    if [ "$XTRACE" = 1 ] ; then
	set -x
    fi
}

#
# Run several commands in one pass over stdin, report output in CMD_OUTPUT array
#
run_cmds() {
    if [ $# -lt 1 ] ; then
        crash "run_cmds needs arguments"
    fi
    CMDS_STRING=
    CRT_CMD=1
    while [ $# -ge 2 ] ; do
        CMDS_STRING="$CMDS_STRING tee-p >(echo \"$CRT_CMD \$($1)\" >&3 ) |"
        CRT_CMD=$(($CRT_CMD + 1))
        shift
    done
    CMDS_STRING="$CMDS_STRING echo \"$CRT_CMD \$($1)\" >&3"
    echo "CMDS_STRING:$CMDS_STRING" >&2
    CMDS_RAW_OUTPUT=$( { eval $CMDS_STRING ; } 3>&1 )
    echo "CMDS_RAW_OUTPUT:$CMDS_RAW_OUTPUT" >&2
    for i in $(seq 1 $CRT_CMD) ; do
        CMD_OUTPUT[$i]=$(echo "$CMDS_RAW_OUTPUT" | grep "^$i " | cut -d " " -f 2-)
    done
}

#
# Check whether to proceed with action
#
ask_confirmation() {
    if [ "${CONF:-}" = a ] ; then return ; fi
    read -p "continue? ([y]es/[a]lways/[s]kip/[N]o) " CONF
    if [ "$CONF" != y -a "$CONF" != a -a "$CONF" != s ] ; then exit ; fi
    if [ "$CONF" = s ] ; then return 1 ; fi
}

#
# Create a file from given prerequisites
#
gen_file() {
    GEN_FILE_SUFFIX=${GEN_FILE_SUFFIX:-.pgen}

    # check prerequisites
    for f in "${INPUT_FILES[@]}" ; do
        if [ ! -r "$f" -a ! -r "${f}${GEN_FILE_SUFFIX}" ] ; then
	    crash "$f[$GEN_FILE_SUFFIX] missing"
	fi
    done

    # recursively regenerate prerequisites, if necessary
    for f in "${INPUT_FILES[@]}" ; do
	if [ -r "${f}${GEN_FILE_SUFFIX}" ] ; then
	    bash "${f}${GEN_FILE_SUFFIX}" >/dev/null
	fi
    done

    # check if we need to run
    NEED_TO_RUN=0
    if [ ! -r "$OUTPUT_FILE" ] ; then
	make_note "$OUTPUT_FILE does not exist"
	NEED_TO_RUN=1
    fi
    if [ $NEED_TO_RUN = 0 ] ; then
	for f in "${INPUT_FILES[@]}" ; do
	    if [ "$f" -nt "$OUTPUT_FILE" ] ; then
		make_note "$f is newer than $OUTPUT_FILE"
		NEED_TO_RUN=1
		break
	    fi
	done
    fi

    if [ $NEED_TO_RUN = 1 ] ; then
	make_note "need to run"
	typeset -f COMMAND >&2
	if ask_confirmation ; then
	    make_note "starting"
	    ( set -o pipefail ; COMMAND | tee-p "$OUTPUT_FILE" ; ) \
		|| { rm "$OUTPUT_FILE" ; crash "failed" ; }
	    make_note "done"
	else
	    make_note "skipping"
	    cat "$OUTPUT_FILE"
	fi
    else
	cat "$OUTPUT_FILE"
    fi
}
