#!/bin/bash
set -u

NCPU=${NCPU:-4}
if which pigz >/dev/null 2>&1 ; then
    ZIP=${ZIP:-"pigz -9 -p $NCPU"}
    UNZIP=${UNZIP:-"pigz -d -p $NCPU"}
else
    ZIP=${ZIP:-"gzip -9"}
    UNZIP=${UNZIP:-gunzip}
fi

if which pv >/dev/null 2>&1; then
    PV=${PV:-"pv -f -i 60"}
else
    PV=${PV:-cat}
fi

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
    DELIM="${DELIM-"\\t"}"
    NEWLINE="${NEWLINE-"\\n"}"
    OPTIND=1
    while getopts "d:n" OPT "$@" ; do
	case $OPT in
	    d)
		DELIM="$OPTARG"
		;;
            n)
                NEWLINE=""
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
    printf "$NEWLINE"
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

run_stage() {
    echo "--------------------" >&2
    make_note "Stage $STAGE_NUM: $STAGE_NAME"

    if [ "${STAGE_NUM%%.*}" -lt "$START_STAGE" ] ; then
	make_note "skipping"
	return
    elif [ "${STAGE_NUM%%.*}" -gt "$END_STAGE" ] ; then
	make_note "done because END_STAGE=$END_STAGE"
	exit
    fi

    for f in $INPUT_FILES ; do
	if [ ! -r $f ] ; then crash "$f missing" ; fi
    done

    NEED_TO_RUN=0
    for f in $OUTPUT_FILES ; do
	if [ ! -r "$f" ] ; then
	    make_note "$f does not exist"
	    NEED_TO_RUN=1
	    break
	fi
    done
    if [ $NEED_TO_RUN = 0 ] ; then
	for f in $OUTPUT_FILES ; do
	    for g in $INPUT_FILES ; do
		if [ "$g" -nt "$f" ] ; then
		    make_note "$g is newer than $f"
		    NEED_TO_RUN=1
		    break
		fi
	    done
	    if [ $NEED_TO_RUN = 1 ] ; then break ; fi
	done
    fi
    if [ $NEED_TO_RUN = 0 ] ; then
	make_note "$OUTPUT_FILES up to date"
    else
	make_note "need to run"
	typeset -f stage_command >&2
	if ask_confirmation ; then
	    make_note "starting"
	    force_errexit stage_command || { rm -f $OUTPUT_FILES ; crash "failed" ; }
	    make_note "done"
	else
	    make_note "skipping"
	fi
    fi

    if [ "$STATS_FILE" ] ; then
	NEED_TO_RUN=0
	if [ ! -r "$STATS_FILE" ] ; then
	    make_note "$STATS_FILE does not exist"
	    NEED_TO_RUN=1
	fi
	if [ $NEED_TO_RUN = 0 ] ; then
	    for g in $OUTPUT_FILES ; do
		if [ "$g" -nt "$STATS_FILE" ] ; then
		    make_note "$g is newer than $STATS_FILE"
		    NEED_TO_RUN=1
		    break
		fi
	    done
	fi
	if [ $NEED_TO_RUN = 0 ] ; then
	    make_note "using existing stats"
	else
	    make_note "collecting stats"
	    typeset -f stage_stats_command >&2
	    if ask_confirmation ; then
	    make_note "starting"
		force_errexit stage_stats_command >$STATS_FILE || { rm -f $STATS_FILE ; crash "failed" ; }
		make_note "done"
	    else
		make_note "skipping"
	    fi
	fi
	cat $STATS_FILE >&2
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
    CMD_OUTPUT=()
    for _i in $(seq 1 $CRT_CMD) ; do
        CMD_OUTPUT[$_i]=$(echo "$CMDS_RAW_OUTPUT" | grep "^$_i " | cut -d " " -f 2-)
    done
    echo "CMD_OUTPUT:${CMD_OUTPUT[@]}" >&2
}

#
# Check whether to proceed with action
#
ask_confirmation() {
    if [ "${CONF:-}" = a ] ; then return ; fi
    read -p "continue? ([y]es/[a]lways/[s]kip/[N]o) " CONF
    if [ "$CONF" != y -a "$CONF" != a -a "$CONF" != s ] ; then exit 1 ; fi
    if [ "$CONF" = s ] ; then return 1 ; fi
}

#
# Create a file from given prerequisites
#
gen_file() {
    suspend_xtrace
    GEN_FILE_SUFFIX="${GEN_FILE_SUFFIX:-.pgen}"
    PRINT_OUTPUT="${PRINT_OUTPUT:-}"

    # check prerequisites
    for f in "${INPUT_FILES[@]:+${INPUT_FILES[@]}}" ; do
        if [ ! -r "$f" -a ! -r "${f}${GEN_FILE_SUFFIX}" ] ; then
	    crash "$f[$GEN_FILE_SUFFIX] missing"
	fi
    done

    # recursively regenerate prerequisites, if necessary
    for f in "${INPUT_FILES[@]:+${INPUT_FILES[@]}}" ; do
	if [ -r "${f}${GEN_FILE_SUFFIX}" ] ; then
	    PRINT_OUTPUT= ./"${f}${GEN_FILE_SUFFIX}" || crash "./${f}${GEN_FILE_SUFFIX}: failed"
	fi
    done

    # check if we need to run
    NEED_TO_RUN=0
    if [ ! -r "$OUTPUT_FILE" ] ; then
	make_note "$OUTPUT_FILE does not exist"
	NEED_TO_RUN=1
    fi
    if [ $NEED_TO_RUN = 0 ] ; then
	for f in "${INPUT_FILES[@]:+${INPUT_FILES[@]}}" ; do
	    if [ "$f" -nt "$OUTPUT_FILE" ] ; then
		make_note "$f is newer than $OUTPUT_FILE"
		NEED_TO_RUN=1
		break
	    fi
	done
    fi
    if [ $NEED_TO_RUN = 0 ] ; then
	if [ "$OUTPUT_FILE$GEN_FILE_SUFFIX" -nt "$OUTPUT_FILE" ] ; then
	    make_note "$OUTPUT_FILE$GEN_FILE_SUFFIX newer than $OUTPUT_FILE"
	    NEED_TO_RUN=1
	fi
    fi

    if [ $NEED_TO_RUN = 1 ] ; then
	make_note "need to run"
	typeset -f COMMAND >&2
	if ask_confirmation ; then
	    make_note "starting"
	    if [ "$PRINT_OUTPUT" ]; then
		exec 3>&1
	    else
		exec 3>/dev/null
	    fi

	    ( set -o pipefail ; COMMAND 3>&- | tee-p "$OUTPUT_FILE" >&3 ; ) \
		|| { rm "$OUTPUT_FILE" ; crash "failed" ; }

	    exec 3>&-
	    make_note "done"
	else
	    make_note "skipping"
	    if [ "$PRINT_OUTPUT" ]; then
		cat "$OUTPUT_FILE"
	    fi
	fi
    else
	if [ "$PRINT_OUTPUT" ]; then
	    cat "$OUTPUT_FILE"
	fi
    fi
    restore_xtrace
}

#
# Get relative path to file 1 from dir 2
#
rel_path() {
    [ $# -ne 2 ] && { make_note "use: rel_path <dest_file> <src_dir>"; return 1; }
    dest="$(cd -P $(dirname "$1") ; pwd)"
    src="$(cd -P "$2" ; pwd)"
    dest_array=($(echo "$dest" | sed 's/\// /g'))
    src_array=($(echo "$src" | sed 's/\// /g'))
    i=0
    while true; do
	[ $i -ge ${#dest_array[@]} ] && break
	[ $i -ge ${#src_array[@]} ] && break
	[ ! ${dest_array[i]} = ${src_array[i]} ] && break
	i=$(($i + 1))
    done
    j=$i
    res=""
    while [ $j -lt ${#src_array[@]} ]; do
	res="${res}../"
	j=$(($j + 1))
    done
    while [ $i -lt ${#dest_array[@]} ]; do
	res="${res}${dest_array[$i]}/"
	i=$(($i + 1))
    done
    echo "${res}$(basename "$1")"
}

#
# execute a command and ignore a return code of 141 = SIGPIPE
#
no_sigpipe () {
    exec 3>&1
    local exit_code=$({ force_errexit "$@" >&3 3>&-; echo $?; } || true)
    exec 3>&-
    if [ $exit_code -eq 141 ]; then
	return 0
    else
	return $exit_code
    fi
}

#
# run in fresh bash, thus enabling -e
#
set_explicit_errtrap () {
    trap 'echo $0: line $LINENO: exit code $?' ERR
}

quote () { 
    echo \'${1//\'/\'\\\'\'}\'
}

force_errexit () {
#    typeset | grep -vE "^(UID|EUID|PPID|BASH[A-Z_]*|SHELLOPTS)=" | grep -n '^' >&2
    bash <(
	typeset | grep -vE "^(UID|EUID|PPID|BASH[A-Z_]*|SHELLOPTS)="
	echo "set_explicit_errtrap"
	echo "BASH_XTRACEFD=\${BASH_XTRACEFD:-}"
	echo "set -eEux -o pipefail"
	for arg in "$@"; do
	    echo -n "$(quote "$arg") ";
	done
    )
}

get_new_fd () {
    local i=3
    while [ $i -lt 256 ] && [ -e /proc/$BASHPID/fd/$i ]; do
	let i+=1
    done
    [ $i -lt 256 ] || crash "no fd is available"
    echo $i
}

double_pipe () {
#    make_note "double_pipe: using elements"
#    typeset -f main_pipe alt_pipe splitter joiner >&2
    (
	coproc { alt_pipe | drain; }
	exec 4>&${COPROC[1]}-
	exec 3<&${COPROC[0]}-

	splitter 3<&- \
	    | { main_pipe | drain; } 3<&- 4>&- \
	    | joiner 4>&- &

	exec 4>&-
	exec 3<&-
	wait
    ) 3<&- 4>&-
}

#
# Actually executed
#
set_explicit_errtrap
if [ "${BASH_XTRACEFD:-}" ]; then
    export BASH_XTRACEFD=$BASH_XTRACEFD
else
    [ ! -e /proc/$$/fd/4 ] || crash "fd 4 needs to be closed"
    exec 4>${XTRACE:-/dev/null}
    export BASH_XTRACEFD=4
fi
