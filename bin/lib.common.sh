#!/bin/bash

crash () {
    echo "error: $@" >&2
    exit 1
}

set_explicit_errtrap () {
    trap 'echo "$0: line $LINENO: exit code $?" >&2' ERR
}

cmd_exists () {
    which "$1" >/dev/null 2>&1
}

check_files_readable () {
    for f in "$@"; do
	[ -r "$f" ] || crash "file not readable [$f]"
    done
}

check_files_not_exist () {
    for f in "$@"; do
	[ ! -e "$f" ] || crash "file exists [$f]"
    done
}

quote () { 
    echo \'${1//\'/\'\\\'\'}\'
}

# use: add_to_path /a/new/dir AWKPATH
add_to_path () {
    local dir=$1
    local var_name=$2
    [[ ! $(eval echo \$$var_name) =~ "$dir" ]] || return 0
    eval $var_name=\"$dir\${$var_name:+:}\$$var_name\"
}

get_unused_fd () {
    local i=3
    while [ $i -lt 256 ] && [ -e /proc/$1/fd/$i ]; do
	let i+=1
    done
    [ $i -lt 256 ] || crash "no fd is available"
    echo $i
}

my_log_level=0
make_note () {
    suspend_xtrace
    local msg_level=0
    OPTIND=1
    while getopts "l:" OPT "$@"; do
	case $OPT in
	    l)
		msg_level=$OPTARG
		;;
	esac
    done
    shift $(($OPTIND - 1))
    if [ "$msg_level" = 0 ] ; then
	echo "note: $@" >&2
    elif [ "$msg_level" -le "$my_log_level" ] ; then
	echo "note [$msg_level]: $@" >&2
    fi
    restore_xtrace
}

find_my_name_and_dir () {
    suspend_xtrace
    [ $# -eq 1 ] || crash "use: $0 \${BASH_SOURCE[0]}"
    local source=$1
    while [ -h "$source" ] ; do
	source=$(readlink "$source")
    done
    MY_NAME=$(basename "$source")
    MY_DIR=$(cd -P "$( dirname "$source" )" && pwd)
    restore_xtrace
}

make_temp_dir () {
    suspend_xtrace
    TDIR="$(mktemp -qd -t "${1:-tmp}"-XXXXXXXXXX)"
    [ ! -z "$TDIR" ] || crash "could not make temporary dir!"
    [ -d "$TDIR" ] || crash "temporary dir not a dir!"
    [ -w "$TDIR" ] || crash "temporary dir not writable!"
    restore_xtrace
}

suspend_xtrace () {
    local opts=$SHELLOPTS
    set +x
    if [[ "$opts" =~ xtrace ]] ; then
	XTRACE=1
    else
	XTRACE=0
    fi
}

restore_xtrace () {
    if [ "$XTRACE" = 1 ] ; then
	set -x
    fi
}

run_stage () {
    echo "--------------------" >&2
    make_note "Stage $STAGE_NUM: $STAGE_NAME"

    if [ "${STAGE_NUM%%.*}" -lt "$START_STAGE" ] ; then
	make_note "skipping"
	return
    elif [ "${STAGE_NUM%%.*}" -gt "$END_STAGE" ] ; then
	make_note "done because END_STAGE=$END_STAGE"
	exit
    fi

    for f in $INPUT_FILES; do
	[ -r $f ] || crash "$f missing"
    done

    local need_to_run=
    for f in $OUTPUT_FILES; do
	if [ ! -r "$f" ] ; then
	    make_note "$f does not exist"
	    need_to_run=1
	    break
	fi
    done
    if [ ! $need_to_run ]; then
	for f in $OUTPUT_FILES; do
	    for g in $INPUT_FILES; do
		if [ "$g" -nt "$f" ]; then
		    make_note "$g is newer than $f"
		    need_to_run=1
		    break
		fi
	    done
	    [ ! $need_to_run ] || break
	done
    fi
    if [ ! $need_to_run ]; then
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
}

#
# Run several commands in one pass over stdin, report output in CMD_OUTPUT array
#
run_cmds() {
    [ $# -ge 1 ] || crash "run_cmds needs arguments"
    local pid=$BASHPID
    local fd=$(get_unused_fd $pid)
    local cmds_string=
    local crt_cmd=1
    while [ $# -ge 2 ]; do
        cmds_string="$cmds_string tee-p >(echo \"$crt_cmd \$($1)\" >&$fd) |"
        let crt_cmd+=1
        shift
    done
    cmds_string="$cmds_string echo \"$crt_cmd \$($1)\" >&$fd"
    echo "cmds_string:$cmds_string" >&2
    local cmds_raw_output=$( { eval $cmds_string ; } {fd}>&1 )
    echo "cmds_raw_output:$cmds_raw_output" >&2
    CMD_OUTPUT=()
    local i
    for i in $(seq 1 $crt_cmd); do
        CMD_OUTPUT[$i]=$(echo "$cmds_raw_output" | grep "^$i " | cut -d " " -f 2-)
    done
    echo "CMD_OUTPUT:${CMD_OUTPUT[@]}" >&2
}

#
# Check whether to proceed with action
#
ask_confirmation () {
    [ ! "${CONF:-}" = a ] || return 0
    read -p "${1:-"continue?"} ([y]es/[a]lways/[s]kip/[N]o) " CONF
    ([ "$CONF" = y ] || [ "$CONF" = a ] || [ "$CONF" = s ]) || exit 1
    [ ! "$CONF" = s ] || return 1
}

#
# Create a file from given prerequisites
#
gen_file () {
    suspend_xtrace
    GEN_FILE_SUFFIX=${GEN_FILE_SUFFIX:-.pgen}
    PRINT_OUTPUT=${PRINT_OUTPUT:-}

    # check prerequisites
    for f in $INPUT_FILES; do
        ([ -r "$f" ] || [ -r "${f}${GEN_FILE_SUFFIX}" ]) ||
	crash "$f[$GEN_FILE_SUFFIX] missing"
    done

    # recursively regenerate prerequisites, if necessary
    for f in $INPUT_FILES; do
	if [ -r "${f}${GEN_FILE_SUFFIX}" ] ; then
	    PRINT_OUTPUT= ./"${f}${GEN_FILE_SUFFIX}" ||
	    crash "./${f}${GEN_FILE_SUFFIX}: failed"
	fi
    done

    # check if we need to run
    local need_to_run=
    if [ ! -r "$OUTPUT_FILE" ]; then
	make_note "$OUTPUT_FILE does not exist"
	need_to_run=1
    fi
    if [ ! $need_to_run ]; then
	for f in $INPUT_FILES; do
	    if [ "$f" -nt "$OUTPUT_FILE" ]; then
		make_note "$f is newer than $OUTPUT_FILE"
		need_to_run=1
		break
	    fi
	done
    fi
    if [ ! $need_to_run ] ; then
	if [ "$OUTPUT_FILE$GEN_FILE_SUFFIX" -nt "$OUTPUT_FILE" ] ; then
	    make_note "$OUTPUT_FILE$GEN_FILE_SUFFIX newer than $OUTPUT_FILE"
	    need_to_run=1
	fi
    fi

    if [ $need_to_run = 1 ] ; then
	make_note "need to run"
	typeset -f COMMAND >&2
	if ask_confirmation ; then
	    make_note "starting"
	    (
		set -o pipefail
		COMMAND |
		{
		    if [ "$PRINT_OUTPUT" ]; then
			exec cat >"$OUTPUT_FILE"
		    else
			exec tee-p "$OUTPUT_FILE"
		    fi
		}
	    ) ||
	    { rm "$OUTPUT_FILE" ; crash "failed" ; }
	    make_note "done"
	else
	    make_note "skipping"
	    [ ! "$PRINT_OUTPUT" ] || cat "$OUTPUT_FILE"
	fi
    else
	[ ! "$PRINT_OUTPUT" ] || cat "$OUTPUT_FILE"
    fi
    restore_xtrace
}

#
# Get relative path to file 1 from dir 2
#
rel_path () {
    [ $# -ne 2 ] && { make_note "use: rel_path <dest_file> <src_dir>"; return 1; }
    local dest=$(cd -P $(dirname "$1") ; pwd)
    local src=$(cd -P "$2" ; pwd)
    local dest_array=($(echo "$dest" | sed 's/\// /g'))
    local src_array=($(echo "$src" | sed 's/\// /g'))
    local i=0
    while true; do
	[ $i -lt ${#dest_array[@]} ] || break
	[ $i -lt ${#src_array[@]} ] || break
	[ ${dest_array[i]} = ${src_array[i]} ] || break
	let i+=1
    done
    local j=$i
    local res=
    while [ $j -lt ${#src_array[@]} ]; do
	res="${res}../"
	let j+=1
    done
    while [ $i -lt ${#dest_array[@]} ]; do
	res="${res}${dest_array[$i]}/"
	let i+=1
    done
    echo "${res}$(basename "$1")"
}

#
# execute a command and ignore a return code of 141 = SIGPIPE
#
no_sigpipe () {
    local fd
    exec {fd}>&1
    local exit_code=$({ force_errexit "$@" >&$fd {fd}>&-; echo $?; } || true)
    exec {fd}>&-
    if [ $exit_code -eq 141 ]; then
	return 0
    else
	return $exit_code
    fi
}

#
# run in fresh bash, thus enabling -e
#
force_errexit () {
    bash <(
	{
	    typeset | grep -vE "^(UID|EUID|PPID|BASH[A-Z_]*|SHELLOPTS)="
	    echo "set_explicit_errtrap"
	    echo "BASH_XTRACEFD=\${BASH_XTRACEFD:-}"
	    echo "set -eEux -o pipefail"
	    for arg in "$@"; do
		echo -n "$(quote "$arg") ";
	    done 
	} |
	tee >(grep -n '^' >&$BASH_XTRACEFD)
    )
}

#
# stdin -> { double pipe } -> stdout
#
# main_pipe, alt_pipe: bash func processing stdin -> stdout; can be undefined
# splitter: bash function stdin -> { stdout, $out_fd }; can be undefined
# joiner: bash func { stdin, $in_fd } -> stdout; must be defined
#
double_pipe () (
#    make_note "double_pipe: using elements"
#    typeset -f main_pipe alt_pipe splitter joiner >&2
    (set +o pipefail; type -t splitter | grep -q function) ||
    {
	make_note "$0: undefined splitter; using default"
	splitter () {
	    exec tee-p >(exec cat >$out_fd)
	}
    }
    (set +o pipefail; type -t joiner | grep -q function) || crash "$0: joiner undefined"

    pipe_drain_wrapper () {
	if (set +o pipefail; type -t $1 | grep -q function); then
#	    make_note "$1 not empty"
	    $1 | drain
	else
#	    make_note "$1 empty"
	    exec drain
	fi
    }

    coproc pipe_drain_wrapper alt_pipe
    local alt_pid=$COPROC_PID
    local out_fd
    local in_fd
    exec {out_fd}>&${COPROC[1]}-
#   make_note "using out_fd: $out_fd"
    exec {in_fd}<&${COPROC[0]}-
#   make_note "using in_fd: $in_fd"

    out_fd=$out_fd splitter {in_fd}<&- |
    pipe_drain_wrapper main_pipe {in_fd}<&- {out_fd}>&- |
    in_fd=$in_fd joiner {out_fd}>&- &
    local main_pid=$!

    exec {out_fd}>&-
    exec {in_fd}<&-
    wait $main_pid
    wait $alt_pid
)
