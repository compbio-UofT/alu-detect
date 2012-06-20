#!/bin/bash
#set -x
set -o nounset
#set -o errexit
BASH_XTRACEFD="${BASH_XTRACEFD:-}"

MY_LOG_LEVEL=0
make_note() {
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
}

crash() {
    echo "error: $@" 1>&2
    exit 1
}

find_my_name_and_dir() {
    if [ $# -ne 1 ] ; then crash "use: $0 \${BASH_SOURCE[0]}" ; fi
    SOURCE="$1"
    while [ -h "$SOURCE" ] ; do
	SOURCE="$(readlink "$SOURCE")"
    done
    MY_NAME="$(basename "$SOURCE")"
    MY_DIR="$( cd -P "$( dirname "$SOURCE" )" && pwd )"
}

make_temp_dir() {
    TDIR="$(mktemp -qd -t "${1:-tmp}"-XXXXXXXXXX)"
    if [ -z "$TDIR" ] ; then crash "could not make temporary dir!" ; fi
    if [ ! -d "$TDIR" ] ; then crash "temporary dir not a dir!" ; fi
    if [ ! -w "$TDIR" ] ; then crash "temporary dir not writable!" ; fi
}

print_with_delim() {
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
}
