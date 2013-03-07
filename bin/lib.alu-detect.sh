#!/bin/bash
. lib.common.sh
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

set -ux
