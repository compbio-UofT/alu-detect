#!/bin/bash

DIR="$( cd -P "$( dirname "${BASH_SOURCE[0]}" )" && pwd )"
if [[ ! -re "$DIR/settings/settings.hg18" ]]; then
	echo "Please run setup!"
	exit 1
fi
. $DIR/settings/settings.hg18
$DIR/bin/alu-detect
