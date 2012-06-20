#!/bin/bash

DIR="$( cd -P "$( dirname "${BASH_SOURCE[0]}" )" && pwd )"
if [[ ! -re "$DIR/settings/settings.hg18" ]]; then
	echo "Please run setup!"
	exit 1
fi
. $DIR/data/settings.hg18
#RGID=00 NCPU=24 ORIG_MAPPINGS=/data/matei/experiment/SRX079579-WXS-NA18507/map.sam.gz CONF=a $DIR/bin/alu-detect 2> 00.alu-detect.log
RGID=00 NCPU=24 ORIG_MAPPINGS=/data/matei/2012.05.28/sim2/mappings/map.sam.gz CONF=a $DIR/bin/alu-detect 2> 00.alu-detect.log
