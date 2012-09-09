#!/usr/bin/env python2.6+

import sys
import argparse
import functools
from my import *
import my_sam as sam

default_min_mqv = 5
default_min_len = 20
default_min_tail_match_len = 5
default_min_tail_insert_size = 15


parser = argparse.ArgumentParser(description=sys.argv[0])
parser.add_argument('-v', '--verbose', action='append_const', const=1, default=[], dest='verbose')
parser.add_argument('-p', '--pairing', action='store', default='paired=0', dest='pairing')
parser.add_argument('-q', '--min-mqv', action='store', type=int,
                    default=default_min_mqv, dest='min_mqv')
parser.add_argument('-l', '--min-len', action='store', type=int,
                    default=default_min_len, dest='min_len')
parser.add_argument('-m', '--min-tail-match-len', action='store', type=int,
                    default=default_min_tail_match_len, dest='min_tail_match_len')
parser.add_argument('-i', '--min-tail-insert-size', action='store', type=int,
                    default=default_min_tail_insert_size, dest='min_tail_insert_size')
parser.add_argument('--cid-parser', action='store', dest='cid_parser')
parser.add_argument('input', action='store', nargs='?')
args = parser.parse_args()
set_log_level(len(args.verbose))
sam.set_pairing(args.pairing)

if args.cid_parser:
    _temp = __import__(args.cid_parser, globals(), locals(), ['cid_parser'])
    cid_parser = _temp.cid_parser
else:
    cid_parser = sam.default_cid_parser
mapping_parser = functools.partial(sam.mapping_parser, cid_parser=cid_parser)

if args.input:
    in_fd = gzopen(args.input)
else:
    in_fd = gzopen('-')

for ms in sam.get_mapping_set_gen(in_fd,
                                  header_line_hook=sys.stdout.write,
                                  mapping_parser=mapping_parser):
    note('ms: ' + str(ms), 2)
    if len(ms) > 2:
        crash('ms too large')

    for d in ms:
        if d['len'] < args.min_len:
            d['flag'] |= 0x10000
        if d['mapped']:
            if d['mqv'] >= args.min_mqv:
                d['flag'] |= 0x1000
            if args.min_tail_insert_size == 0:
                continue
            ops = sam.parse_cigar_string(d['cigar'])
            tails = sam.get_tail_insert_size(ops, min_tail_match_len=args.min_tail_match_len)
            if tails[0] >= args.min_tail_insert_size:
                d['flag'] |= 0x2000
            if tails[1] >= args.min_tail_insert_size:
                d['flag'] |= 0x4000

    if (sam.pairing['paired'] == 1 and len(ms) == 2 and ms[0]['mapped'] and ms[1]['mapped']
        and sam.pair_concordant(ms)):
        ms[0]['flag'] |= 0x8000
        ms[1]['flag'] |= 0x8000

    for d in ms:
        d['mapping'][1] = str(d['flag'])
        print '\t'.join(d['mapping'])

gzclose(in_fd)
