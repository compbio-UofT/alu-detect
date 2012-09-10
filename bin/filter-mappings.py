#!/usr/bin/env python2.6+

import sys
import argparse
import functools
from my import *
import my_sam as sam

warned_paired = False
default_min_mqv = 5
default_min_tail_match_len = 5
default_min_tail_insert_size = 15


def parse_bitmask(s):
    res = [0, 0]
    a = s.split('/')
    if len(a) == 1 or len(a) == 2:
        res[0] = int(a[0], 0)
        if len(a) == 2:
            res[1] = int(a[1], 0)
    else:
        crash('invalid bitmask specification: ' + s)
    if res[0] & res[1] != 0:
        crash('invalid bitmask specification: some bits must be both set and unset: ' + s)
    return res

def parse_filter(s):
    d = {}
    [conditions, dest_file] = e.split(':', 1)
    cond_list = conditions.split(',')
    if sam.pairing['paired'] == 0:
        if len(cond_list) == 1 or len(cond_list) == 2:
            [d['must_have'], d['must_not_have']] = parse_bitmask(cond_list[0])
            if len(cond_list) == 2:
                d['stop_on_hit'] = not (int(cond_list[1]) == 0)
            else:
                d['stop_on_hit'] = True
        else:
            crash('invalid condition specification: ' + conditions)
    else:
        if len(cond_list) == 2 or len(cond_list) == 3:
            d['must_have'] = [0, 0]
            d['must_not_have'] = [0, 0]
            [d['must_have'][0], d['must_not_have'][0]] = parse_bitmask(cond_list[0])
            [d['must_have'][1], d['must_not_have'][1]] = parse_bitmask(cond_list[1])
            if len(cond_list) == 3:
                d['stop_on_hit'] = not (int(cond_list[2]) == 0)
            else:
                d['stop_on_hit'] = True
        else:
            crash('invalid condition specification: ' + conditions)
    d['dest_file'] = dest_file
    return d


parser = argparse.ArgumentParser(description=sys.argv[0])
parser.add_argument('-v', '--verbose', action='append_const', const=1, default=[], dest='verbose')
parser.add_argument('-p', '--pairing', action='store', default='paired=0', dest='pairing')
parser.add_argument('--cid-parser', action='store', dest='cid_parser')
parser.add_argument('--check-pairing', action='store_true', default=False, dest='check_pairing')
parser.add_argument('--check-num-mappings', action='store_true', default=False, dest='check_num_mappings')
parser.add_argument('-f', '--filter', action='append', dest='filter', required=True)
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

if sam.pairing['paired'] == 1:
    args.check_num_mappings = True

if args.input:
    in_fd = gzopen(args.input)
else:
    in_fd = gzopen('-')

filters = []
for e in args.filter:
    d = parse_filter(e)
    d['file_object'] = open_file_or_fd(d['dest_file'], 'w')
    filters.append(d)
note('using filters: ' + str(filters), 1)

for ms in sam.get_mapping_set_gen(in_fd,
                                  header_line_hook=sys.stdout.write,
                                  mapping_parser=mapping_parser,
                                  check_pairing=args.check_pairing,
                                  check_num_mappings=args.check_num_mappings):
    note('ms: ' + str(ms), 2)

    for d in filters:
        if sam.pairing['paired'] == 0:
            if ((ms[0]['flag'] & d['must_have']) == d['must_have']
                and (ms[0]['flag'] & d['must_not_have']) == 0):
                print >> d['file_object'], '\t'.join(ms[0]['mapping'])
                if d['stop_on_hit']:
                    break
        else:
            if ((ms[0]['flag'] & d['must_have'][0]) == d['must_have'][0]
                and (ms[0]['flag'] & d['must_not_have'][0]) == 0
                and (ms[1]['flag'] & d['must_have'][1]) == d['must_have'][1]
                and (ms[1]['flag'] & d['must_not_have'][1]) == 0):
                print >> d['file_object'], '\t'.join(ms[0]['mapping'])
                print >> d['file_object'], '\t'.join(ms[1]['mapping'])
                if d['stop_on_hit']:
                    break

for d in filters:
    d['file_object'].close()
gzclose(in_fd)
