#!/usr/bin/env python
import signal
signal.signal(signal.SIGPIPE, signal.SIG_DFL)

import sys
import argparse
import functools
from lib_alu_detect import *
import lib_alu_detect_sam as sam

cid_parser = None
incomplete_fd = None
save_rgid = False
rgid_dict = dict()

def print_read_from_mapping(d):
    seq = d['seq']
    qual = d['qual']
    tail=''
    if d['paired'] and d['name'] == d['clone_id']:
        tail = '/' + ['1', '2'][d['nip']]
    if sam.is_reversed(d):
        seq = sam.reverse_complement(seq)
        qual = sam.reverse(qual)
    print '\n'.join(['@' + d['name'] + tail, seq, '+' + d['rgid'], qual])

def print_missing_read(d):
    if incomplete_fd:
        print >> incomplete_fd, d['name']
    else:
        print '\n'.join(['@' + d['name'], '*', '+' + d['rgid'], '*'])


parser = argparse.ArgumentParser(description='Given SAM mappings, produce fastq entries.')
parser.add_argument('-v', '--verbose', action='append_const', const=1, default=[], dest='verbose')
parser.add_argument('-l', '--pairing-file', action='store', dest='pairing_file')
parser.add_argument('-i', '--incomplete', action='store', dest='incomplete')
parser.add_argument('-s', '--save-rgid', action='store_true', default=False, dest='save_rgid')
parser.add_argument('-g', '--default-rg', action='store', default='00', dest='default_rg')
parser.add_argument('--cid-parser', action='store', dest='cid_parser')
parser.add_argument('input', action='store', default='-', nargs='?')
args = parser.parse_args()
set_log_level(len(args.verbose))

if args.pairing_file:
    sam.set_pairing('paired=1')
    pairing_file_fd = gzopen(args.pairing_file)
    for l in pairing_file_fd:
        l = l.strip().split('\t')
        rgid_dict[l[0]] = l[1]
    gzclose(pairing_file_fd)
    note('using rgid_dict=' + str(rgid_dict), 1)
else:
    if args.save_rgid:
        crash('when saving rgid-s, pairing file must be given')
    sam.set_pairing('paired=0')

in_fd = gzopen(args.input)

if args.cid_parser:
    _temp = __import__(args.cid_parser, globals(), locals(), ['cid_parser'])
    cid_parser = _temp.cid_parser
else:
    cid_parser = sam.default_cid_parser
mapping_parser = functools.partial(sam.mapping_parser, cid_parser=cid_parser)

if args.incomplete:
    incomplete_fd = open(args.incomplete, 'w')
if args.save_rgid:
    save_rgid = True

for ms in sam.get_mapping_set_gen(in_fd,
                                  mapping_parser=mapping_parser,
                                  check_pairing=True, check_num_mappings=True):

    for m in ms:
        if save_rgid:
            if 'RG' in m:
                m['rgid'] = rgid_dict[m['RG']]
            else:
                m['rgid'] = rgid_dict[args.default_rg]
        else:
            m['rgid'] = ''

    if sam.pairing['paired'] == 0:
        if ((not ms[0]['mapped'] or 'H' not in ms[0]['cigar'])
            and ms[0]['seq'] != '*' and len(ms[0]['seq']) == len(ms[0]['qual'])):
            print_read_from_mapping(ms[0])
        else:
            print_missing_read(ms[0])

    else:
        if ((not ms[0]['mapped'] or 'H' not in ms[0]['cigar'])
            and ms[0]['seq'] != '*' and len(ms[0]['seq']) == len(ms[0]['qual'])
            and (not ms[1]['mapped'] or 'H' not in ms[1]['cigar'])
            and ms[1]['seq'] != '*' and len(ms[1]['seq']) == len(ms[1]['qual'])):
            print_read_from_mapping(ms[0])
            print_read_from_mapping(ms[1])
        else:
            print_missing_read(ms[0])
            print_missing_read(ms[1])

gzclose(in_fd)
if incomplete_fd:
    incomplete_fd.close()
