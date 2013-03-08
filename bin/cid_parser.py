#!/usr/bin/env python

import sys

def cid_parser(n):
    d = {}
    a = n.split(':')
    if len(a) < 5:
        sys.exit('could not parse read name: ' + n)
    d['rgid'] = a[0]
    d['clone_id'] = a[1]
    if a[2] in '12':
        d['paired'] = True
        d['nip'] = int(a[2]) - 1
        d['len'] = int(a[3 + d['nip']])
        d['mp_len'] = int(a[4 - d['nip']])
    else:
        d['paired'] = False
        d['len'] = int(a[3])
    if a[6] == '1':
        d['orig_seq'] = a[7]
        d['orig_qual'] = ':'.join(a[8:])[:d['len']]
    return d

def remove_seq_from_name(n):
    a = n.split(':')
    if a[6] != '1':
        return n
    if a[2] in '12':
        l = int(a[2 + int(a[2])])
    else:
        l = int(a[3])
    return (':'.join(a[:6])
            + ':0:'
            + ':'.join(a[7:])[(2 * (l + 1)):])
