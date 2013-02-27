#!/usr/bin/env python

import sys

def cid_parser(n):
    d = {}
    a = n.split(':')
    if len(a) < 4:
        sys.exit('could not parse read name: ' + n)
    d['clone_id'] = a[0]
    if a[1] in '12':
        d['paired'] = True
        d['nip'] = int(a[1]) - 1
        d['len'] = int(a[2 + d['nip']])
        d['mp_len'] = int(a[3 - d['nip']])
    else:
        d['paired'] = False
        d['len'] = int(a[2])
    if a[5] == '1':
        d['orig_seq'] = a[6]
        d['orig_qual'] = ':'.join(a[7:])[:d['len']]
    return d

def remove_seq_from_name(n):
    a = n.split(':')
    if a[5] != '1':
        return n
    if a[1] in '12':
        l = int(a[1 + int(a[1])])
    else:
        l = int(a[2])
    return (':'.join(a[:5])
            + ':0:'
            + ':'.join(a[6:])[(2 * (l + 1)):])
