#!/usr/bin/env python

import sys
import operator
from lib_alu_detect import *

IS_PAIRED = 0x1
ALL_PROPERLY_ALIGNED = 0x2
IS_UNMAPPED = 0x4
MP_UNMAPPED = 0x8
IS_REVERSED = 0x10
MP_REVERSED = 0x20
IS_FIRST = 0x40
IS_LAST = 0x80
SECONDARY_ALIGNMENT = 0x100
FAILED_QC = 0x200
IS_DUPLICATE = 0x400

cigar_qr_op = 'HSIM=X'
cigar_db_op = 'NDM=X'
cigar_clip_op = 'HS'

rc_dict = {'A': 'T',
           'C': 'G',
           'G': 'C',
           'T': 'A',
           'R': 'Y',
           'Y': 'R',
           'S': 'W',
           'W': 'S',
           'K': 'M',
           'M': 'K',
           'B': 'V',
           'V': 'B',
           'D': 'H',
           'H': 'D',
           'N': 'N'}

note_debug_level = 3

def is_paired(m):
    return flag(m, IS_PAIRED)
def all_properly_aligned(m):
    return flag(m, ALL_PROPERLY_ALIGNED)
def is_unmapped(m):
    return flag(m, IS_UNMAPPED)
def mp_unmapped(m):
    return flag(m, MP_UNMAPPED)
def is_reversed(m):
    return flag(m, IS_REVERSED)
def mp_reversed(m):
    return flag(m, MP_REVERSED)
def is_first(m):
    return flag(m, IS_FIRST)
def is_last(m):
    return flag(m, IS_LAST)
def secondary_alignment(m):
    return flag(m, SECONDARY_ALIGNMENT)
def failed_qc(m):
    return flag(m, FAILED_QC)
def is_duplicate(m):
    return flag(m, IS_DUPLICATE)

def flag(m, f):
    if isinstance(m, dict):
        val = m['flag']
    else:
        val = int(m[1])
    return val & f != 0

def default_cid_parser(n):
    d = {}
    if n[-2] == '/' and n[-1] in '12':
        d['clone_id'] = n[:-2]
    else:
        d['clone_id'] = n
    return d

def mapping_parser(m, cid_parser=default_cid_parser):
    if isinstance(m, str):
        m = m.strip().split('\t')
    d = {}
    d['mapping'] = m
    d['name'] = m[0]         # read name
    d['flag'] = int(m[1])    # flags
    d['paired'] = is_paired(m)
    d['mapped'] = not is_unmapped(m)
    d['chr'] = m[2]         # chr
    d['pos'] = int(m[3])    # pos
    d['mqv'] = int(m[4])    # mapq
    d['cigar'] = m[5]         # cigar string
    if d['paired']:
        d['mp_chr'] = m[6]         # mp chr
        if d['mp_chr'] == '=':
            d['mp_chr'] = d['chr']
        d['mp_pos'] = int(m[7])    # mp pos
        d['t_len'] = int(m[8])    # template len
    d['seq'] = m[9]         # sequence
    d['qual'] = m[10]       # qual

    if d['seq'] != '*':
        d['len'] = len(d['seq'])
    if d['paired']:
        d['nip'] = int(is_last(m))
    if d['mapped']:
        d['st'] = int(is_reversed(m))

    i = 11
    while i < len(m):
        key = m[i][0:2]
        if m[i][3] == 'i':
            d[key] = int(m[i][5:])
        else:
            d[key] = m[i][5:]
        i += 1
    d.update(cid_parser(m[0]))
    return d

def parse_cigar_string(s):
    note('parsing cigar [%s]' % s, note_debug_level)
    res = []
    crt_len = ''
    i = 0
    while i < len(s):
        if str.isdigit(s[i]):
            crt_len += s[i]
        else:
            res.append([s[i], int(crt_len)])
            crt_len = ''
        i += 1
    return res

def get_len_from_cigar(s):
    op = parse_cigar_string(s)
    qr_len = sum(map(operator.itemgetter(1), filter(lambda e: e[0] in cigar_qr_op, op)))
    db_len = sum(map(operator.itemgetter(1), filter(lambda e: e[0] in cigar_db_op, op)))
    clip_len = [0, 0]
    i = 0
    while i < len(op) and op[i][0] in cigar_clip_op:
        clip_len[0] += op[i][1]
        i += 1
    i = len(op) - 1
    while i >= 0 and op[i][0] in cigar_clip_op:
        clip_len[1] += op[i][1]
        i -= 1

    return [qr_len, db_len, clip_len]

def get_pos(m):
    if isinstance(m, dict):
        cigar = m['cigar']
        left_pos = m['pos']
    else:
        cigar = m[5]
        left_pos = int(m[3])
    [qr_len, db_len, clip_len] = get_len_from_cigar(cigar)
    pos = [left_pos, left_pos + db_len - 1]
    return [pos, clip_len]

pairing = {
    'paired': 1,
    'st_diff': 1,
    'min': 100,
    'max': 1000,
    'mean': 300,
    'stddev': 60,
    'r1_len': 100,
    'r2_len': 100
    }

def set_pairing(s):
    global pairing
    if s != None:
        settings_list = s.split(',')
        for setting in settings_list:
            [lhs, rhs] = setting.split('=')
            if lhs not in pairing:
                note('invalid pairing keyword [%s]' % lhs)
            else:
                pairing[lhs] = int(rhs)
    note('set pairing mode: ' + str(pairing) + ' ("' + get_pairing_string() + '")', note_debug_level)

def get_pairing_string():
    if pairing['st_diff'] == 1:
        # different strands
        if pairing['mean'] > 0:
            return 'opp-in'
        else:
            return 'opp-out'
    else:
        # same strand
        if pairing['mean'] > 0:
            return 'col-fw'
        else:
            return 'col-bw'

def get_bowtie_pairing():
    if pairing['paired'] == 0:
        print ''
        return
    if pairing['st_diff'] == 0:
        bowtie_orientation = 'ff'
        if pairing['min'] >= 0 and pairing['max'] >= 0:
            bowtie_min_frag = max(pairing['r1_len'], pairing['min'] + pairing['r2_len'])
            bowtie_max_frag = max(pairing['r1_len'], pairing['max'] + pairing['r2_len'])
        elif pairing['min'] < 0 and pairing['max'] >= 0:
            bowtie_min_frag = max(pairing['r1_len'], pairing['r2_len'])
            bowtie_max_frag = max(pairing['r1_len'],
                                  pairing['r2_len'],
                                  pairing['r1_len'] - pairing['min'],
                                  pairing['max'] + pairing['r2_len'])
        else: # both <0
            bowtie_min_frag = max(pairing['r1_len'] - pairing['max'],
                                  pairing['r2_len'])
            bowtie_max_frag = max(pairing['r1_len'] - pairing['min'],
                                  pairing['r2_len'])
    else: #st_diff==1
        left = [pairing['min'] - (pairing['r2_len'] - 1),
                pairing['mean'] - (pairing['r2_len'] - 1),
                pairing['max'] - (pairing['r2_len'] - 1)]
        if left[1] >= 0:
            bowtie_orientation = 'fr'
        else:
            bowtie_orientation = 'rf'
        if left[0] >= 0:
            bowtie_min_frag = max(pairing['r1_len'], pairing['min'])
            bowtie_max_frag = max(pairing['r1_len'], pairing['max'])
        elif left[1] >= 0: # left[0] < 0 <= left[1] <= left[2]
            bowtie_min_frag = max(pairing['r1_len'], pairing['r2_len'])
            bowtie_max_frag = max(pairing['r1_len'], pairing['r2_len'], pairing['max'])
        else:
            bowtie_min_frag = max(pairing['r1_len'], pairing['r2_len'])
            bowtie_max_frag = max(pairing['r1_len'], pairing['r2_len'], - pairing['min'])
    print ('--%s --minins %d --maxins %d' % (bowtie_orientation, bowtie_min_frag, bowtie_max_frag))


def get_mp_st(st):
    return (st + pairing['st_diff']) % 2

def is_mp_downstream(nip, st):
    mp_st = get_mp_st(st)
    if nip == 0:
        sign_5p_difference = [+1, -1][st]
    else:
        sign_5p_difference = -1 * [+1, -1][mp_st]
    return (pairing['mean'] * sign_5p_difference > 0)

def get_mp_pos(nip, st, pos, clip_len, rlen_mp):
    mp_st = get_mp_st(st)
    pos_5p = pos[st] + [-1, +1][st] * clip_len[st]
    #note('pos_5p=' + str(pos_5p), 1)
    if nip == 0:
        sign_5p_difference = [+1, -1][st]
    else:
        sign_5p_difference = -1 * [+1, -1][mp_st]
    #note('sign_5p_difference=' + str(sign_5p_difference), 1)
    mp_pos_5p = [pos_5p + sign_5p_difference * pairing['min'],
                 pos_5p + sign_5p_difference * pairing['max']]
    mp_pos_5p.sort()
    #note('mp_pos_5p=' + str(mp_pos_5p), 1)
    sign_mp_st = [+1, -1][mp_st]
    res = [[mp_pos_5p[0], mp_pos_5p[0] + sign_mp_st * (rlen_mp - 1)],
           [mp_pos_5p[1], mp_pos_5p[1] + sign_mp_st * (rlen_mp - 1)]]
    res[0].sort()
    res[1].sort()
    return res

def pair_concordant(ms):
    [pos, clip_len] = get_pos(ms[0])
    mp_st = get_mp_st(ms[0]['st'])
    mp_pos = get_mp_pos(ms[0]['nip'], ms[0]['st'], pos, clip_len, ms[1]['len'])
    return ms[1]['st'] == mp_st and ms[1]['pos'] >= mp_pos[0][0] and ms[1]['pos'] <= mp_pos[1][0]

def check_pair_flags(p):
    if (is_first(p[0]) == is_first(p[1])
        or is_last(p[0]) == is_last(p[1])
        or is_unmapped(p[0]) != mp_unmapped(p[1])
        or is_unmapped(p[1]) != mp_unmapped(p[0])):
        return False
    elif is_unmapped(p[0]) or is_unmapped(p[1]):
        return True
    elif (is_reversed(p[0]) != mp_reversed(p[1])
          or is_reversed(p[1]) != mp_reversed(p[0])):
        return False
    return True

input_phred = 0
def autodetect_phred(s):
    min_mqv = -5
    max_mqv = 45
    common_phred_vals = [33, 64]
    for x in common_phred_vals:
        valid = True
        for c in s:
            if ord(c) - x < min_mqv or ord(c) - x > max_mqv:
                note('c=[%s] with ord=[%d] contradicts phred value [%d]' % (c, ord(c), x))
                valid = False
                break
        if valid:
            note('detected phred value [%d]' % x, note_debug_level)
            return x
    crash('could not detect input phred value')

def convert_phred(s, output_phred):
    global input_phred
    if input_phred == 0:
        input_phred = autodetect_phred(s)
    if output_phred == input_phred:
        return s
    else:
        return ''.join(map(lambda x: chr(ord(x) - input_phred + output_phred), s))

def reverse(s):
    return s[::-1]

def complement(s):
    res = map(lambda e: rc_dict[e.upper()], s)
    if isinstance(s, str):
        return ''.join(res)
    else:
        return res

def reverse_complement(s):
    return reverse(complement(s))

# given cigar ops, return possible insert sizes near left and right tails
def get_tail_insert_size(ops, min_tail_match_len=5):
    res = [0, 0]

    for [start_idx, end_idx, step, dest_idx] in [[0, len(ops), +1, 0],
                                                 [len(ops) - 1, -1, -1, 1]]:
        i = start_idx
        while i != end_idx:
            #note('considering op: ' + str(ops[i]))
            if ops[i][0] in 'HSI':
                res[dest_idx] += ops[i][1]
                i += step
            elif ops[i][0] in 'M=X':
                #note('starting match stretch')
                tmp = ops[i][1]
                j = i + step
                while j != end_idx and ops[j][0] in 'M=X':
                    note('including op in match stretch: ' + str(ops[j]), note_debug_level)
                    tmp += ops[j][1]
                    j += step
                #note('stretch len: ' + str(tmp))
                if tmp < min_tail_match_len:
                    res[dest_idx] += tmp
                    i = j
                else:
                    break
            else: # 'DN' ops
                break
    return res

#
# 0: garbage/irrelevant
# 1: unmapped (can clip&remap)
# 2: unreliable (cannot clip&remap)
# 3: reliable + big insertion at tail (can clip&remap)
# 4: reliable + solid tails
#
def classify_mapping(d, min_mqv=10, min_tail_insert_size=15, min_tail_match_len=5):
    if (is_duplicate(d)
        or failed_qc(d)
        or secondary_alignment(d)):
        return 0
    if is_unmapped(d):
        return 1
    if d['mqv'] < min_mqv:
        return 2
    if min_insert_size > 0:
        ops = parse_cigar_string(d['cigar'])
        tails = get_tail_insert_size(ops, min_tail_match_len=min_tail_match_len)
        if tails[0] >= min_tail_insert_size or tails[1] >= min_tail_insert_size:
            return 3
    return 4

def get_mapping(in_fd, header_line_hook=None):
    line = in_fd.readline()
    while len(line) != 0 and line[0] == '@':
        if header_line_hook != None:
            header_line_hook(line)
        line = in_fd.readline()
    if len(line) == 0:
        return None
    return line.strip().split('\t')

def get_mapping_set_gen(in_fd, mapping_parser=mapping_parser, header_line_hook=None,
                        ignore_duplicate=True, ignore_failed_qc=True, ignore_secondary=True,
                        check_pairing=False, check_num_mappings=False):
    map_set = []
    while True:
        tmp = get_mapping(in_fd, header_line_hook=header_line_hook)
        if tmp == None:
            if len(map_set) > 0:
                if (check_num_mappings and
                    ((pairing['paired'] == 0 and len(map_set) != 1)
                     or (pairing['paired'] == 1 and len(map_set) != 2))):
                    crash('incorrect number of mappings (check clone_ids): ' + str(map_set))
                yield map_set
            break
        d = mapping_parser(tmp)
        if 'clone_id' not in d:
            crash('mapping parser produced no clone_id: [%s]' % tmp[0])
#        if ((ignore_duplicate and is_duplicate(d))
#            or (ignore_failed_qc and failed_qc(d))
#            or (ignore_secondary and secondary_alignment(d))):
#        sys.stderr.write(str(d['clone_id'])+"\n")
        if ignore_secondary and secondary_alignment(d):
            note('ignoring mapping: ' + str(d), note_debug_level)
            continue
        if check_pairing and (d['paired'] != (pairing['paired'] == 1)):
            crash('got mapping with bad pairing: ' + str(d))
        if len(map_set) > 0 and d['clone_id'] != map_set[0]['clone_id']:
            if (check_num_mappings and
                ((pairing['paired'] == 0 and len(map_set) != 1)
                 or (pairing['paired'] == 1 and len(map_set) != 2))):
                crash('incorrect number of mappings (check clone_ids?): ' + str(map_set))
            if ((ignore_duplicate and any(map(is_duplicate, map_set)))
                or (ignore_failed_qc and any(map(failed_qc, map_set)))):
                note('ignoring mapping: ' + str(d), note_debug_level)
            else:
                yield map_set
            map_set = []
        map_set.append(d)
