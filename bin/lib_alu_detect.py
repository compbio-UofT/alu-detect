#!/usr/bin/env python

import sys
import os
import subprocess
import atexit

def crash(msg, exit_code = 1):
    print >> sys.stderr, 'error: %s' % msg
    sys.exit(exit_code)

log_level = 0

def note(msg, level = 0):
    if level == 0:
        print >> sys.stderr, 'note: %s' % msg
    elif level <= log_level:
        print >> sys.stderr, 'note [%d]: %s' % (level, msg)

def set_log_level(v):
    global log_level
    if v != log_level:
        log_level = v
        note('set log level to [%d]' % v)

def must_readline(fd, msg = 'unexpected EOF'):
    line = fd.readline()
    if not line:
        crash(msg)
    return line

gzfile_tuples = []
gzcleanup_registered = 0

def gzopen(file_name):
    global gzfile_tuples, gzcleanup_registered

    if not gzcleanup_registered:
        atexit.register(gzcleanup)
        gzcleanup_registered = 1
    if file_name == '-':
        inner_file = sys.stdin
    else:
        inner_file = open(file_name, 'r')
    inner_proc = subprocess.Popen('zc', stdin=inner_file, stdout=subprocess.PIPE, bufsize=1048576)
    outter_file = inner_proc.stdout
    gzfile_tuples.append([outter_file, inner_proc, inner_file])
    return outter_file

def gzclose(outter_file):
    global gzfile_tuples
    i = 0
    while i < len(gzfile_tuples):
        if gzfile_tuples[i][0] == outter_file:
            gzfile_tuples[i][1].kill()
            gzfile_tuples[i][2].close()
            del gzfile_tuples[i]
        else:
            i += 1
    outter_file.close()

def gzcleanup():
    global gzfile_tuples
    for e in gzfile_tuples:
        e[1].kill()

opened_fd = {}

def open_file_or_fd(file_name, mode):
    global opened_fd
    if file_name[0] == '&':
        fd = int(file_name[1:])
        if fd == sys.stdin.fileno():
            res = sys.stdin
        elif fd == sys.stdout.fileno():
            res = sys.stdout
        elif fd == sys.stderr.fileno():
            res = sys.stderr
        elif (fd, mode) in opened_fd:
            res = opened_fd[(fd, mode)]
        else:
            res = os.fdopen(fd, mode)
            opened_fd[(fd, mode)] = res
    else:
        res = open(file_name, mode)
    return res

def get_next(gen):
    try:
        result = gen.next()
    except StopIteration:
        result = None
    return result
