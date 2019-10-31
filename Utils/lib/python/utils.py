from __future__ import print_function
from __future__ import with_statement

import sys
import os
import gzip
import datetime

from subprocess import check_output, CalledProcessError

def execute_cmd(cmd, cwd=None, printCmdOnly=False, verbose=True):
    '''
    execute a command
    '''
    if verbose or printCmdOnly:
        warning("EXECUTING: ", ' '.join(cmd))
        if printCmdOnly: return
    try:
        output = check_output(cmd, cwd=cwd)
        warning(output)
    except CalledProcessError as e:
        die(e)



def gzip_file(filename, removeOriginal):
    '''
    gzip a file
    '''
    with open(filename) as f_in, gzip.open(filename + '.gz', 'wb') as f_out:
        f_out.writelines(f_in)
    if removeOriginal:
        os.remove(filename)


def truncate(s, length):
    '''
    if string s is > length, return truncated string
    with .. added to end
    '''
    return (s[:(length - 2)] + '..') if len(s) > length else s


def xstrN(s):
    '''
    wrapper for str() that handles Nones -> 'NULL'
    '''
    if s is None:
        return 'NULL'
    else:
        return str(s)

def xstr(s):
    '''
    wrapper for str() that handles Nones
    '''
    if s is None:
        return ""
    else:
        return str(s)


def warning(*objs, **kwargs):
    '''
    print messages to stderr
    '''
    fh = sys.stderr
    if kwargs:
        if 'file' in kwargs: fh = kwargs['file']

    print('[' + str(datetime.datetime.now()) + ']\t', *objs, file=fh)


def qw(s, returnTuple=False):
    '''
    mimics perl's qw function
    usage: qw('a b c') will yield ['a','b','c']
    returnTuple: return a tuple if true, otherwise return list
    '''
    if returnTuple:
        return tuple(s.split())
    else:
        return s.split()


def create_dir(dirName):
    '''
    check if directory exists in the path, if not create
    '''
    try:
        os.stat(dirName)
    except OSError:
        os.mkdir(dirName)

    return dirName


def verify_path(fileName, isDir=False):
    '''
    verify that a file exists
    if isDir is True, just verifies that the path
    exists
    '''
    if isDir:
        return os.path.exists(fileName)
    else:
        return os.path.isfile(fileName)
        

def die(message):
    '''
    mimics Perl's die function
    '''
    warning(message)
    sys.exit(1)

