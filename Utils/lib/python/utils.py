from __future__ import print_function
from __future__ import with_statement

import sys
import os
import gzip
import datetime
from urllib import urlencode

from subprocess import check_output, CalledProcessError


def build_url(base_url , *res, **params):
    ''' 
    assemble a url 
    from https://stackoverflow.com/questions/15799696/library-to-build-urls-in-python
    e.g., make_url('http://example.com', 'user', 'ivan', aloholic='true', age=18)
    http://example.com/user/ivan?age=18&aloholic=true
    '''
    url = base_url
    for r in res:
        url = '{}/{}'.format(url, r)
    if params:
        url = '{}?{}'.format(url, urlencode(params))
    return url


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
    flush = False
    if kwargs:
        if 'file' in kwargs: fh = kwargs['file']
        if 'flush' in kwargs: flush = kwargs['flush']

    print('[' + str(datetime.datetime.now()) + ']\t', *objs, file=fh)
    if flush:
        fh.flush()


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

