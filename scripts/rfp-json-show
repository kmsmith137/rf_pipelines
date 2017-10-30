#!/usr/bin/env python

import sys
import argparse

# ArgumentParser subclass, to customize the error message
class MyParser(argparse.ArgumentParser):
    def error(self, message=None):
        print >>sys.stderr, 'Usage: rfp-json-show.py [-d DEPTH] file.json [x y ...]'
        print >>sys.stderr, 'The arguments "x y ..." are either string field names, or integer list indices'
        print >>sys.stderr, 'The -d flag expands output to specified depth (default 1)'

        if message is not None:
            print >>sys.stderr, '\nError:', message

        sys.exit(2)


####################################################################################################
#
# Argument parsing, checking


parser = MyParser()

if len(sys.argv) <= 1:
    parser.error()

parser.add_argument('json_filename')
parser.add_argument('modifiers', nargs='*')
parser.add_argument('-d', dest='depth', type=int, default=1, help='expands output to specified depth (default1)')

args = parser.parse_args()


####################################################################################################


import json
import rf_pipelines.utils

f = open(args.json_filename, 'r')
j = json.load(f)

for m in args.modifiers:
    if isinstance(j, dict):
        if not j.has_key(m):
            raise RuntimeError("couldn't apply modifier '%s': json object doesn't contain member with given name" % m)
        j = j[m]
    elif isinstance(j, list):
        try:
            i = int(m)
        except:
            raise RuntimeError("couldn't apply modifier '%s': expected json object, got array instead" % m)
        if (i < 0) or (i >= len(j)):
            raise RuntimeError("couldn't apply modifier '%s': index out-of-range (array length is %d)" % (m,i))
        j = j[i]
    else:
        raise RuntimeError("couldn't apply modifier '%s': expected json object/array, got %s instead" % (m,j.__class__))

rf_pipelines.utils.json_show(j, depth=args.depth)
