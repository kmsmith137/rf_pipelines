#!/usr/bin/env python

import sys
import argparse

# ArgumentParser subclass, to customize the error message
class MyParser(argparse.ArgumentParser):
    def error(self, message=None):
        print >>sys.stderr, 'Usage: rfp-run.py [-n] [-w run_name] file1.json [file2.json file3.json ...]'
        print >>sys.stderr, '    -n: runs the pipeline with no output directory'
        print >>sys.stderr, '    -w: runs the pipeline in a directory which is indexed by the web viewer (frb1 only)'

        if message is not None:
            print >>sys.stderr, '\nError:', message

        sys.exit(2)


####################################################################################################
#
# Argument parsing, checking


parser = MyParser()

parser.add_argument('json_filenames', nargs='*')
parser.add_argument('-n', action = 'store_true', help = 'runs the pipeline with no output directory')
parser.add_argument('-w', dest = 'wv_name', help = 'runs the pipeline in a directory which is indexed by the web viewer (frb1 only)')

args = parser.parse_args()

ocount = 0
if args.n:
    ocount += 1
if args.wv_name is not None:
    ocount += 1

if len(sys.argv) == 1:
    parser.error()
if ocount != 1:
    parser.error('exactly one of the flags -n,-w must be specified')
if len(args.json_filenames) == 0:
    parser.error('at least one json filename must be specified')

for j in args.json_filenames:
    if not j.endswith('.json'):
        parser.error("filename '%s' does not end in .json, currently treated as an error" % j)

if args.wv_name is not None:
    if args.wv_name.startswith('_') or args.wv_name.endswith('.json'):
        parser.error("invalid web viewer run name '%s" % args.wv_name)
    for s in args.wv_name:
        if s.isspace() or (s == '/'):
            parser.error("invalid web viewer run name '%s" % args.wv_name)


####################################################################################################
#
# Create pipeline object and run pipeline.


import json
import rf_pipelines

p = [ ]
for filename in args.json_filenames:
    f = open(filename, 'r')
    j = json.load(f)
    x = rf_pipelines.pipeline_object.from_json(j)
    p.append(x)

if len(p) > 1:
    p = rf_pipelines.pipeline(p)



if args.wv_name is not None:
    rf_pipelines.utils.run_for_web_viewer(args.wv_name, p)
else:
    assert args.nflag
    print >>sys.stderr, 'Running pipeline with no output directory'
    p.run(outdir=None)
