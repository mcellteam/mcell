#!/usr/bin/env python

import argparse
import os
import sys

def define_console():
    parser = argparse.ArgumentParser(description='MCellR runner/postprocessor')
    parser.add_argument('-v', '--version',            action='store_true',            help='print program version and exit')
    parser.add_argument('-f', '--fullversion',            action='store_true', help='print detailed version report and exit')
    parser.add_argument('-s', '--seed',           type=str,            help='choose random sequence number (default: 1)')
    parser.add_argument('-i', '--iterations',   type=str,            help='override iterations in MDL file')
    parser.add_argument('-l', '--logfile', type=str,            help='send output log to file (default: stdout)')
    parser.add_argument('-x', '--logfreq', type=str,        help='output log frequency')
    parser.add_argument('-e', '--errfile', type=str,            help='send errors to log file (default: stderr)')
    parser.add_argument('-p', '--checkpoint_infile', type=str,            help='read checkpoint file')
    parser.add_argument('-o', '--checkpoint_outfile', type=str,            help='write checkpoint file')
    parser.add_argument('-q', '--quiet', action='store_true',            help='suppress all unrequested output except for errors')
    parser.add_argument('-w', '--with_checks', type=str,            help='yes/no : default yes   perform check of the geometry for coincident walls')
    parser.add_argument('-r', '--rules',            type=str,            help='rules file',              required=True)
    parser.add_argument('-m', '--mdl_infile',            type=str,            help='MDL main input file',              required=True)
    return parser

if __name__ == "__main__":
    parser = define_console()
    args = parser.parse_args()
  
    print(args)

    cmd_args = ''

    if args.version:
       cmd_args = cmd_args + ' -version'
    if args.fullversion:
       cmd_args = cmd_args + ' -fullversion'
    if args.seed:
       cmd_args = cmd_args + ' -seed %s' % (args.seed)
    else:
       cmd_args = cmd_args + ' -seed 1'
       args.seed = '1'
    if args.iterations:
       cmd_args = cmd_args + ' -iterations %s' % (args.iterations)
    if args.logfile:
       cmd_args = cmd_args + ' -logfile %s' % (args.logfile)
    if args.logfreq:
       cmd_args = cmd_args + ' -logfreq %s' % (args.logfreq)
    if args.errfile:
       cmd_args = cmd_args + ' -errfile %s' % (args.errfile)
    if args.checkpoint_infile:
       cmd_args = cmd_args + ' -checkpoint_infile %s' % (args.checkpoint_infile)
    if args.checkpoint_outfile:
       cmd_args = cmd_args + ' -checkpoint_outfile %s' % (args.checkpoint_outfile)
    if args.quiet:
       cmd_args = cmd_args + ' -quiet'
    if args.with_checks:
       cmd_args = cmd_args + ' -with_checks %s' % (args.with_checks)
    if args.rules:
       cmd_args = cmd_args + ' -rules %s' % (args.rules)
    if args.mdl_infile:
       cmd_args = cmd_args + ' %s' % (args.mdl_infile)

    print(cmd_args)

    script_path = os.path.dirname(os.path.realpath(__file__))

    mcell_cmd = os.path.join(script_path, 'mcell')

    mcell_cmd = mcell_cmd + ' ' + cmd_args

    postproc_cmd = os.path.join(script_path, 'postprocess_mcell3r.py')

    postproc_cmd = postproc_cmd + ' ' + args.seed + ' ' + args.rules

    os.system(mcell_cmd)
    os.system(postproc_cmd)

