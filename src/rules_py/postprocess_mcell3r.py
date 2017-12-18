#!/usr/bin/env python

import os
import sys
import shutil

def postprocess(mcellr_react_filename, run_seed):
  print('')
  print ( '  Postprocessing: %s' % (mcellr_react_filename) )

  mcellr_react_dir = './'
  react_dir = os.path.join(mcellr_react_dir, 'react_data')

#  if os.path.exists(react_dir):
#      shutil.rmtree(react_dir,ignore_errors=True)
#  if not os.path.exists(react_dir):
#      os.makedirs(react_dir)

  seed_dir = 'seed_%05d' % run_seed

  react_seed_dir = os.path.join(react_dir, seed_dir)

#  if os.path.exists(react_seed_dir):
#      shutil.rmtree(react_seed_dir,ignore_errors=True)
  if not os.path.exists(react_seed_dir):
      os.makedirs(react_seed_dir)

  # Read the MCellR data file and split into a list of rows where each row is a list of columns
  mcellr_react_file = open (mcellr_react_filename,'r')
  all_react_data = mcellr_react_file.read()
  react_data_all = [ [t.strip() for t in s.split(',') if len(t.strip()) > 0] for s in all_react_data.split('\n') if len(s) > 0 ]
  react_data_header = []
  react_data_rows = []
  if len(react_data_all) > 0:
    react_data_header = react_data_all[0]
    react_data_rows = react_data_all[1:]

  for col in range(1,len(react_data_header)):
    out_file_name = os.path.join ( react_seed_dir, react_data_header[col] + '.dat' )
    print ( '    Writing data to ' + out_file_name )
    f = open(out_file_name,'w')
    for row in react_data_rows:
#      print ( '  ' + row[0] + ' ' + row[col] )
      f.write ( row[0] + ' ' + row[col] + '\n' )
    f.close()


if __name__ == '__main__':

  if (len(sys.argv) < 3):
    print('\nUsage: %s run_seed mcellr_rules_file\n' % (sys.argv[0]))
    print('         Split contents of MCellR gdat reaction output into')
    print('         separate output files in react_data/seed_nnnnn/ subdirs\n')
    exit(1)

  run_seed = int(sys.argv[1])
  rules_file = sys.argv[2]

  mols_data_file = rules_file + '.seed_%05d.gdat' % (run_seed)
  reactions_data_file = rules_file + '_reactions.seed_%05d.gdat' % (run_seed)

  print('')
  print ( 'Postprocessing MCellR Reaction Output...')
  postprocess(mols_data_file, run_seed)
  postprocess(reactions_data_file, run_seed)

  print('')
  print ( 'Done Postprocessing MCellR Reaction Output' )
  print('')

