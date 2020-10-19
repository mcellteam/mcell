#!/usr/bin/env python
import sys
import os
import math


if __name__ == '__main__':

  if (len(sys.argv)<5):
    print ('\nUsage: %s seed_list_spec seeds_per_node rules_file mdl_file\n' % (sys.argv[0]))
    print ('\n    Example: %s "1:100:2" 8 Scene.mdlr_rules.xml Scene.main.mdl\n' % (sys.argv[0]))
    print ('        This would submit mcell jobs to simulate Scene.main.mdl')
    print ('         using seeds 1 through 100 in steps of 2,')
    print ('         with 8 individual simulations running concurrently in each job \n\n')
    sys.exit()
  else:
    seed_list_spec = sys.argv[1]
    seeds_per_node = int(sys.argv[2])
    rules_file = sys.argv[3]
    mdl_file = sys.argv[4]
   
    seed_list = []
    seed_list_toks = seed_list_spec.split(',') 
    for seed_list_tok in seed_list_toks:
     
      seed_range_toks = seed_list_tok.split(':') 
      if (len(seed_range_toks) > 1):
        if (len(seed_range_toks) == 2):
          seed_list.extend(range(int(seed_range_toks[0]),int(seed_range_toks[1])+1))
        else:
          seed_list.extend(range(int(seed_range_toks[0]),int(seed_range_toks[1])+1,int(seed_range_toks[2])))

      else:
        seed_list.append(int(seed_list_tok))


    mcell_exec = '$HOME/bin/mcell3.5/mcell3r.py'
    mcell_cmd = 'source ~/.bashrc ; %s -s $seed -r %s -m %s' % (mcell_exec,rules_file,mdl_file)


    job_script_template = '%s; for seed in ${seeds[@]}; do { %s & }; done; wait;'

    job_name_template = 'mcell_%04d_%04d'

    qsub_option_template = '-l nodes=1 -l walltime=48:00:00 -m n -A kennedy-lab -d ./ -S /bin/bash -N %s -o qsub_%s.out -e qsub_%s.err'

  for start_i in range(0,len(seed_list),seeds_per_node):
    end_i = start_i + seeds_per_node - 1
    if end_i > len(seed_list)-1:
      end_i = len(seed_list)-1
   
    seed_chunk = 'seeds=('
    for seed in seed_list[start_i:end_i+1]:
      seed_chunk = seed_chunk + ' ' + str(seed)
    seed_chunk = seed_chunk + ' )'
    
    job_script = job_script_template % (seed_chunk, mcell_cmd)
    job_name = job_name_template % (seed_list[start_i], seed_list[end_i])
    qsub_options = qsub_option_template % (job_name, job_name, job_name)

    cmd = "echo '%s' | qsub %s -" % (job_script, qsub_options)
    print(cmd)
    os.system(cmd)

