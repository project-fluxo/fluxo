import os
import sys

#from helpers import modify_prm  # include via file link to ../../../tests/helpers.py

jobs ={}
#SUCCESSFUL!
jobs['master_mhd_para']  ={'exec' :'../../build_master_mhd_para/bin/fluxo',
                   'param':'parameter_mhd.ini'}
#jobs['mortarmaster_ns_para']  ={'exec' :'../../build_mortarmaster_ns_para/bin/fluxo',
#                   'param':'parameter_ns.ini'}
#jobs['linadv_nofix']  ={'exec' :'../../build_linadv_nofix/bin/fluxo',
#                   'param':'parameter_linadv.ini'}
#jobs['linadv_fix']  ={'exec' :'../../build_linadv_fix/bin/fluxo',
#                   'param':'parameter_linadv.ini'}
#jobs['ns_nopara']={'exec' :'../../build_master_ns_nopara/bin/fluxo',
#                   'param':'parameter_ns.ini'}
#jobs['ns_conf_para']  ={'exec' :'../../build_master_ns_para/bin/fluxo',
#                   'param':'parameter_ns_conf.ini'}
#jobs['ns_conf_paracons']  ={'exec' :'../../build_master_ns_paracons/bin/fluxo',
#                   'param':'parameter_ns_conf.ini'}
#jobs['ns_para_prim_fix']  ={'exec' :'../../build_ns_para_prim_fix/bin/fluxo',
#                   'param':'parameter_ns.ini'}
#jobs['ns_para_br2cons_fix']  ={'exec' :'../../build_master_ns_para_br2cons/bin/fluxo',
#                   'param':'parameter_ns.ini'}
#jobs['ns_para_fix']  ={'exec' :'../../build_master_ns_para/bin/fluxo',
#                   'param':'parameter_ns.ini'}
#unsuccessfull...
#jobs['jesse_ns_para']  ={'exec' :'../../build_jesse_ns_para/bin/fluxo',
#                   'param':'parameter_ns.ini'}
#jobs['mortarmaster_off_ns_para']  ={'exec' :'../../build_mortarmaster_off_ns_para/bin/fluxo',
#                   'param':'parameter_ns.ini'}
#jobs['ns_paracons']  ={'exec' :'../../build_master_ns_paracons/bin/fluxo',
#                   'param':'parameter_ns.ini'}
#jobs['ns_para_br2cons']  ={'exec' :'../../build_master_ns_para_br2cons/bin/fluxo',
#                           'param':'parameter_ns.ini'}

procs=["1","2","3"]


       
for j_name,job in jobs.items():
  for proc in procs:
      print('running job %s on procs %s ' %(j_name,proc))
      logfile='log_'+j_name+'_np'+proc
      os.system('mpirun -np '+proc+' '+ job['exec'] +' '+ job['param'] + ' > '+logfile)

    
