import os
import sys
#import helpers


###########################################################################
# helper routines
###########################################################################

#==========================================================================
# open the file (default "std.out") and parse for `findstr`, 
# in the last `n_tail` lines, 
# returns if found and line of the file
#==========================================================================
def find_last_occurence(findstr,n_tail=40,filepath="std.out",**kwargs):
   with open(filepath,'r') as fp: 
      lines=fp.readlines()
   for line in reversed(lines[-n_tail:]) :
      if findstr in line : 
         return True,line
   return False, ' '  


###########################################################################
#  ROUTINES FOR `test_opts`: 
#    * MUST have only keyword arguments (not positional!). 
#    * MUST always have the container argument `**kwargs` as last argument!
#    * MUST return status (True/False) and message
#    * ASSUMES current directory is the one where the code was executed!
###########################################################################

#==========================================================================
# find the line of the string `whichError` and check all entries after :  to be < err_tol
#==========================================================================
def check_error(whichError='L_inf ', err_tol=1.0e-12,**kwargs ):
   [found,line]=find_last_occurence(whichError,**kwargs)
   assert found, ('check_error: did not find "%s" in std.out' % (whichError))
   errors=[float(x) for x in (line.split(":")[1]).split()]
   check=all([e < err_tol for e in errors])
   msg=('===> check_error: "%s" < %e = %s' % (line.strip(),err_tol,check))
   return check,msg

#==========================================================================
# find occurence of the final_str
#==========================================================================
def run_finalized(final_str='FLUXO FINISHED!',**kwargs):
   [found,line]=find_last_occurence(final_str,n_tail=5,**kwargs)
   return found, ('%s found' % final_str)


    
###########################################################################
# DEFINITION OF THE DICTIONARY FOR THE JOBS
# Description:
# all jobs are stored as a dictionary, each job has a UNIQUE name and is a sub-dictionary
# with mandatory keys
# * 'case': unique case index,
# * 'tags': [list of values], 
# * build_opts : {sub-dictionary of key: value pairs, key is cmake variable which set to the value}
#                would be appended to global build opts and  overwrite existing entries
# and an OPTIONAL KEY
# * 'run_opts':  set-up runs of the build, 
#    - name of the run refers to the path to the test folder, which must be unique.
#      the test folder MUST contain a `parameter.ini` file!
#    - tags can be used
#   and an OPTIONAL KEY 
#    - 'test_opts': check the result of the run, again a sub-dictionary, each test must have a unique name
#      sub-dictionary must have
#       - 'func' :  python function name, without positional arguments, only key=value arguments!
#      and if the function needs additional parameters (like tolerances),
#       - 'f_kwargs': keyword arguments of the function 'func'
#
# DEFINE NEW JOBS HERE!!!
#####################################################################################
def job_definition():
   jobs ={}
   #============================================================================
   caseID=400
   baseopts={ 'FLUXO_EQNSYSNAME' :'navierstokes',
              'FLUXO_TESTCASE'   :'default'
            }
   
   
   #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
   caseID=caseID+1
   jobs['navierstokes_type1_br1_GL']={
          'case': caseID ,
          'tags': ['navierstokes','standardDG','br1','GL'] ,
          'build_opts':{**baseopts,
                        'FLUXO_DISCTYPE'         :'1',
                        'FLUXO_DISC_NODETYPE'    :'GAUSS-LOBATTO',
                        'FLUXO_PARABOLIC'        :'ON',
                        'FLUXO_PARABOLIC_LIFTING':'br1',
                        'FLUXO_PARABOLIC_LIFTING_VAR':'cons_var',
                       },
          'run_opts': {'runs/navst/freestream/conforming':
                         {'tags': ['freestream','conforming'] ,
                          'test_opts':{'finished':{'func': run_finalized}  ,
                                       'err_Linf':{'func': check_error ,'f_kwargs': {'whichError':'L_inf ','err_tol': 1e-11} } ,
                                       'err_L2'  :{'func': check_error ,'f_kwargs': {'whichError':'L_2 '  ,'err_tol': 1e-11} } ,
                                      },
                         },
                       'runs/navst/freestream/nonconforming':
                         {'tags': ['freestream','nonconforming'],
                          'test_opts':{'finished':{  'func': run_finalized}  ,
                                       'err_Linf':{  'func': check_error ,'f_kwargs': {'whichError':'L_inf ','err_tol': 1e-11} } ,
                                       'err_L2'  :{  'func': check_error ,'f_kwargs': {'whichError':'L_2 '  ,'err_tol': 1e-11} },
                                      },
                         },
                      },
         }
   #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
   caseID=caseID+1
   jobs['navierstokes_type2_br2_entropy_var']={
          'case': caseID,
          'tags': [ 'navierstokes','split-form','GL','br2','entropy_var'],
          'build_opts':{**baseopts,
                        'FLUXO_DISCTYPE'         :'2',
                        'FLUXO_DISC_NODETYPE'    :'GAUSS-LOBATTO',
                        'FLUXO_PARABOLIC'        :'ON',
                        'FLUXO_PARABOLIC_LIFTING':'br2',
                        'FLUXO_PARABOLIC_LIFTING_VAR':'entropy_var',
                        'FLUXO_EQN_VOLFLUX'      :'5',
                       },
          'run_opts': {'runs/navst/freestream/conforming': 
                         {'tags': ['freestream'],
                          'test_opts':{'finished':{'func': run_finalized },
                                       'err_Linf':{'func': check_error ,'f_kwargs': {'whichError':'L_inf ','err_tol': 1e-11} },
                                      },
                         },
                       'runs/navst/freestream/nonconforming': 
                         {'tags': ['freestream'],
                         },
                      },
         }
   #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

   return jobs

###test
#jobs=job_definition()
#ff= jobs['navierstokes_type1_br1_GL']['run_opts']['runs/navst/freestream/conforming']['test_opts']['err_Linf']
#
#[stat,msg]=ff['func'](**ff['f_kwargs'])
#print(stat,msg)
