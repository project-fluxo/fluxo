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
def find_last_occurence(findstr,n_tail=40,stdout_filepath="std.out",**kwargs):
   with open(stdout_filepath,'r') as fp: 
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
   msg=('check_error: "%s" < %e = %s' % (line.strip(),err_tol,check))
   return check,msg

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
# NOTE: only python >3.7 guarantees that the order of the dicitionary corresponds to the initialization order
#
# DEFINE NEW JOBS HERE!!! 
# => THEN CHECK THE CORRECT DEFINITION OF THE JOBS, running job_init from job_check.py with 
#    ` python -c 'from job_check import * ; job_init(list_all=True)' `
# => job_init is also used to run the build tests.
#####################################################################################
def job_definition():
   jobs ={}
   #============================================================================
   #============================================================================
   #first group, LINADV, 100 < caseID <200
   #============================================================================
   #============================================================================
   caseID=100 
   baseopts={"FLUXO_EQNSYSNAME"       :"linearscalaradvection",
             "_BUILD_FLUXO_POST"      :"ON"}
   run_opt_fsp={'runs/linadv/freestream': 
                 {'tags': ['linadv','freestream'],
                  'test_opts':{ 'err_Linf':{'func': check_error ,
                                            'f_kwargs': {'whichError':'L_inf ','err_tol': 1e-12} },
                              },
                 },
             }
   #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
   caseID=caseID+1
   jobs['linadv_type1_GL_br1']={
          'case': caseID,
          'tags': ['linadv','standardDG','GL','br1'],
          'build_opts':{**baseopts,
                        'FLUXO_DISCTYPE'         :'1',
                        'FLUXO_DISC_NODETYPE'    :'GAUSS-LOBATTO',
                        'FLUXO_PARABOLIC'        :'ON',
                        'FLUXO_PARABOLIC_LIFTING':'br1',
                       },
          'run_opts': {**run_opt_fsp,}
         }
      
      
   #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
   caseID=caseID+1
   jobs['linadv_type1_gauss_br1']={
          'case': caseID,
          'tags': ['linadv','standardDG','gauss','br1'],
          'build_opts':{**baseopts,
                        'FLUXO_DISCTYPE'         :'1',
                        'FLUXO_DISC_NODETYPE'    :'GAUSS',
                        'FLUXO_PARABOLIC'        :'ON',
                        'FLUXO_PARABOLIC_LIFTING':'br1'
                       },
          'run_opts': {**run_opt_fsp,}
         }
   #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
   caseID=caseID+1
   jobs['linadv_split_br1']={
          'case': caseID,
          'tags': ['linadv','split-form','GL','br1'],
          'build_opts':{**baseopts,
                        'FLUXO_DISCTYPE'         :'2',
                        'FLUXO_DISC_NODETYPE'    :'GAUSS-LOBATTO',
                        'FLUXO_PARABOLIC'        :'ON',
                        'FLUXO_PARABOLIC_LIFTING':'br1',
                       },
          'run_opts': {**run_opt_fsp,}
         }
   #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
   caseID=caseID+1
   jobs['linadv_nopara_type1_gauss']={
          'case': caseID,
          'tags': ['linadv','standardDG','GL'],
          'build_opts':{**baseopts,
                        'FLUXO_DISCTYPE'         :'1',
                        'FLUXO_DISC_NODETYPE'    :'GAUSS',
                        'FLUXO_PARABOLIC'        :'OFF',
                       },
          'run_opts': {**run_opt_fsp,}
         }
   #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
   caseID=caseID+1
   jobs['linadv_type1_GL_br2']={
          'case': caseID,
          'tags': ['linadv','standardDG','GL','br2'],
          'build_opts':{**baseopts,
                        'FLUXO_DISCTYPE'         :'1',
                        'FLUXO_DISC_NODETYPE'    :'GAUSS-LOBATTO',
                        'FLUXO_PARABOLIC'        :'ON',
                        'FLUXO_PARABOLIC_LIFTING':'br2',
                       },
          'run_opts': {**run_opt_fsp,}
         }
   #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
   caseID=caseID+1
   jobs['linadv_type1_gauss_br2']={
          'case': caseID,
          'tags': ['linadv','standardDG','gauss','br2'],
          'build_opts':{**baseopts,
                        'FLUXO_DISCTYPE'         :'1',
                        'FLUXO_DISC_NODETYPE'    :'GAUSS',
                        'FLUXO_PARABOLIC'        :'ON',
                        'FLUXO_PARABOLIC_LIFTING':'br2',
                       },
          'run_opts': {**run_opt_fsp,}
         }
   #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
   volfluxes=['-1','0','1']
   for  vvv  in range(0,len(volfluxes)):
      volflux=volfluxes[vvv]
      caseID=caseID+1
      jobs['linadv_split_nopara_volflux_'+volflux]={
          'case': caseID,
          'tags': ['linadv','split-form','GL'],
          'build_opts':{**baseopts,
                        'FLUXO_DISCTYPE'         :'2',
                        'FLUXO_DISC_NODETYPE'    :'GAUSS-LOBATTO',
                        'FLUXO_PARABOLIC'        :'OFF',
                        'FLUXO_EQN_VOLFLUX'      :volflux,
                       },
          'run_opts': {**run_opt_fsp,}
         }
   #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
   #============================================================================
   #============================================================================
   #next group, MAXWELL, 200 < caseID <300
   #============================================================================
   #============================================================================
   caseID=200 
   baseopts={ 'FLUXO_EQNSYSNAME'       :'maxwell',
              'FLUXO_PARABOLIC'        :'OFF',
            }
   run_opt_fsp={'runs/maxwell/freestream': 
                {'tags': ['maxwell','freestream'],
                 'test_opts':{ 'err_Linf':{'func': check_error ,
                                           'f_kwargs': {'whichError':'L_inf ','err_tol': 1e-12} },
                             },
                },
              }
   #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
   caseID=caseID+1
   jobs['maxwell_type1_GL']={
          'case': caseID,
          'tags': ['maxwell','standardDG','GL'],
          'build_opts':{**baseopts,
                        'FLUXO_DISCTYPE'         :'1',
                        'FLUXO_DISC_NODETYPE'    :'GAUSS-LOBATTO',
                       },
          'run_opts': {**run_opt_fsp,}
         }
   #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
   caseID=caseID+1
   jobs['maxwell_split']={
          'case': caseID,
          'tags': ['maxwell','split-form','GL'],
          'build_opts':{**baseopts,
                        'FLUXO_DISCTYPE'         :'2',
                        'FLUXO_DISC_NODETYPE'    :'GAUSS-LOBATTO',
                       },
          'run_opts': {**run_opt_fsp,}
         }
   #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
   caseID=caseID+1
   jobs['maxwell_type1_gauss']={
          'case': caseID,
          'tags': ['maxwell','standardDG','gauss'],
          'build_opts':{**baseopts,
                        'FLUXO_DISCTYPE'         :'1',
                        'FLUXO_DISC_NODETYPE'    :'GAUSS',
                       },
          'run_opts': {**run_opt_fsp,}
         }
   #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
   volfluxes=['-1','0','1']
   for  vvv  in range(0,len(volfluxes)):
      volflux=volfluxes[vvv]
      jobs['maxwell_split_volflux_'+volflux]={
          'case': caseID,
          'tags': ['maxwell','split-form','GL'],
          'build_opts':{**baseopts,
                        'FLUXO_DISCTYPE'         :'2',
                        'FLUXO_DISC_NODETYPE'    :'GAUSS-LOBATTO',
                        'FLUXO_EQN_VOLFLUX'      :volflux,
                       },
          'run_opts': {**run_opt_fsp,}
         }
   #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
   #============================================================================
   #============================================================================
   #third group, MHD, 300 < caseID <400
   #============================================================================
   #============================================================================
   caseID=300 
   baseopts={ 'FLUXO_EQNSYSNAME' :'mhd',
            }
   
   run_opt_fsp_conf={'runs/mhd/freestream/conforming':
         {'tags': ['mhd','freestream','conforming'] ,
          'test_opts':{'err_Linf':{'func': check_error ,
                                   'f_kwargs': {'whichError':'L_inf ','err_tol': 1e-11} } ,
                      },
         },
      }
   run_opt_fsp_nonconf={'runs/mhd/freestream/nonconforming':
         {'tags': ['mhd','freestream','nonconforming'] ,
          'test_opts':{'err_Linf':{'func': check_error ,
                                   'f_kwargs': {'whichError':'L_inf ','err_tol': 1e-11} } ,
                      },
         },
      }
   #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
   caseID=caseID+1
   jobs['mhd_type1_br1_cons_GL']={
          'case': caseID ,
          'tags': ['mhd','standardDG','br1','GL','GLM'] ,
          'build_opts':{**baseopts,
                        'FLUXO_DISCTYPE'         :'1',
                        'FLUXO_DISC_NODETYPE'    :'GAUSS-LOBATTO',
                        'FLUXO_EQN_GLM'          :'ON',
                        'FLUXO_PARABOLIC'        :'ON',
                        'FLUXO_PARABOLIC_LIFTING':'br1',
                        'FLUXO_PARABOLIC_LIFTING_VAR':'cons_var',
                       },
          'run_opts': {**run_opt_fsp_conf, 
                       **run_opt_fsp_nonconf,
                      }
         }
   #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
   caseID=caseID+1
   jobs['mhd_split_noglm_noncons_br1entr_ecvolflux']={
          'case': caseID ,
          'tags': ['mhd','split-form','br1','GL','NONCONS','ECflux'] ,
          'build_opts':{**baseopts,
                        'FLUXO_DISCTYPE'         :'2',
                        'FLUXO_DISC_NODETYPE'    :'GAUSS-LOBATTO',
                        'FLUXO_EQN_GLM'          :'ON',
                        'FLUXO_EQN_NONCONS'      :'ON',
                        'FLUXO_EQN_NONCONS_GLM'  :'ON',
                        'FLUXO_PARABOLIC'        :'ON',
                        'FLUXO_PARABOLIC_LIFTING':'br1',
                        'FLUXO_PARABOLIC_LIFTING_VAR':'entropy_var',
                        'FLUXO_EQN_VOLFLUX'      :'10',
                       },
          'run_opts': {**run_opt_fsp_conf, 
                       **run_opt_fsp_nonconf,
                      }
         }
   #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
   caseID=caseID+1
   jobs['mhd_split_glm_noncons_br1entr_ecvolflux_testcase_mhdeq']={
          'case': caseID ,
          'tags': ['mhd','split-form','br1','GL','GLM','NONCONS','ECflux'] ,
          'build_opts':{**baseopts,
                        'FLUXO_DISCTYPE'         :'2',
                        'FLUXO_DISC_NODETYPE'    :'GAUSS-LOBATTO',
                        'FLUXO_EQN_GLM'          :'ON',
                        'FLUXO_EQN_NONCONS'      :'ON',
                        'FLUXO_EQN_NONCONS_GLM'  :'ON',
                        'FLUXO_PARABOLIC'        :'ON',
                        'FLUXO_PARABOLIC_LIFTING':'br1',
                        'FLUXO_PARABOLIC_LIFTING_VAR':'entropy_var',
                        'FLUXO_EQN_VOLFLUX'      :'12',
                        'FLUXO_TESTCASE'         :'mhd_equilibrium'
                       },
          'run_opts': {**run_opt_fsp_conf, 
                       **run_opt_fsp_nonconf,
                      }
         }
   #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
   caseID=caseID+1
   jobs['mhd_split_noglm_noncons_nopara_testcase_tgv']={
          'case': caseID ,
          'tags': ['mhd','split-form','GL','NONCONS','ECflux'] ,
          'build_opts':{**baseopts,
                        'FLUXO_DISCTYPE'         :'2',
                        'FLUXO_DISC_NODETYPE'    :'GAUSS-LOBATTO',
                        'FLUXO_EQN_GLM'          :'OFF',
                        'FLUXO_EQN_NONCONS'      :'ON',
                        'FLUXO_PARABOLIC'        :'OFF',
                        'FLUXO_EQN_VOLFLUX'      :'10',
                        'FLUXO_TESTCASE'         :'taylorgreenvortex',
                       },
          'run_opts': {**run_opt_fsp_conf, 
                       **run_opt_fsp_nonconf,
                      }
         }
   #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
   caseID=caseID+1
   jobs['mhd_type1_br2prim_gauss']={
          'case': caseID ,
          'tags': ['mhd','standardDG','br2','gauss','GLM'] ,
          'build_opts':{**baseopts,
                        'FLUXO_DISCTYPE'         :'1',
                        'FLUXO_DISC_NODETYPE'    :'GAUSS',
                        'FLUXO_EQN_GLM'          :'ON',
                        'FLUXO_PARABOLIC'        :'ON',
                        'FLUXO_PARABOLIC_LIFTING':'br2',
                        'FLUXO_PARABOLIC_LIFTING_VAR':'prim_var',
                       },
          'run_opts': {**run_opt_fsp_conf, 
                       **run_opt_fsp_nonconf,
                      }
         }
   #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
   caseID=caseID+1
   jobs['mhd_split_glm_no-noncons_br2cons']={
          'case': caseID ,
          'tags': ['mhd','split-form','br2','GL','GLM'] ,
          'build_opts':{**baseopts,
                        'FLUXO_DISCTYPE'         :'2',
                        'FLUXO_DISC_NODETYPE'    :'GAUSS-LOBATTO',
                        'FLUXO_EQN_GLM'          :'ON',
                        'FLUXO_EQN_NONCONS'      :'OFF',
                        'FLUXO_PARABOLIC'        :'ON',
                        'FLUXO_PARABOLIC_LIFTING':'br2',
                        'FLUXO_PARABOLIC_LIFTING_VAR':'cons_var',
                       },
          'run_opts': {**run_opt_fsp_conf, 
                       **run_opt_fsp_nonconf,
                      }
         }
   #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
   caseID=caseID+1
   jobs['mhd_split_glm_nonconsglmbrack_br1prim']={
          'case': caseID ,
          'tags': ['mhd','split-form','br1','GL','GLM','NONCONS'] ,
          'build_opts':{**baseopts,
                        'FLUXO_DISCTYPE'         :'2',
                        'FLUXO_DISC_NODETYPE'    :'GAUSS-LOBATTO',
                        'FLUXO_EQN_GLM'          :'ON',
                        'FLUXO_EQN_NONCONS'      :'ON',
                        'FLUXO_EQN_NONCONS_TYPE' :'Brackbill',
                        'FLUXO_EQN_NONCONS_GLM'  :'ON',
                        'FLUXO_PARABOLIC'        :'ON',
                        'FLUXO_PARABOLIC_LIFTING':'br1',
                        'FLUXO_PARABOLIC_LIFTING_VAR':'prim_var',
                       },
          'run_opts': {**run_opt_fsp_conf, 
                       **run_opt_fsp_nonconf,
                      }
         }
   #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
   caseID=caseID+1
   jobs['mhd_type1_GL_br1entr']={
          'case': caseID ,
          'tags': ['mhd','standardDG','br1','GL','GLM'] ,
          'build_opts':{**baseopts,
                        'FLUXO_DISCTYPE'         :'1',
                        'FLUXO_DISC_NODETYPE'    :'GAUSS-LOBATTO',
                        'FLUXO_EQN_GLM'          :'ON',
                        'FLUXO_PARABOLIC'        :'ON',
                        'FLUXO_PARABOLIC_LIFTING':'br1',
                        'FLUXO_PARABOLIC_LIFTING_VAR':'entropy_var',
                       },
          'run_opts': {**run_opt_fsp_conf, 
                       **run_opt_fsp_nonconf,
                      }
         }
   #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
   caseID=caseID+1
   jobs['mhd_type1_GL_br2prim_TCmhd']={
          'case': caseID ,
          'tags': ['mhd','standardDG','br2','GL','GLM'] ,
          'build_opts':{**baseopts,
                        'FLUXO_DISCTYPE'         :'1',
                        'FLUXO_DISC_NODETYPE'    :'GAUSS-LOBATTO',
                        'FLUXO_EQN_GLM'          :'ON',
                        'FLUXO_PARABOLIC'        :'ON',
                        'FLUXO_PARABOLIC_LIFTING':'br2',
                        'FLUXO_PARABOLIC_LIFTING_VAR':'prim_var',
                        'FLUXO_TESTCASE'         :'mhd_equilibrium'
                       },
          'run_opts': {**run_opt_fsp_conf, 
                       **run_opt_fsp_nonconf,
                      }
         }
   #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
   volfluxes=['0','10','12']
   for  vvv  in range(0,len(volfluxes)):
      volflux=volfluxes[vvv]
      caseID=caseID+1
      jobs['mhd_split_noglm_noncons_nopara_volflux_'+volflux]={
          'case': caseID ,
          'tags': ['mhd','split-form','GL','NONCONS'] ,
          'build_opts':{**baseopts,
                        'FLUXO_DISCTYPE'         :'2',
                        'FLUXO_DISC_NODETYPE'    :'GAUSS-LOBATTO',
                        'FLUXO_EQN_GLM'          :'OFF',
                        'FLUXO_EQN_NONCONS'      :'ON',
                        'FLUXO_PARABOLIC'        :'OFF',
                        'FLUXO_EQN_VOLFLUX'      : volflux,
                       },
          'run_opts': {**run_opt_fsp_conf, 
                       **run_opt_fsp_nonconf,
                      }
         }
   #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
   caseID=caseID+1
   jobs['mhd_split_glm_no-noncons_nopara']={
          'case': caseID ,
          'tags': ['mhd','split-form','GL','GLM'] ,
          'build_opts':{**baseopts,
                        'FLUXO_DISCTYPE'         :'2',
                        'FLUXO_DISC_NODETYPE'    :'GAUSS-LOBATTO',
                        'FLUXO_EQN_GLM'          :'ON',
                        'FLUXO_EQN_NONCONS'      :'OFF',
                        'FLUXO_PARABOLIC'        :'OFF',
                       },
          'run_opts': {**run_opt_fsp_conf, 
                       **run_opt_fsp_nonconf,
                      }
         }
   #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
   run_opt_fsp_SC_first={'runs/mhd/freestream/SC_firstOrder':
         {'tags': ['mhd','freestream','conforming','SC','firstorder'] ,
          'test_opts':{'err_Linf':{'func': check_error ,
                                   'f_kwargs': {'whichError':'L_inf ','err_tol': 1e-11} } ,
                      },
         },
      }
   #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
   caseID=caseID+1
   jobs['mhd_split_glm_noncons_nopara_SC']={
          'case': caseID ,
          'tags': ['mhd','split-form','SC','GL','GLM','NONCONS'] ,
          'build_opts':{**baseopts,
                        'FLUXO_DISCTYPE'         :'2',
                        'FLUXO_DISC_NODETYPE'    :'GAUSS-LOBATTO',
                        'FLUXO_EQN_GLM'          :'ON',
                        'FLUXO_EQN_NONCONS'      :'ON',
                        'FLUXO_EQN_NONCONS_GLM'  :'ON',
                        'FLUXO_PARABOLIC'        :'OFF',
                        'FLUXO_SHOCKCAPTURE'     :'ON',
                        'FLUXO_SHOCKCAP_NFVSE'   :'ON',
                        'FLUXO_SHOCKINDICATOR'   :'custom',
                       },
          'run_opts': {**run_opt_fsp_conf,
                       **run_opt_fsp_nonconf,
                       **run_opt_fsp_SC_first,
                      },

         }
   #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
   run_opt_fsp_SC_TVD_ES={'runs/mhd/freestream/SC_TVD-ES_Fjordholm':
         {'tags': ['mhd','freestream','conforming','SC','TVD_ES'] ,
          'test_opts':{'err_Linf':{'func': check_error ,
                                   'f_kwargs': {'whichError':'L_inf ','err_tol': 1e-11} } ,
                      },
         },
      }
   #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
   caseID=caseID+1
   jobs['mhd_split_glm_noncons_br1_entr_vars_SC']={
          'case': caseID ,
          'tags': ['mhd','split-form','br1','SC','GL','GLM','NONCONS'] ,
          'build_opts':{**baseopts,
                        'FLUXO_DISCTYPE'         :'2',
                        'FLUXO_DISC_NODETYPE'    :'GAUSS-LOBATTO',
                        'FLUXO_EQN_GLM'          :'ON',
                        'FLUXO_EQN_NONCONS'      :'ON',
                        'FLUXO_EQN_NONCONS_GLM'  :'ON',
                        'FLUXO_SHOCKCAPTURE'     :'ON',
                        'FLUXO_SHOCKCAP_NFVSE'   :'ON',
                        'FLUXO_SHOCKINDICATOR'   :'custom',
                        'FLUXO_PARABOLIC'        :'ON',
                        'FLUXO_PARABOLIC_LIFTING':'br1',
                        'FLUXO_PARABOLIC_LIFTING_VAR':'entropy_var',
                       },
          'run_opts': {**run_opt_fsp_conf,
                       **run_opt_fsp_nonconf,
                       **run_opt_fsp_SC_first,
                       **run_opt_fsp_SC_TVD_ES,
                      },
      }
   #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

   #============================================================================
   #============================================================================
   #fourth group, Navierstokes, 400 < caseID <500
   #============================================================================
   #============================================================================
   caseID=400
   baseopts={ 'FLUXO_EQNSYSNAME' :'navierstokes',
              'FLUXO_TESTCASE'   :'default'
            }
   
   run_opt_fsp_conf={'runs/navst/freestream/conforming':
         {'tags': ['navierstokes','freestream','conforming'] ,
          'test_opts':{'err_Linf':{'func': check_error ,'f_kwargs': {'whichError':'L_inf ','err_tol': 1e-11} } ,
                      },
         },
      }
   run_opt_fsp_nonconf={'runs/navst/freestream/nonconforming':
         {'tags': ['navierstokes','freestream','nonconforming'] ,
          'test_opts':{'err_Linf':{'func': check_error ,'f_kwargs': {'whichError':'L_inf ','err_tol': 1e-11} } ,
                      },
         },
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
          'run_opts': {**run_opt_fsp_conf, 
                       **run_opt_fsp_nonconf,
                      }
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
          'run_opts': {**run_opt_fsp_conf, 
                       **run_opt_fsp_nonconf,
                      }
         }
   #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
   caseID=caseID+1
   jobs['build_navierstokes_type1_br1_entropy_var']={
          'case': caseID,
          'tags': [ 'navierstokes','standardDG','GL','br1','entropy_var'],
          'build_opts':{**baseopts,
                        'FLUXO_DISCTYPE'         :'1',
                        'FLUXO_DISC_NODETYPE'    :'GAUSS-LOBATTO',
                        'FLUXO_PARABOLIC'        :'ON',
                        'FLUXO_PARABOLIC_LIFTING':'br1',
                        'FLUXO_PARABOLIC_LIFTING_VAR':'entropy_var',
                       },
          'run_opts': {**run_opt_fsp_conf, 
                       **run_opt_fsp_nonconf,
                      }
         }
   #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
   caseID=caseID+1
   jobs['build_navierstokes_type2_br1_prim_var']={
          'case': caseID,
          'tags': [ 'navierstokes','split-form','GL','br1','prim_var'],
          'build_opts':{**baseopts,
                        'FLUXO_DISCTYPE'         :'2',
                        'FLUXO_DISC_NODETYPE'    :'GAUSS-LOBATTO',
                        'FLUXO_PARABOLIC'        :'ON',
                        'FLUXO_PARABOLIC_LIFTING':'br1',
                        'FLUXO_PARABOLIC_LIFTING_VAR':'prim_var',
                       },
          'run_opts': {**run_opt_fsp_conf, 
                       **run_opt_fsp_nonconf,
                      }
         }
   #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
   caseID=caseID+1
   jobs['build_navierstokes_type1_Gauss_br1_viscsuth']={
          'case': caseID,
          'tags': [ 'navierstokes','standardDG','gauss','br1','cons_var',"sutherland"],
          'build_opts':{**baseopts,
                        'FLUXO_DISCTYPE'         :'1',
                        "FLUXO_EQN_VISCOSITY"    :"sutherland",
                        "FLUXO_DISC_NODETYPE"    :"GAUSS",
                        "FLUXO_PARABOLIC"        :"ON",
                        "FLUXO_PARABOLIC_LIFTING":"br1",
                        "FLUXO_PARABOLIC_LIFTING_VAR":"cons_var",
                       },
          'run_opts': {**run_opt_fsp_conf, 
                       **run_opt_fsp_nonconf,
                      }
         }
   #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
   caseID=caseID+1
   jobs['build_navierstokes_type1_GL_br2_viscpow']={
          'case': caseID,
          'tags': [ 'navierstokes','standardDG','GL','br2','cons_var',"powerlaw"],
          'build_opts':{**baseopts,
                        'FLUXO_DISCTYPE'         :'1',
                        "FLUXO_EQN_VISCOSITY"    :"powerlaw",
                        "FLUXO_DISC_NODETYPE"    :"GAUSS-LOBATTO",
                        "FLUXO_PARABOLIC"        :"ON",
                        "FLUXO_PARABOLIC_LIFTING":"br2",
                        "FLUXO_PARABOLIC_LIFTING_VAR":"cons_var",
                       },
          'run_opts': {**run_opt_fsp_conf, 
                       **run_opt_fsp_nonconf,
                      }
         }
   #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
   caseID=caseID+1
   jobs['build_navierstokes_type1_GL_nopara']={
          'case': caseID,
          'tags': [ 'navierstokes','standardDG','GL'],
          'build_opts':{**baseopts,
                        'FLUXO_DISCTYPE'         :'1',
                        "FLUXO_DISC_NODETYPE"    :"GAUSS-LOBATTO",
                        "FLUXO_PARABOLIC"        :"OFF",
                       },
          'run_opts': {**run_opt_fsp_conf, 
                       **run_opt_fsp_nonconf,
                      }
         }
   #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
   volfluxes=['-1','0','5','8']
   for  vvv  in range(0,len(volfluxes)):
      volflux=volfluxes[vvv]
      caseID=caseID+1
   
      jobs['build_navierstokes_type2_nopara_volFlux_'+volflux]={
            'case': caseID,
            'tags': [ 'navierstokes','split-form','GL'],
            'build_opts':{**baseopts,
                          'FLUXO_DISCTYPE'         :'2',
                          "FLUXO_DISC_NODETYPE"    :"GAUSS-LOBATTO",
                          "FLUXO_PARABOLIC"        :"OFF",
                          'FLUXO_EQN_VOLFLUX'      : volflux,
                         },
            'run_opts': {**run_opt_fsp_conf, 
                         **run_opt_fsp_nonconf,
                        }
           }
   #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
   caseID=caseID+1
   jobs['build_navierstokes_type1_GL_nopara_TC_angmom']={
          'case': caseID,
          'tags': [ 'navierstokes','standardDG','GL'],
          'build_opts':{**baseopts,
                        'FLUXO_DISCTYPE'         :'1',
                        "FLUXO_DISC_NODETYPE"    :"GAUSS-LOBATTO",
                        "FLUXO_PARABOLIC"        :"OFF",
                        "FLUXO_TESTCASE"         :"ns_angularmomentum",
                       },
          'run_opts': {**run_opt_fsp_conf, 
                       **run_opt_fsp_nonconf,
                      }
         }
   #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
   # New test: A shock at a very low density (rho=1e-14) that is captured by the positivity limiter
   run_opt_shock_posit={'runs/navst/shock/positivity':
         {'tags': ['navierstokes','shock','conforming','positivity','SC','firstorder'] ,
          'test_opts':{'err_Linf':{'func': check_error ,
                                   'f_kwargs': {'whichError':'L_inf ','err_tol': 1e-11} } ,
                      },
         },
      }
   #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
   caseID=caseID+1
   jobs['build_navierstokes_type2_nopara_SC_posit']={
          'case': caseID,
          'tags': [ 'navierstokes','split-form','GL','SC','positivity'],
          'build_opts':{**baseopts,
                        'FLUXO_DISCTYPE'         :'2',
                        "FLUXO_DISC_NODETYPE"    :"GAUSS-LOBATTO",
                        "FLUXO_PARABOLIC"        :"OFF",
                        "FLUXO_SHOCKCAPTURE"     :"ON",
                        "FLUXO_SHOCKCAP_NFVSE"   :"ON",
                        "FLUXO_SHOCK_NFVSE_CORR" :"ON",
                       },
          'run_opts': {**run_opt_fsp_conf, 
                       **run_opt_fsp_nonconf,
                       **run_opt_shock_posit,
                      }
         }
   #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
   # New test: A shock at a very low density (rho=1e-14) that is captured by the shock indicator
   run_opt_shock_SC={'runs/navst/shock/SC_firstOrder':
         {'tags': ['navierstokes','shock','conforming','SC','firstorder'] ,
          'test_opts':{'err_Linf':{'func': check_error ,
                                   'f_kwargs': {'whichError':'L_inf ','err_tol': 1e-11} } ,
                      },
         },
      }
   #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
   caseID=caseID+1
   jobs['build_navierstokes_type2_nopara_SC_firstOrder']={
          'case': caseID,
          'tags': [ 'navierstokes','split-form','GL','SC'],
          'build_opts':{**baseopts,
                        'FLUXO_DISCTYPE'         :'2',
                        "FLUXO_DISC_NODETYPE"    :"GAUSS-LOBATTO",
                        "FLUXO_PARABOLIC"        :"OFF",
                        "FLUXO_SHOCKCAPTURE"     :"ON",
                        "FLUXO_SHOCKCAP_NFVSE"   :"ON",
                       },
          'run_opts': {**run_opt_fsp_conf, 
                       **run_opt_fsp_nonconf,
                       **run_opt_shock_SC,
                      }
         }
   #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
   # New test: Free stream for the TVD SC (prim vars, ReconsBoundaries=1)
   run_opt_fsp_SC_TVD={'runs/navst/freestream/SC_TVD':
         {'tags': ['navierstokes','freestream','conforming','SC','TVD'] ,
          'test_opts':{'err_Linf':{'func': check_error ,
                                   'f_kwargs': {'whichError':'L_inf ','err_tol': 1e-11} } ,
                      },
         },
      }
   #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
   # New test: Free stream for the TVD-ES SC (prim vars with entropy fix, ReconsBoundaries=2)
   run_opt_fsp_SC_TVD_ES_fix={'runs/navst/freestream/SC_TVD_ES_fix':
         {'tags': ['navierstokes','freestream','conforming','SC','TVD_ES_fix'] ,
          'test_opts':{'err_Linf':{'func': check_error ,
                                   'f_kwargs': {'whichError':'L_inf ','err_tol': 1e-11} } ,
                      },
         },
      } 
   #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
   # New test: Free stream for the TVD-ES SC à la Fjordholm (ReconsBoundaries=3)
   run_opt_fsp_SC_TVD_ES_Fjordholm={'runs/navst/freestream/SC_TVD_ES_Fjordholm':
         {'tags': ['navierstokes','freestream','conforming','SC','TVD_ES_Fjod'] ,
          'test_opts':{'err_Linf':{'func': check_error ,
                                   'f_kwargs': {'whichError':'L_inf ','err_tol': 1e-11} } ,
                      },
         },
      } 
   #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
   caseID=caseID+1
   jobs['build_navierstokes_type2_br1_entr_vars_SC_TVD']={
          'case': caseID,
          'tags': [ 'navierstokes','split-form','GL','br1','entropy_var','SC'],
          'build_opts':{**baseopts,
                        'FLUXO_DISCTYPE'         :'2',
                        "FLUXO_DISC_NODETYPE"    :"GAUSS-LOBATTO",
                        "FLUXO_PARABOLIC"        :"ON",
                        "FLUXO_PARABOLIC_LIFTING":"br1",
                        "FLUXO_PARABOLIC_LIFTING_VAR":"entropy_var",
                        "FLUXO_SHOCKCAPTURE"     :"ON",
                        "FLUXO_SHOCKCAP_NFVSE"   :"ON",
                       },
          'run_opts': {**run_opt_fsp_conf, 
                       **run_opt_fsp_nonconf,
                       **run_opt_fsp_SC_TVD,
                       **run_opt_fsp_SC_TVD_ES_fix,
                       **run_opt_fsp_SC_TVD_ES_Fjordholm,
                      }
         }
   #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
   # New test: Free stream for the TVD-ES SC à la Fjordholm (ReconsBoundaries=3)
   run_opt_fsp_p4est={'runs/navst/freestream/p4est_1':
         {'tags': ['navierstokes','freestream','nonconforming','amr','p4est'] ,
          'test_opts':{'err_Linf':{'func': check_error ,
                                   'f_kwargs': {'whichError':'L_inf ','err_tol': 1e-11} } ,
                      },
         },
      } 
   #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
   caseID=caseID+1
   jobs['build_navierstokes_p4est']={
          'case': caseID,
          'tags': [ 'navierstokes','split-form','GL','amr','p4est'],
          'build_opts':{**baseopts,
                        'FLUXO_DISCTYPE'         :'2',
                        "FLUXO_DISC_NODETYPE"    :"GAUSS-LOBATTO",
                        "FLUXO_PARABOLIC"        :"OFF",
                        "FLUXO_AMR"              :"ON",
                        "FLUXO_BUILD_P4EST"      :"OFF",
                       },
          'run_opts': {**run_opt_fsp_conf, 
                       **run_opt_fsp_nonconf,
                       **run_opt_fsp_p4est,
                      }
         }
   return jobs
