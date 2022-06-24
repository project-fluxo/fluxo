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
#==========================================================================
# open the file (default "std.out") and parse for `findstr`, 
# returns all the lines of the file that contain `findstr`
#==========================================================================
def find_all_occurences(findstr,stdout_filepath="std.out",**kwargs):
   with open(stdout_filepath,'r') as fp: 
      lines=fp.readlines()
   linesWithStr=[]
   for line in reversed(lines) :
      if findstr in line : 
         linesWithStr.append(line)
   return linesWithStr


###########################################################################
#  ROUTINES FOR `test_opts`: 
#    * MUST have only keyword arguments (not positional!). 
#    * MUST always have the container argument `**kwargs` as last argument!
#    * MUST return status (True/False) and message
#    * ASSUMES current directory is the one where the code was executed!
###########################################################################

#==========================================================================
# find the line of the string `whichError` and check all entries after :  to_be < err_tol
#==========================================================================
def check_error(whichError='L_inf ', err_tol=1.0e-12,err_abs=True, to_be=[0.0],**kwargs ):
   [found,line]=find_last_occurence(whichError,**kwargs)
   assert found, ('check_error: did not find "%s" in std.out' % (whichError))
   errors=[float(x) for x in (line.split(":")[1]).split()]
   if(len(to_be)==1):
     errors_ref=[e-to_be[0] for e in errors]
   else:
     assert len(errors)==len(to_be), ('to_be must either be a scalar of have the same length as errors')
     errors_ref=[e-ref for (e,ref) in zip(errors,to_be)]
   if err_abs:
     errors_ref = [abs(e) for e in errors_ref]
   check=all([e < err_tol for e in errors_ref])
   msg=('check_error: "%s" - to_be < %e = %s' % (line.strip(),err_tol,check))
   return check,msg
   
#==========================================================================
# find all the lines of the string `whichError` and check all entries after :  to be < err_tol
#==========================================================================
def check_all_errors(whichError='dSdU*Ut', err_tol=1.0e-12,err_abs=True,**kwargs ):
   lines=find_all_occurences(whichError,**kwargs)
   assert len(lines)>=1, ('check_all_errors: did not find "%s" in std.out' % (whichError))
   errors = []
   check = True
   for line in lines:
      newErrors=[float(x) for x in (line.split(":")[1]).split()]
      if err_abs:
        newErrors = [abs(e) for e in newErrors]
      errors.append (newErrors)
      check=all([e < err_tol for e in newErrors]) and check
   msg=('check_all_errors: "%s" < %e = %s' % (max(errors),err_tol,check))
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
   run_opt_fsp_conf={'runs/linadv/freestream/conforming': 
                 {'tags': ['linadv','freestream','curved','conforming'],
                  'test_opts':{ 'err_Linf':{'func': check_error ,
                                            'f_kwargs': {'whichError':'L_inf ','err_tol': 1e-12} },
                              },
                 },
             }
   run_opt_fsp_nonconf_coll={'runs/linadv/freestream/nonconforming_collmortar': 
                 {'tags': ['linadv','freestream','curved','nonconforming','collocation-mortar'],
                  'test_opts':{ 'err_Linf':{'func': check_error ,
                                            'f_kwargs': {'whichError':'L_inf ','err_tol': 1e-12} },
                              },
                 },
             }
   run_opt_fsp_nonconf_proj={'runs/linadv/freestream/nonconforming_projmortar': 
                 {'tags': ['linadv','freestream','curved','nonconforming','projection-mortar'],
                  'test_opts':{ 'err_Linf':{'func': check_error ,
                                            'f_kwargs': {'whichError':'L_inf ','err_tol': 1e-12} },
                              },
                 },
             }
   run_opt_linfunc_conf={'runs/linadv/linfunc/conforming': 
                 {'tags': ['linadv','linfunc','curved','conforming'],
                  'test_opts':{ 'err_Linf':{'func': check_error ,
                                            'f_kwargs': {'whichError':'L_inf ','err_tol': 1e-12} },
                              },
                 },
             }
   run_opt_linfunc_nonconf_coll={'runs/linadv/linfunc/nonconforming_collmortar': 
                 {'tags': ['linadv','linfunc','curved','nonconforming','collocation-mortar'],
                  'test_opts':{ 'err_Linf':{'func': check_error ,
                                            'f_kwargs': {'whichError':'L_inf ','err_tol': 5e-12} },
                              },
                 },
             }
   run_opt_linfunc_nonconf_proj={'runs/linadv/linfunc/nonconforming_projmortar': 
                 {'tags': ['linadv','linfunc','curved','nonconforming','projection-mortar'],
                  'test_opts':{ 'err_Linf':{'func': check_error ,
                                            'f_kwargs': {'whichError':'L_inf ','err_tol': 5e-12} },
                              },
                 },
             }
   #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
   # TEST FOR JESSE_MORTAR: (CENTRAL FLUX+CURVED+MORTAR ONLY FSP with JESSE_MORTAR)
   run_opt_fsp_nonconf_coll_central={'runs/linadv/freestream/nonconforming_collmortar_central': 
                 {'tags': ['linadv','freestream','curved','nonconforming','collocation-mortar','centralflux'],
                  'test_opts':{ 'err_Linf':{'func': check_error ,
                                            'f_kwargs': {'whichError':'L_inf ','err_tol': 1e-10} },
                              },
                 },
             }
   run_opt_linfunc_nonconf_coll_central={'runs/linadv/linfunc/nonconforming_collmortar_central': 
                 {'tags': ['linadv','linfunc','curved','nonconforming','collocation-mortar','centralflux'],
                  'test_opts':{ 'err_Linf':{'func': check_error ,
                                            'f_kwargs': {'whichError':'L_inf ','err_tol': 1e-10} },
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
          'run_opts': {**run_opt_fsp_conf, 
                       **run_opt_fsp_nonconf_proj,
                       **run_opt_linfunc_conf,
                       **run_opt_linfunc_nonconf_proj,
                      }
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
          'run_opts': {**run_opt_fsp_conf, 
                       **run_opt_fsp_nonconf_coll,
                       **run_opt_linfunc_conf,
                       **run_opt_linfunc_nonconf_coll,
                      }
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
          'run_opts': {**run_opt_fsp_conf, 
                       **run_opt_fsp_nonconf_coll,
                       **run_opt_linfunc_conf,
                       **run_opt_linfunc_nonconf_coll,
                      }
         }
   #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
   caseID=caseID+1
   jobs['linadv_nopara_type1_gauss']={
          'case': caseID,
          'tags': ['linadv','standardDG','gauss'],
          'build_opts':{**baseopts,
                        'FLUXO_DISCTYPE'         :'1',
                        'FLUXO_DISC_NODETYPE'    :'GAUSS',
                        'FLUXO_PARABOLIC'        :'OFF',
                       },
          'run_opts': {**run_opt_fsp_conf, 
                       **run_opt_fsp_nonconf_coll,
                       **run_opt_linfunc_conf,
                       **run_opt_linfunc_nonconf_coll,
                      }
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
          'run_opts': {**run_opt_fsp_conf, 
                       **run_opt_fsp_nonconf_proj,
                       **run_opt_linfunc_conf,
                       **run_opt_linfunc_nonconf_proj,
                      }
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
          'run_opts': {**run_opt_fsp_conf, 
                       **run_opt_fsp_nonconf_coll,
                       **run_opt_linfunc_conf,
                       **run_opt_linfunc_nonconf_coll,    
                      }
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
          'run_opts': {**run_opt_fsp_conf, 
                       **run_opt_linfunc_conf,
                      }
         }
   #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
   caseID=caseID+1
   jobs['linadv_split_jesse_mortar']={
          'case': caseID,
          'tags': ['linadv','split-form','GL','jesse-mortar'],
          'build_opts':{**baseopts,
                        'FLUXO_DISCTYPE'         :'2',
                        'FLUXO_DISC_NODETYPE'    :'GAUSS-LOBATTO',
                        'FLUXO_JESSE_MORTAR'     :'ON',
                       },
          'run_opts': {**run_opt_fsp_conf, 
                       **run_opt_fsp_nonconf_coll,
                       **run_opt_fsp_nonconf_coll_central,
                       **run_opt_linfunc_conf,
                       **run_opt_linfunc_nonconf_coll,
                       **run_opt_linfunc_nonconf_coll_central,
                      }
         }
   #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
   caseID=caseID+1
   jobs['linadv_split_br1_jesse_mortar']={
          'case': caseID,
          'tags': ['linadv','split-form','GL','br1','jesse-mortar'],
          'build_opts':{**baseopts,
                        'FLUXO_DISCTYPE'         :'2',
                        'FLUXO_DISC_NODETYPE'    :'GAUSS-LOBATTO',
                        'FLUXO_JESSE_MORTAR'     :'ON',
                        'FLUXO_PARABOLIC'        :'ON',
                        'FLUXO_PARABOLIC_LIFTING':'br1',
                       },
          'run_opts': {**run_opt_fsp_conf, 
                       **run_opt_fsp_nonconf_coll,
                       **run_opt_fsp_nonconf_coll_central,
                       **run_opt_linfunc_conf,
                       **run_opt_linfunc_nonconf_coll,
                       **run_opt_linfunc_nonconf_coll_central,
                      }
         }
   #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
   caseID=caseID+1
   jobs['linadv_split_gauss_br1_jesse_mortar']={
          'case': caseID,
          'tags': ['linadv','split-form','gauss','br1','jesse-mortar'],
          'build_opts':{**baseopts,
                        'FLUXO_DISCTYPE'         :'2',
                        'FLUXO_DISC_NODETYPE'    :'GAUSS',
                        'FLUXO_JESSE_MORTAR'     :'ON',
                        'FLUXO_PARABOLIC'        :'ON',
                        'FLUXO_PARABOLIC_LIFTING':'br1',
                       },
          'run_opts': {**run_opt_fsp_conf, 
                       **run_opt_fsp_nonconf_coll,
                       **run_opt_fsp_nonconf_coll_central,
                       **run_opt_linfunc_conf,
                       **run_opt_linfunc_nonconf_coll,
                       **run_opt_linfunc_nonconf_coll_central,
                      }
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
   run_opt_fsp_conf={'runs/maxwell/freestream/conforming': 
                    {'tags': ['maxwell','freestream','curved','conforming'],
                     'test_opts':{ 'err_Linf':{'func': check_error ,
                                               'f_kwargs': {'whichError':'L_inf ','err_tol': 1e-11} },
                                 },
                    },
                  }
   run_opt_fsp_conf_central={'runs/maxwell/freestream/conforming_central': 
                    {'tags': ['maxwell','freestream','curved','conforming','centralflux'],
                     'test_opts':{ 'err_Linf':{'func': check_error ,
                                               'f_kwargs': {'whichError':'L_inf ','err_tol': 1e-11} },
                                 },
                    },
                  }
   run_opt_fsp_nonconf_proj={'runs/maxwell/freestream/nonconforming_projmortar': 
                    {'tags': ['maxwell','freestream','curved','nonconforming','projection-mortar'],
                     'test_opts':{ 'err_Linf':{'func': check_error ,
                                               'f_kwargs': {'whichError':'L_inf ','err_tol': 1e-12} },
                                 },
                    },
                  }
   run_opt_fsp_nonconf_coll={'runs/maxwell/freestream/nonconforming_collmortar': 
                    {'tags': ['maxwell','freestream','curved','nonconforming','collocation-mortar'],
                     'test_opts':{ 'err_Linf':{'func': check_error ,
                                               'f_kwargs': {'whichError':'L_inf ','err_tol': 1e-11} },
                                 },
                    },
                  }
   #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
   # TEST FOR JESSE_MORTAR: (CENTRAL FLUX+CURVED+MORTAR ONLY FSP with JESSE_MORTAR)
   run_opt_fsp_nonconf_coll_central={'runs/maxwell/freestream/nonconforming_collmortar_central': 
                    {'tags': ['maxwell','freestream','curved','nonconforming','collocation-mortar','centralflux'],
                     'test_opts':{ 'err_Linf':{'func': check_error ,
                                               'f_kwargs': {'whichError':'L_inf ','err_tol': 1e-11} },
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
          'run_opts': {**run_opt_fsp_conf,
                       **run_opt_fsp_nonconf_proj,
                      }
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
          'run_opts': {**run_opt_fsp_conf,
                       **run_opt_fsp_conf_central,
                       **run_opt_fsp_nonconf_coll,
                      }
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
          'run_opts': {**run_opt_fsp_conf,
                       **run_opt_fsp_nonconf_coll,
                      }
         }
   #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
   volfluxes=['-1','0','1']
   for  vvv  in range(0,len(volfluxes)):
      caseID=caseID+1
      volflux=volfluxes[vvv]
      jobs['maxwell_split_volflux_'+volflux]={
          'case': caseID,
          'tags': ['maxwell','split-form','GL'],
          'build_opts':{**baseopts,
                        'FLUXO_DISCTYPE'         :'2',
                        'FLUXO_DISC_NODETYPE'    :'GAUSS-LOBATTO',
                        'FLUXO_EQN_VOLFLUX'      :volflux,
                       },
          'run_opts': {**run_opt_fsp_conf,}
         }
   #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
   caseID=caseID+1
   jobs['maxwell_split_jesse_mortar']={
          'case': caseID,
          'tags': ['maxwell','split-form','GL','jesse-mortar'],
          'build_opts':{**baseopts,
                        'FLUXO_DISCTYPE'         :'2',
                        'FLUXO_DISC_NODETYPE'    :'GAUSS-LOBATTO',
                        'FLUXO_JESSE_MORTAR'     :'ON',
                       },
          'run_opts': {**run_opt_fsp_conf, 
                       **run_opt_fsp_conf_central,
                       **run_opt_fsp_nonconf_coll,
                       **run_opt_fsp_nonconf_coll_central,
                      }
         }
   #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
   caseID=caseID+1
   jobs['maxwell_split_gauss_jesse_mortar']={
          'case': caseID,
          'tags': ['maxwell','split-form','gauss','jesse-mortar'],
          'build_opts':{**baseopts,
                        'FLUXO_DISCTYPE'         :'2',
                        'FLUXO_DISC_NODETYPE'    :'GAUSS',
                        'FLUXO_JESSE_MORTAR'     :'ON',
                       },
          'run_opts': {**run_opt_fsp_conf, 
                       **run_opt_fsp_conf_central,
                       **run_opt_fsp_nonconf_coll,
                       **run_opt_fsp_nonconf_coll_central,
                      }
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
         {'tags': ['mhd','freestream','curved','conforming'] ,
          'test_opts':{'err_Linf':{'func': check_error ,
                                   'f_kwargs': {'whichError':'L_inf ','err_tol': 1e-11} } ,
                      },
         },
      }
   run_opt_fsp_conf_multicomponent={'runs/mhd/freestream/conforming_multicomponent':
         {'tags': ['mhd','multicomponent','freestream','curved','conforming'] ,
          'test_opts':{'err_Linf':{'func': check_error ,
                                   'f_kwargs': {'whichError':'L_inf ','err_tol': 1e-11} } ,
                      },
         },
      }
   run_opt_fsp_nonconf_proj={'runs/mhd/freestream/nonconforming_projmortar':
         {'tags': ['mhd','freestream','curved','nonconforming','projection-mortar'] ,
          'test_opts':{'err_Linf':{'func': check_error ,
                                   'f_kwargs': {'whichError':'L_inf ','err_tol': 1e-11} } ,
                      },
         },
      }
   run_opt_fsp_nonconf_coll={'runs/mhd/freestream/nonconforming_collmortar':
         {'tags': ['mhd','freestream','curved','nonconforming','collocation-mortar'] ,
          'test_opts':{'err_Linf':{'func': check_error ,
                                   'f_kwargs': {'whichError':'L_inf ','err_tol': 1e-11} } ,
                      },
         },
      }
   # Entropy conservation test with EC flux (with and without shock-capturing)
   run_opt_entropyCons={'runs/mhd/softBlast/entropyCons':
         {'tags': ['mhd','entropyCons','curved','conforming'] ,
          'test_opts':{'abs(dSdU*Ut)':{'func': check_all_errors ,
                                       'f_kwargs': {'whichError':'dSdU*Ut','err_tol': 1e-13,'err_abs':True} } ,
                      },
         },
      }
   # Entropy conservation test with EC flux for multi-component MHD (with and without shock-capturing)
   run_opt_entropyCons_multicomponent={'runs/mhd/softBlast/entropyCons_multicomponent':
         {'tags': ['mhd','multicomponent','entropyCons','curved','conforming'] ,
          'test_opts':{'abs(dSdU*Ut)':{'func': check_all_errors ,
                                       'f_kwargs': {'whichError':'dSdU*Ut','err_tol': 1e-11,'err_abs':True} } ,
                      },
         },
      }
   # Entropy conservation test with EC flux (needs coll. + Jesse_mortars!)
   run_opt_entropyCons_nonconf={'runs/mhd/softBlast/entropyCons_nonconf':
         {'tags': ['mhd','entropyCons','curved','nonconforming'] ,
          'test_opts':{'abs(dSdU*Ut)':{'func': check_all_errors ,
                                       'f_kwargs': {'whichError':'dSdU*Ut','err_tol': 1e-13,'err_abs':True} } ,
                      },
         },
      } 
   # Entropy stability test with EC flux and LLF (with and without shock-capturing)
   run_opt_entropyStab={'runs/mhd/softBlast/entropyStab':
         {'tags': ['mhd','entropyStab','curved','conforming'] ,
          'test_opts':{'dSdU*Ut':{'func': check_all_errors ,
                                  'f_kwargs': {'whichError':'dSdU*Ut','err_tol': 1e-13,'err_abs':False} } ,
                      },
         },
      }  
   run_opt_entropyStab_nonconf={'runs/mhd/softBlast/entropyStab_nonconf':
         {'tags': ['mhd','entropyStab','curved','nonconforming'] ,
          'test_opts':{'dSdU*Ut':{'func': check_all_errors ,
                                  'f_kwargs': {'whichError':'dSdU*Ut','err_tol': 1e-13,'err_abs':False} } ,
                      },
         },
      }  
   # Entropy stability test with EC flux and LLF (with and without shock-capturing)
   run_opt_entropyStab_FloGor9waves={'runs/mhd/softBlast/entropyStab_FloGor9waves':
         {'tags': ['mhd','entropyStab','curved','conforming','FloGor'] ,
          'test_opts':{'dSdU*Ut':{'func': check_all_errors ,
                                  'f_kwargs': {'whichError':'dSdU*Ut','err_tol': 1e-13,'err_abs':False} } ,
                      },
         },
      }  
   # Entropy conservation test with EC flux, SC, and AMR (needs coll. + Jesse_mortars!)
   run_opt_entropyCons_AMR_SC={'runs/mhd/softBlast/entropyCons_AMR':
         {'tags': ['mhd','entropyCons','curved','nonconforming','amr','SC'] ,
          'test_opts':{'abs(dSdU*Ut)':{'func': check_all_errors ,
                                       'f_kwargs': {'whichError':'dSdU*Ut','err_tol': 1e-13,'err_abs':True} } ,
                      },
         },
      }
   # Entropy stability test with EC flux + full wave dissipation, reconstructed FV-SC, and AMR (needs coll. + Jesse_mortars!)
   run_opt_entropyStab_AMR_SC={'runs/mhd/softBlast/entropyStab_AMR':
         {'tags': ['mhd','entropyStab','curved','nonconforming','amr','SC'] ,
          'test_opts':{'abs(dSdU*Ut)':{'func': check_all_errors ,
                                       'f_kwargs': {'whichError':'dSdU*Ut','err_tol': 1e-13,'err_abs':False} } ,
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
                       **run_opt_fsp_nonconf_proj,
                      }
         }
   #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
   caseID=caseID+1
   jobs['mhd_split_glm_noncons_br1entr_ecvolflux']={
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
                        'FLUXO_EQN_VOLFLUX'      :'10',
                       },
          'run_opts': {**run_opt_fsp_conf, 
                       **run_opt_fsp_nonconf_coll,
                       **run_opt_entropyCons,
                       **run_opt_entropyStab,
                      }
         }
   #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
   caseID=caseID+1
   jobs['mhd_split_glm_noncons_br1entr_jesse']={
          'case': caseID ,
          'tags': ['mhd','split-form','br1','GL','GLM','NONCONS','ECflux','jesse-mortar'] ,
          'build_opts':{**baseopts,
                        'FLUXO_DISCTYPE'         :'2',
                        'FLUXO_DISC_NODETYPE'    :'GAUSS-LOBATTO',
                        'FLUXO_EQN_GLM'          :'ON',
                        'FLUXO_EQN_NONCONS'      :'ON',
                        'FLUXO_EQN_NONCONS_GLM'  :'ON',
                        'FLUXO_PARABOLIC'        :'ON',
                        'FLUXO_PARABOLIC_LIFTING':'br1',
                        'FLUXO_PARABOLIC_LIFTING_VAR':'entropy_var',
                        'FLUXO_JESSE_MORTAR'     :'ON',
                       },
          'run_opts': {**run_opt_fsp_conf, 
                       **run_opt_fsp_nonconf_coll,
                       **run_opt_entropyCons,
                       **run_opt_entropyCons_nonconf,
                       **run_opt_entropyStab,
                       **run_opt_entropyStab_nonconf,
                      }
         }
   #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
   caseID=caseID+1
   jobs['mhd_split_glm_noncons_br1entr_ecvolflux_testcase_mhdeq']={
          'case': caseID ,
          'tags': ['mhd','split-form','br1','GL','GLM','NONCONS','ECflux','FloGor'] ,
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
          'run_opts': {
                       **run_opt_fsp_conf, 
                       **run_opt_fsp_nonconf_coll,
                       **run_opt_entropyCons,
                       **run_opt_entropyStab,
                       **run_opt_entropyStab_FloGor9waves,
                      }
         }
   #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
   caseID=caseID+1
   jobs['mhd_split_noglm_noncons_nopara_testcase_tgv']={
          'case': caseID ,
          'tags': ['mhd','split-form','GL','NONCONS','ECflux','jesse-mortar'] ,
          'build_opts':{**baseopts,
                        'FLUXO_DISCTYPE'         :'2',
                        'FLUXO_DISC_NODETYPE'    :'GAUSS-LOBATTO',
                        'FLUXO_EQN_GLM'          :'OFF',
                        'FLUXO_EQN_NONCONS'      :'ON',
                        'FLUXO_PARABOLIC'        :'OFF',
                        'FLUXO_EQN_VOLFLUX'      :'10',
                        'FLUXO_TESTCASE'         :'taylorgreenvortex',
                        'FLUXO_JESSE_MORTAR'     :'ON',                 # We need Jesse mortars to achieve FSP with LGL+noncons(no GLM)
                       },
          'run_opts': {**run_opt_entropyCons,
                       **run_opt_fsp_conf, 
                       **run_opt_fsp_nonconf_coll,
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
                       **run_opt_fsp_nonconf_coll,
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
                       **run_opt_fsp_nonconf_coll,
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
          'run_opts': {**run_opt_entropyCons,
                       **run_opt_entropyStab,
                       **run_opt_fsp_conf, 
                       **run_opt_fsp_nonconf_coll,
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
                       **run_opt_fsp_nonconf_proj,
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
                       **run_opt_fsp_nonconf_proj,
                      }
         }
   #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
   volfluxes=['10','12']
   for  vvv  in range(0,len(volfluxes)):
      volflux=volfluxes[vvv]
      caseID=caseID+1
      
      jobs['mhd_split_noglm_noncons_nopara_volflux_'+volflux]={
          'case': caseID ,
          'tags': ['mhd','split-form','GL','NONCONS','jesse-mortar'] ,
          'build_opts':{**baseopts,
                        'FLUXO_DISCTYPE'         :'2',
                        'FLUXO_DISC_NODETYPE'    :'GAUSS-LOBATTO',
                        'FLUXO_EQN_GLM'          :'ON',
                        'FLUXO_EQN_NONCONS'      :'ON',
                        'FLUXO_PARABOLIC'        :'OFF',
                        'FLUXO_EQN_VOLFLUX'      : volflux,
                        'FLUXO_JESSE_MORTAR'     :'ON',
                       },
          'run_opts': {**run_opt_fsp_conf, 
                       **run_opt_fsp_nonconf_coll,
                       **run_opt_entropyCons,
                       **run_opt_entropyCons_nonconf,
                       **run_opt_entropyStab,
                       **run_opt_entropyStab_nonconf,
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
                       **run_opt_fsp_nonconf_coll,
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
          'run_opts': {
                       **run_opt_fsp_conf,
                       **run_opt_fsp_nonconf_coll,
                       **run_opt_fsp_SC_first,
                       **run_opt_entropyCons,
                       **run_opt_entropyStab,
                       **run_opt_entropyStab_FloGor9waves,
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
                       **run_opt_fsp_nonconf_coll,
                       **run_opt_fsp_SC_first,
                       **run_opt_fsp_SC_TVD_ES,
                       **run_opt_entropyCons,
                       **run_opt_entropyStab,
                      },
      }
   #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
   # New test: Soft blast with shock capturing and AMR
   run_opt_p4est_SC={'runs/mhd/softBlast/p4est_SC':
         {'tags': ['mhd','blast','nonconforming','p4est','amr','SC'] ,
          'test_opts':{'max|Ut|':{'func': check_error ,
                                  'f_kwargs': {'whichError':'max|Ut| ',
                                                    'to_be': [1.262018175963E+00,   2.523291459432E-01,   5.001874316192E-01,   2.534612844799E-01,   1.333612422134E-01,   3.295500842832E-02,   3.605822666281E-02,   5.257420500555E-02,   5.354065204360E-02],
                                                  'err_tol': 1e-11} } , # err_tol is limited by the precision of the output..
                      },
         },
      }
   # New test: Soft blast with shock capturing and AMR. This test must be run after run_opt_p4est_SC becasue it restarts from an intermediate state of the other test
   run_opt_p4est_SC_restart={'runs/mhd/softBlast/p4est_SC_restart':
         {'tags': ['mhd','blast','nonconforming','p4est','amr','SC'] ,
          'restartfile': '../p4est_SC/MHD_ENTROPYCONS_State_0000000.400000000.h5',
          'test_opts':{'max|Ut|':{'func': check_error ,
                                  'f_kwargs': {'whichError':'max|Ut| ',
                                                    'to_be': [1.262018175963E+00,   2.523291459432E-01,   5.001874316192E-01,   2.534612844799E-01,   1.333612422134E-01,   3.295500842832E-02,   3.605822666281E-02,   5.257420500555E-02,   5.354065204360E-02],
                                                  'err_tol': 1e-11} } , # err_tol is limited by the precision of the output..
                      },
         },
      }
   #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
   caseID=caseID+1
   jobs['mhd_split_glm_noncons_nopara_p4est_SC']={
          'case': caseID ,
          'tags': ['mhd','split-form','SC','GL','GLM','NONCONS','p4est','amr'] ,
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
                        "FLUXO_AMR"              :"ON",
                        "FLUXO_BUILD_P4EST"      :"OFF",
                       },
          'run_opts': {
                       **run_opt_fsp_conf,
                       **run_opt_fsp_nonconf_coll,
                       **run_opt_fsp_SC_first,
                       **run_opt_p4est_SC,
                       **run_opt_p4est_SC_restart,
                      },

         }
   #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
   # New test: Soft blast with shock capturing, AMR and BR1
   run_opt_p4est_SC_br1={'runs/mhd/softBlast/p4est_SC':
         {'tags': ['mhd','blast','nonconforming','p4est','amr','SC','br1'] ,
          'test_opts':{'max|Ut|':{'func': check_error ,
                                  'f_kwargs': {'whichError':'max|Ut| ',
                                                    'to_be': [1.262012788899E+00,   2.523279414674E-01,   5.001853694600E-01,   2.534600694872E-01,   1.333608373794E-01,   3.295501210640E-02,   3.605822133723E-02,   5.257421412853E-02,   5.354058828700E-02],
                                                  'err_tol': 1e-11} } , # err_tol is limited by the precision of the output..
                      },
         },
      }
   #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
   caseID=caseID+1
   jobs['mhd_split_glm_noncons_br1_p4est_SC']={
          'case': caseID ,
          'tags': ['mhd','split-form','SC','GL','GLM','NONCONS','p4est','amr','br1'] ,
          'build_opts':{**baseopts,
                        'FLUXO_DISCTYPE'         :'2',
                        'FLUXO_DISC_NODETYPE'    :'GAUSS-LOBATTO',
                        'FLUXO_EQN_GLM'          :'ON',
                        'FLUXO_EQN_NONCONS'      :'ON',
                        'FLUXO_EQN_NONCONS_GLM'  :'ON',
                        'FLUXO_PARABOLIC'        :'ON',
                        'FLUXO_PARABOLIC_LIFTING':'br1',
                        'FLUXO_PARABOLIC_LIFTING_VAR':'entropy_var',
                        'FLUXO_SHOCKCAPTURE'     :'ON',
                        'FLUXO_SHOCKCAP_NFVSE'   :'ON',
                        'FLUXO_SHOCKINDICATOR'   :'custom',
                        "FLUXO_AMR"              :"ON",
                        "FLUXO_BUILD_P4EST"      :"OFF",
                       },
          'run_opts': {
                       **run_opt_fsp_conf,
                       **run_opt_fsp_nonconf_coll,
                       **run_opt_fsp_SC_first,
                       **run_opt_p4est_SC_br1,
                      },

         }
    #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
   caseID=caseID+1
   jobs['mhd_split_glm_noncons_nopara_p4est_SC_jesse']={
          'case': caseID ,
          'tags': ['mhd','split-form','SC','GL','GLM','NONCONS','p4est','amr','jesse-mortar'] ,
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
                        "FLUXO_AMR"              :"ON",
                        "FLUXO_BUILD_P4EST"      :"OFF",
                        'FLUXO_JESSE_MORTAR'     :'ON',
                       },
          'run_opts': {
                       **run_opt_fsp_conf,
                       **run_opt_fsp_nonconf_coll,
                       **run_opt_fsp_SC_first,
                       **run_opt_entropyCons_AMR_SC,
                       **run_opt_entropyStab_AMR_SC,
                      },

         }
   #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
   # ES Gauss tests
   #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
   caseID=caseID+1
   jobs['mhd_esgauss_noglm_noncons_nopara_testcase_tgv']={
          'case': caseID ,
          'tags': ['mhd','split-form','gauss','esgauss','NONCONS','ECflux'] ,
          'build_opts':{**baseopts,
                        'FLUXO_DISCTYPE'         :'2',
                        'FLUXO_DISC_NODETYPE'    :'GAUSS',
                        'FLUXO_EQN_GLM'          :'OFF',
                        'FLUXO_EQN_NONCONS'      :'ON',
                        'FLUXO_PARABOLIC'        :'OFF',
                        'FLUXO_EQN_VOLFLUX'      :'10',
                        'FLUXO_TESTCASE'         :'taylorgreenvortex',
                       },
          'run_opts': {**run_opt_entropyCons,
                       **run_opt_fsp_conf, 
                       **run_opt_fsp_nonconf_coll,
                      }
         }
   #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
   volfluxes=['10','12']
   for  vvv  in range(0,len(volfluxes)):
      volflux=volfluxes[vvv]
      caseID=caseID+1
      
      jobs['mhd_esgauss_glm_noncons_nopara_volflux_'+volflux]={
          'case': caseID ,
          'tags': ['mhd','split-form','gauss','esgauss','NONCONS'] ,
          'build_opts':{**baseopts,
                        'FLUXO_DISCTYPE'         :'2',
                        'FLUXO_DISC_NODETYPE'    :'GAUSS',
                        'FLUXO_EQN_GLM'          :'ON',
                        'FLUXO_EQN_NONCONS'      :'ON',
                        'FLUXO_PARABOLIC'        :'OFF',
                        'FLUXO_EQN_VOLFLUX'      : volflux,
                        'FLUXO_JESSE_MORTAR'     :'ON',
                       },
          'run_opts': {**run_opt_fsp_conf, 
                       **run_opt_fsp_nonconf_coll,
                       **run_opt_entropyCons,
                       **run_opt_entropyCons_nonconf,
                       **run_opt_entropyStab,
                       **run_opt_entropyStab_nonconf,
                      }
         }
   #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
   caseID=caseID+1
   jobs['mhd_esgauss_glm_no-noncons_nopara']={
          'case': caseID ,
          'tags': ['mhd','split-form','gauss','esgauss','GLM'] ,
          'build_opts':{**baseopts,
                        'FLUXO_DISCTYPE'         :'2',
                        'FLUXO_DISC_NODETYPE'    :'GAUSS',
                        'FLUXO_EQN_GLM'          :'ON',
                        'FLUXO_EQN_NONCONS'      :'OFF',
                        'FLUXO_PARABOLIC'        :'OFF',
                       },
          'run_opts': {**run_opt_fsp_conf, 
                       **run_opt_fsp_nonconf_coll,
                      }
         }
   #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
   caseID=caseID+1
   jobs['mhd_esgauss_glm_noncons_nopara_p4est_jesse']={
          'case': caseID ,
          'tags': ['mhd','split-form','gauss','esgauss','GLM','NONCONS','p4est','amr','jesse-mortar'] ,
          'build_opts':{**baseopts,
                        'FLUXO_DISCTYPE'         :'2',
                        'FLUXO_DISC_NODETYPE'    :'GAUSS',
                        'FLUXO_EQN_GLM'          :'ON',
                        'FLUXO_EQN_NONCONS'      :'ON',
                        'FLUXO_EQN_NONCONS_GLM'  :'ON',
                        'FLUXO_PARABOLIC'        :'OFF',
                        "FLUXO_AMR"              :"ON",
                        "FLUXO_BUILD_P4EST"      :"OFF",
                        'FLUXO_JESSE_MORTAR'     :'ON',
                       },
          'run_opts': {
                       **run_opt_fsp_conf,
                       **run_opt_fsp_nonconf_coll,
                       **run_opt_entropyCons_AMR_SC,
                       **run_opt_entropyStab_AMR_SC,
                      },

         }
   #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
   # Multi-component MHD tests
   #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
   # Split-form DG
   caseID=caseID+1
   jobs['mhd_split_glm_noncons_br1entr_ecvolflux']={
          'case': caseID ,
          'tags': ['mhd','multicomponent','split-form','br1','GL','GLM','NONCONS','ECflux'] ,
          'build_opts':{**baseopts,
                        'FLUXO_DISCTYPE'             :'2',
                        'FLUXO_EQN_NUM_COMPONENTS'   :'2',
                        'FLUXO_DISC_NODETYPE'        :'GAUSS-LOBATTO',
                        'FLUXO_EQN_GLM'              :'ON',
                        'FLUXO_EQN_NONCONS'          :'ON',
                        'FLUXO_EQN_NONCONS_GLM'      :'ON',
                        'FLUXO_PARABOLIC'            :'ON',
                        'FLUXO_PARABOLIC_LIFTING'    :'br1',
                        'FLUXO_PARABOLIC_LIFTING_VAR':'entropy_var',
                        'FLUXO_EQN_VOLFLUX'          :'10',
                       },
          'run_opts': {**run_opt_fsp_conf_multicomponent, 
                       #**run_opt_fsp_nonconf_coll,
                       **run_opt_entropyCons_multicomponent,
                       #**run_opt_entropyStab,
                      }
         }
   #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
   # Split-form DG and shock capturing
   caseID=caseID+1
   jobs['mhd_split_glm_noncons_nopara_SC']={
          'case': caseID ,
          'tags': ['mhd','multicomponent','split-form','SC','GL','GLM','NONCONS'] ,
          'build_opts':{**baseopts,
                        'FLUXO_DISCTYPE'          :'2',
                        'FLUXO_EQN_NUM_COMPONENTS':'2',
                        'FLUXO_DISC_NODETYPE'     :'GAUSS-LOBATTO',
                        'FLUXO_EQN_GLM'           :'ON',
                        'FLUXO_EQN_NONCONS'       :'ON',
                        'FLUXO_EQN_NONCONS_GLM'   :'ON',
                        'FLUXO_PARABOLIC'         :'OFF',
                        'FLUXO_SHOCKCAPTURE'      :'ON',
                        'FLUXO_SHOCKCAP_NFVSE'    :'ON',
                        'FLUXO_SHOCKINDICATOR'    :'custom',
                       },
          'run_opts': {
                       **run_opt_fsp_conf_multicomponent,
                       #**run_opt_fsp_nonconf_coll,
                       #**run_opt_fsp_SC_first,
                       **run_opt_entropyCons_multicomponent,
                       #**run_opt_entropyStab,
                       #**run_opt_entropyStab_FloGor9waves,
                      },

         }
   #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
   # Standard DG
   caseID=caseID+1
   jobs['mhd_type1_br1_cons_GL']={
          'case': caseID ,
          'tags': ['mhd','multicomponent','standardDG','br1','GL','GLM'] ,
          'build_opts':{**baseopts,
                        'FLUXO_DISCTYPE'             :'1',
                        'FLUXO_EQN_NUM_COMPONENTS'   :'2',
                        'FLUXO_DISC_NODETYPE'        :'GAUSS-LOBATTO',
                        'FLUXO_EQN_GLM'              :'ON',
                        'FLUXO_PARABOLIC'            :'ON',
                        'FLUXO_PARABOLIC_LIFTING'    :'br1',
                        'FLUXO_PARABOLIC_LIFTING_VAR':'cons_var',
                       },
          'run_opts': {**run_opt_fsp_conf_multicomponent, 
                       #**run_opt_fsp_nonconf_proj,
                      }
         }
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
         {'tags': ['navierstokes','freestream','curved','conforming','restart'] ,
          'restartfile': 'NAVIERSTOKES_FREESTREAM_State_0000000.000000000.h5',
          'test_opts':{'err_Linf':{'func': check_error ,'f_kwargs': {'whichError':'L_inf ','err_tol': 1e-11} } ,
                      },
         },
      }
   run_opt_fsp_nonconf_proj={'runs/navst/freestream/nonconforming_projmortar':
         {'tags': ['navierstokes','freestream','curved','nonconforming','projection-mortar'] ,
          'test_opts':{'err_Linf':{'func': check_error ,'f_kwargs': {'whichError':'L_inf ','err_tol': 1e-11} } ,
                      },
         },
      }
   run_opt_fsp_nonconf_coll={'runs/navst/freestream/nonconforming_collmortar':
         {'tags': ['navierstokes','freestream','curved','nonconforming','collocation-mortar'] ,
          'test_opts':{'err_Linf':{'func': check_error ,'f_kwargs': {'whichError':'L_inf ','err_tol': 1e-11} } ,
                      },
         },
      }
   run_opt_fsp_nonconf_coll_central={'runs/navst/freestream/nonconforming_collmortar_central':
         {'tags': ['navierstokes','freestream','curved','nonconforming','collocation-mortar','centralflux'] ,
          'test_opts':{'err_Linf':{'func': check_error ,'f_kwargs': {'whichError':'L_inf ','err_tol': 1e-10} } ,
                      },
         },
      }
   # Entropy conservation test with EC-KEP flux (with and without shock-capturing)
   run_opt_entropyCons={'runs/navst/softBlast/entropyCons':
         {'tags': ['navierstokes','entropyCons','conforming','curved'] ,
          'test_opts':{'abs(dSdU*Ut)':{'func': check_all_errors ,
                                       'f_kwargs': {'whichError':'dSdU*Ut','err_tol': 1e-13,'err_abs':True} } ,
                      },
         },
      } 
   # Entropy stability test with EC-KEP flux and LLF (with and without shock-capturing)
   run_opt_entropyStab={'runs/navst/softBlast/entropyStab':
         {'tags': ['navierstokes','entropyCons','conforming','curved'] ,
          'test_opts':{'dSdU*Ut':{'func': check_all_errors ,
                                  'f_kwargs': {'whichError':'dSdU*Ut','err_tol': 1e-13,'err_abs':False} } ,
                      },
         },
      } 
   #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
   # TEST FOR JESSE_MORTAR: (CURVED+MORTAR ONLY EC with JESSE_MORTAR)
   # Entropy conservation test with EC-KEP flux (with and without shock-capturing)
   run_opt_entropyCons_nonconf={'runs/navst/softBlast/entropyCons_nonconforming':
         {'tags': ['navierstokes','entropyCons','nonconforming','curved'] ,
          'test_opts':{'abs(dSdU*Ut)':{'func': check_all_errors ,
                                       'f_kwargs': {'whichError':'dSdU*Ut','err_tol': 1e-13,'err_abs':True} } ,
                      },
         },
      } 
   # Entropy stability test with EC-KEP flux and LLF (with and without shock-capturing)
   run_opt_entropyStab_nonconf={'runs/navst/softBlast/entropyStab_nonconforming':
         {'tags': ['navierstokes','entropyCons','nonconforming','curved'] ,
          'test_opts':{'dSdU*Ut':{'func': check_all_errors ,
                                  'f_kwargs': {'whichError':'dSdU*Ut','err_tol': 1e-13,'err_abs':False} } ,
                      },
         },
      }
   # Entropy conservation test with Ranocha EC-KEP flux (with and without shock-capturing)
   run_opt_entropyCons_Ranocha={'runs/navst/softBlast/entropyCons_Ranocha':
         {'tags': ['navierstokes','entropyCons','conforming'] ,
          'test_opts':{'abs(dSdU*Ut)':{'func': check_all_errors ,
                                       'f_kwargs': {'whichError':'dSdU*Ut','err_tol': 1e-13,'err_abs':True} } ,
                      },
         },
      }
   
   # Entropy conservation test with Ranocha EC-KEP flux, SC and AMR
   run_opt_entropyCons_AMR={'runs/navst/softBlast/entropyCons_AMR':
         {'tags': ['navierstokes','entropyCons','nonconforming','curved','amr','p4est'] ,
          'test_opts':{'abs(dSdU*Ut)':{'func': check_all_errors ,
                                       'f_kwargs': {'whichError':'dSdU*Ut','err_tol': 1e-13,'err_abs':True} } ,
                      },
         },
      } 
   # Entropy stability test with Ranocha EC-KEP flux, full wave dissipation, reconstructed FV-SC and AMR
   run_opt_entropyStab_AMR={'runs/navst/softBlast/entropyStab_AMR':
         {'tags': ['navierstokes','entropyStab','nonconforming','curved','amr','p4est'] ,
          'test_opts':{'abs(dSdU*Ut)':{'func': check_all_errors ,
                                       'f_kwargs': {'whichError':'dSdU*Ut','err_tol': 1e-13,'err_abs':False} } ,
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
                       **run_opt_fsp_nonconf_proj,
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
                       **run_opt_fsp_nonconf_coll,
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
                       **run_opt_fsp_nonconf_proj,
                      }
         }
   #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
   caseID=caseID+1
   jobs['build_navierstokes_type2_br1_entropy_var_jesse_mortar']={
          'case': caseID,
          'tags': [ 'navierstokes','split-form','GL','br1','entropy_var','jesse-mortar'],
          'build_opts':{**baseopts,
                        'FLUXO_DISCTYPE'         :'2',
                        'FLUXO_DISC_NODETYPE'    :'GAUSS-LOBATTO',
                        'FLUXO_PARABOLIC'        :'ON',
                        'FLUXO_PARABOLIC_LIFTING':'br1',
                        'FLUXO_PARABOLIC_LIFTING_VAR':'entropy_var',
                        'FLUXO_JESSE_MORTAR'     :'ON',
                       },
          'run_opts': {**run_opt_fsp_conf, 
                       **run_opt_fsp_nonconf_coll,
                       **run_opt_entropyCons,
                       **run_opt_entropyStab,
                       **run_opt_entropyCons_nonconf,
                       **run_opt_entropyStab_nonconf,
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
                       **run_opt_fsp_nonconf_coll,
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
                       **run_opt_fsp_nonconf_coll,
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
                       **run_opt_fsp_nonconf_proj,
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
                       **run_opt_fsp_nonconf_proj,
                      }
         }
   #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
   volfluxes=['-1','0','5','8']
   for  vvv  in range(0,len(volfluxes)):
      volflux=volfluxes[vvv]
      caseID=caseID+1
      
      my_run_opts={**run_opt_fsp_conf }
      if vvv!=1:
        my_run_opts={**my_run_opts,
                     ** run_opt_fsp_nonconf_coll }
      if vvv==0 or vvv==2:
        my_run_opts={**my_run_opts,
                     **run_opt_entropyCons,
                     **run_opt_entropyStab,
                     **run_opt_entropyCons_Ranocha,
                    }
      
      jobs['build_navierstokes_type2_nopara_volFlux_'+volflux]={
            'case': caseID,
            'tags': [ 'navierstokes','split-form','GL'],
            'build_opts':{**baseopts,
                          'FLUXO_DISCTYPE'         :'2',
                          "FLUXO_DISC_NODETYPE"    :"GAUSS-LOBATTO",
                          "FLUXO_PARABOLIC"        :"OFF",
                          'FLUXO_EQN_VOLFLUX'      : volflux,
                         },
            'run_opts': my_run_opts
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
                       **run_opt_fsp_nonconf_proj,
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
                       **run_opt_fsp_nonconf_coll,
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
                       **run_opt_fsp_nonconf_coll,
                       **run_opt_shock_SC,
                       **run_opt_entropyCons,
                       **run_opt_entropyStab,
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
                       **run_opt_fsp_nonconf_coll,
                       **run_opt_fsp_SC_TVD,
                       **run_opt_fsp_SC_TVD_ES_fix,
                       **run_opt_fsp_SC_TVD_ES_Fjordholm,
                      }
         }
   #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
   # New test: Free stream for AMR simulations
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
                       **run_opt_fsp_nonconf_coll,
                       **run_opt_fsp_p4est,
                      }
         }
    #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
   caseID=caseID+1
   jobs['build_navierstokes_type2_nopara_Zhang-Shu_posit']={
          'case': caseID,
          'tags': [ 'navierstokes','split-form','GL','positivity','Zhang-Shu'],
          'build_opts':{**baseopts,
                        'FLUXO_DISCTYPE'         :'2',
                        "FLUXO_DISC_NODETYPE"    :"GAUSS-LOBATTO",
                        "FLUXO_PARABOLIC"        :"OFF",
                        "FLUXO_POSITIVITYPRES"   :"ON",
                       },
          'run_opts': {**run_opt_fsp_conf, 
                       **run_opt_fsp_nonconf_coll,
                       **run_opt_shock_posit,
                      }
         }
   #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
   # New test: Soft blast with SC, BR2, and non-conforming
   run_opt_nonConf_parabolic={'runs/navst/softBlast/nonconfParabolic':
         {'tags': ['navierstokes','blast','curved','nonconforming','projection-mortar','br2','SC'] ,
          'test_opts':{'max|Ut|':{'func': check_error ,
                                  'f_kwargs': {'whichError':'max|Ut| ',
                                                    'to_be': [5.945798911280E-01,   2.065011304993E-01,   2.003460892178E-01,   2.207665371548E-01,   2.168581765236E+00],
                                                  'err_tol': 1e-11} } , # err_tol is limited by the precision of the output..
                      },
         },
      }
   #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
   caseID=caseID+1
   jobs['build_navierstokes_type2_SC_br2_entr_vars']={
          'case': caseID,
          'tags': [ 'navierstokes','split-form','GL','br2','entropy_var','nonconforming','SC'],
          'build_opts':{**baseopts,
                        'FLUXO_DISCTYPE'         :'2',
                        "FLUXO_DISC_NODETYPE"    :"GAUSS-LOBATTO",
                        "FLUXO_PARABOLIC"        :"ON",
                        "FLUXO_PARABOLIC_LIFTING":"br2",
                        "FLUXO_PARABOLIC_LIFTING_VAR":"entropy_var",
                        "FLUXO_SHOCKCAPTURE"     :"ON",
                        "FLUXO_SHOCKCAP_NFVSE"   :"ON",
                       },
          'run_opts': {
                       **run_opt_nonConf_parabolic,
                      }
         }
   #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
   caseID=caseID+1
   jobs['build_navierstokes_p4est_SC_jesse_mortar']={
          'case': caseID,
          'tags': [ 'navierstokes','split-form','GL','amr','p4est','SC','jesse-mortar'],
          'build_opts':{**baseopts,
                        'FLUXO_DISCTYPE'         :'2',
                        "FLUXO_DISC_NODETYPE"    :"GAUSS-LOBATTO",
                        "FLUXO_PARABOLIC"        :"OFF",
                        "FLUXO_AMR"              :"ON",
                        "FLUXO_BUILD_P4EST"      :"OFF",
                        'FLUXO_JESSE_MORTAR'     :'ON',
                        "FLUXO_SHOCKCAPTURE"     :"ON",
                        "FLUXO_SHOCKCAP_NFVSE"   :"ON",
                       },
          'run_opts': {**run_opt_fsp_conf, 
                       **run_opt_fsp_nonconf_coll,
                       **run_opt_fsp_p4est,
                       **run_opt_entropyCons_AMR,
                       **run_opt_entropyStab_AMR,
                      }
         }
   #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
   # ES Gauss tests
   #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
   volfluxes=['-1','5']
   for  vvv  in range(0,len(volfluxes)):
      volflux=volfluxes[vvv]
      caseID=caseID+1
      
      jobs['build_navierstokes_esgauss_nopara_volFlux_'+volflux]={
            'case': caseID,
            'tags': [ 'navierstokes','split-form','gauss','esgauss'],
            'build_opts':{**baseopts,
                          'FLUXO_DISCTYPE'         :'2',
                          "FLUXO_DISC_NODETYPE"    :"GAUSS",
                          "FLUXO_PARABOLIC"        :"OFF",
                          'FLUXO_EQN_VOLFLUX'      : volflux,
                         },
            'run_opts': {**run_opt_fsp_conf,
                         **run_opt_entropyCons,
                         **run_opt_entropyStab,
                         **run_opt_entropyCons_Ranocha,
                        }
           }
   #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
   caseID=caseID+1
   jobs['build_navierstokes_p4est_esgauss']={
          'case': caseID,
          'tags': [ 'navierstokes','split-form','gauss','esgauss','amr','p4est','jesse-mortar'],
          'build_opts':{**baseopts,
                        'FLUXO_DISCTYPE'         :'2',
                        "FLUXO_DISC_NODETYPE"    :"GAUSS",
                        "FLUXO_PARABOLIC"        :"OFF",
                        "FLUXO_AMR"              :"ON",
                        "FLUXO_BUILD_P4EST"      :"OFF",
                        'FLUXO_JESSE_MORTAR'     :'ON',
                       },
          'run_opts': {**run_opt_fsp_conf, 
                       **run_opt_fsp_nonconf_coll,
                       # #**run_opt_fsp_p4est,     # TODO: This is failing with L_inf = 3.417000016270E-11 > 1.e-11 (is this problematic?)
                       **run_opt_entropyCons_AMR,
                       **run_opt_entropyStab_AMR,
                      }
         }
   #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
   
   #============================================================================
   #============================================================================
   #PERFORMANCE, MHD, 3000 < caseID
   #============================================================================
   #============================================================================
   caseID=3000
   baseopts={ 'FLUXO_EQNSYSNAME' :'mhd',
              'FLUXO_EQN_GLM'    :'ON',
            }
   
   perf_run_opt_ec_conf={'runs/mhd/performance/EC_conforming':
         {'tags': ['mhd','curved','conforming'] ,
         },
      }
   perf_run_opt_ec_nonconf_2to1={'runs/mhd/performance/EC_nonconforming_2to1':
         {'tags': ['mhd','curved','nonconforming','2to1-mortar'] ,
         },
      }
   perf_run_opt_ec_nonconf_4to1={'runs/mhd/performance/EC_nonconforming_4to1':
         {'tags': ['mhd','curved','nonconforming','4to1-mortar'] ,
         },
      }
   #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
   caseID=caseID+1
   jobs['build_mhd_type2_nopara_no-noncons_EC_performance']={
          'case': caseID,
          'tags': [ 'mhd','split-form','GL','performance'],
          'build_opts':{**baseopts,
                        'FLUXO_DISCTYPE'         :'2',
                        "FLUXO_DISC_NODETYPE"    :"GAUSS-LOBATTO",
                        "FLUXO_PARABOLIC"        :"OFF",
                        'FLUXO_EQN_NONCONS'      :'OFF',
                       },
          'run_opts': {**perf_run_opt_ec_conf, 
                       **perf_run_opt_ec_nonconf_2to1, 
                       **perf_run_opt_ec_nonconf_4to1, 
                      }
         }
   #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
   caseID=caseID+1
   jobs['build_mhd_type2_jesse_nopara_no-noncons_EC_performance']={
          'case': caseID,
          'tags': [ 'mhd','split-form','GL','performance','jesse-mortar'],
          'build_opts':{**baseopts,
                        'FLUXO_DISCTYPE'         :'2',
                        "FLUXO_DISC_NODETYPE"    :"GAUSS-LOBATTO",
                        "FLUXO_PARABOLIC"        :"OFF",
                        'FLUXO_EQN_NONCONS'      :'OFF',
                        'FLUXO_JESSE_MORTAR'     :'ON',
                       },
          'run_opts': {**perf_run_opt_ec_conf, 
                       **perf_run_opt_ec_nonconf_2to1, 
                       **perf_run_opt_ec_nonconf_4to1, 
                      }
         }
   #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
   caseID=caseID+1
   jobs['build_mhd_type2_nopara_noncons_EC_performance']={
          'case': caseID,
          'tags': [ 'mhd','split-form','GL','performance'],
          'build_opts':{**baseopts,
                        'FLUXO_DISCTYPE'         :'2',
                        "FLUXO_DISC_NODETYPE"    :"GAUSS-LOBATTO",
                        "FLUXO_PARABOLIC"        :"OFF",
                        'FLUXO_EQN_NONCONS'      :'ON',
                       },
          'run_opts': {**perf_run_opt_ec_conf, 
                       **perf_run_opt_ec_nonconf_2to1, 
                       **perf_run_opt_ec_nonconf_4to1, 
                      }
         }
   #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
   caseID=caseID+1
   jobs['build_mhd_type2_jesse_nopara_noncons_EC_performance']={
          'case': caseID,
          'tags': [ 'mhd','split-form','GL','performance','jesse-mortar'],
          'build_opts':{**baseopts,
                        'FLUXO_DISCTYPE'         :'2',
                        "FLUXO_DISC_NODETYPE"    :"GAUSS-LOBATTO",
                        "FLUXO_PARABOLIC"        :"OFF",
                        'FLUXO_EQN_NONCONS'      :'ON',
                        'FLUXO_JESSE_MORTAR'     :'ON',
                       },
          'run_opts': {**perf_run_opt_ec_conf, 
                       **perf_run_opt_ec_nonconf_2to1, 
                       **perf_run_opt_ec_nonconf_4to1, 
                      }
         }
   #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
   caseID=caseID+1
   jobs['build_mhd_type2_noncons_br1_entropy_vars_EC_performance']={
          'case': caseID,
          'tags': [ 'mhd','split-form','GL','performance','br1'],
          'build_opts':{**baseopts,
                        'FLUXO_DISCTYPE'         :'2',
                        "FLUXO_DISC_NODETYPE"    :"GAUSS-LOBATTO",
                        'FLUXO_EQN_NONCONS'      :'ON',
                        "FLUXO_PARABOLIC"        :"ON",
                        "FLUXO_PARABOLIC_LIFTING":"br1",
                        "FLUXO_PARABOLIC_LIFTING_VAR":"entropy_var",
                       },
          'run_opts': {**perf_run_opt_ec_conf, 
                       **perf_run_opt_ec_nonconf_2to1, 
                       **perf_run_opt_ec_nonconf_4to1, 
                      }
         }
   #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
   caseID=caseID+1
   jobs['build_mhd_type2_noncons_jesse_br1_entropy_vars_EC_performance']={
          'case': caseID,
          'tags': [ 'mhd','split-form','GL','performance','br1','jesse-mortar'],
          'build_opts':{**baseopts,
                        'FLUXO_DISCTYPE'         :'2',
                        "FLUXO_DISC_NODETYPE"    :"GAUSS-LOBATTO",
                        'FLUXO_EQN_NONCONS'      :'ON',
                        "FLUXO_PARABOLIC"        :"ON",
                        "FLUXO_PARABOLIC_LIFTING":"br1",
                        "FLUXO_PARABOLIC_LIFTING_VAR":"entropy_var",
                        'FLUXO_JESSE_MORTAR'     :'ON',
                       },
          'run_opts': {**perf_run_opt_ec_conf, 
                       **perf_run_opt_ec_nonconf_2to1, 
                       **perf_run_opt_ec_nonconf_4to1, 
                      }
         }
   #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
   #============================================================================
   #============================================================================
   #PERFORMANCE, Navierstokes, 4000 < caseID
   #============================================================================
   #============================================================================
   caseID=4000
   baseopts={ 'FLUXO_EQNSYSNAME' :'navierstokes',
              'FLUXO_TESTCASE'   :'default'
            }
   
   perf_run_opt_ec_conf={'runs/navst/performance/EC_conforming':
         {'tags': ['navierstokes','curved','conforming'] ,
         },
      }
   perf_run_opt_ec_nonconf_2to1={'runs/navst/performance/EC_nonconforming_2to1':
         {'tags': ['navierstokes','curved','nonconforming','2to1-mortar'] ,
         },
      }
   perf_run_opt_ec_nonconf_4to1={'runs/navst/performance/EC_nonconforming_4to1':
         {'tags': ['navierstokes','curved','nonconforming','4to1-mortar'] ,
         },
      }
   #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
   caseID=caseID+1
   jobs['build_navierstokes_type2_nopara_EC_performance']={
          'case': caseID,
          'tags': [ 'navierstokes','split-form','GL','performance'],
          'build_opts':{**baseopts,
                        'FLUXO_DISCTYPE'         :'2',
                        "FLUXO_DISC_NODETYPE"    :"GAUSS-LOBATTO",
                        "FLUXO_PARABOLIC"        :"OFF",
                       },
          'run_opts': {**perf_run_opt_ec_conf, 
                       **perf_run_opt_ec_nonconf_2to1, 
                       **perf_run_opt_ec_nonconf_4to1, 
                      }
         }
   #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
   caseID=caseID+1
   jobs['build_navierstokes_type2_jesse_nopara_EC_performance']={
          'case': caseID,
          'tags': [ 'navierstokes','split-form','GL','performance','jesse-mortar'],
          'build_opts':{**baseopts,
                        'FLUXO_DISCTYPE'         :'2',
                        "FLUXO_DISC_NODETYPE"    :"GAUSS-LOBATTO",
                        "FLUXO_PARABOLIC"        :"OFF",
                        'FLUXO_JESSE_MORTAR'     :'ON',
                       },
          'run_opts': {**perf_run_opt_ec_conf, 
                       **perf_run_opt_ec_nonconf_2to1, 
                       **perf_run_opt_ec_nonconf_4to1, 
                      }
         }
   #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
   caseID=caseID+1
   jobs['build_navierstokes_type2_br1_entropy_vars_EC_performance']={
          'case': caseID,
          'tags': [ 'navierstokes','split-form','GL','performance','br1'],
          'build_opts':{**baseopts,
                        'FLUXO_DISCTYPE'         :'2',
                        "FLUXO_DISC_NODETYPE"    :"GAUSS-LOBATTO",
                        "FLUXO_PARABOLIC"        :"ON",
                        "FLUXO_PARABOLIC_LIFTING":"br1",
                        "FLUXO_PARABOLIC_LIFTING_VAR":"entropy_var",
                       },
          'run_opts': {**perf_run_opt_ec_conf, 
                       **perf_run_opt_ec_nonconf_2to1, 
                       **perf_run_opt_ec_nonconf_4to1, 
                      }
         }
   #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
   caseID=caseID+1
   jobs['build_navierstokes_type2_jesse_br1_entropy_vars_EC_performance']={
          'case': caseID,
          'tags': [ 'navierstokes','split-form','GL','performance','br1','jesse-mortar'],
          'build_opts':{**baseopts,
                        'FLUXO_DISCTYPE'         :'2',
                        "FLUXO_DISC_NODETYPE"    :"GAUSS-LOBATTO",
                        "FLUXO_PARABOLIC"        :"ON",
                        "FLUXO_PARABOLIC_LIFTING":"br1",
                        "FLUXO_PARABOLIC_LIFTING_VAR":"entropy_var",
                        'FLUXO_JESSE_MORTAR'     :'ON',
                       },
          'run_opts': {**perf_run_opt_ec_conf, 
                       **perf_run_opt_ec_nonconf_2to1, 
                       **perf_run_opt_ec_nonconf_4to1, 
                      }
         }
   #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~






   return jobs
