#!/usr/bin/python
# -*- coding: utf-8 -*-

import os
import sys
import math
import shutil
import tempfile
import time
import argparse

# use helpers.py
from helpers import execute, read_prm
from helpers import get_last_L2_error, get_last_Linf_error, get_cpu_per_dof

########################################################################################################

def test_fluxo( buildopts=None , case=0, project="test" , ntail = 0 , 
                run_test=None , mpi_procs = 1, keepdir=0, 
                err= None ):
   if ( len(buildopts) < 2 ) :
      print "error, nobuild options given"
      err.extend(["caseID=%6d ,project= %s" % (case,project)])
      return False
   # build directory
   builddir=("dirx_%d_%s" % (case,project))
   log_path=("../log.%d_%s" % (case,project))

   os.system('rm -rf '+builddir) #destroy directory from previous run
   os.system('mkdir '+builddir )
   cwd = os.getcwd()  #current working directory
   os.chdir(builddir)

   allopts=" "
   for i in range((len(buildopts)/2)) :
      allopts=allopts+" -D"+buildopts[2*i]+"="+buildopts[2*i+1]+" " 

   print "===> OPTIONS:"
   print allopts

   logf = open(log_path, 'w')
   logf.write("===> OPTIONS: \n %s \n" % (allopts) )
 
   cmdconfig = "cmake ../../."+allopts
   print "  "
   print "===> configure..."

   os.system(cmdconfig+" 2>stderr 1>stdout")
   stdout=open("stdout",'r').readlines()
   stderr=open("stderr",'r').readlines()

   success= True
   for line in stderr :
      if "ERROR" in line :
         success= False
         sys.stdout.write("%s" % line)

   logf.write('CONFIG: STDOUT FILE:\n')
   for line in stdout :
      logf.write(line)

   logf.write('CONFIG: STDERR FILE:\n')
   for line in stderr :
      logf.write(line)

   logf.close()
   if (ntail > 0) :
      nlines=len(stdout)
      i=0
      for line in stdout :
         i=i+1
         if ( i > (nlines -ntail) ) :
            sys.stdout.write("%s" % line)
   if(not success) :
      print "===============================================  "
      print " !!!! ERROR IN CMAKE, NOT FINISHED CORRECTLY!!!!! "
      print "===============================================  "
      print "  "
      err.extend(["caseID=%6d ,project= %s <=Cmake" % (case,project)])
      os.chdir(cwd)
      return success #=False
   else : 
      print "===============================================  "
      print " Cmake finished successfully."
      print "===============================================  "
      print "  "
   #MAKE
   print "  "
   print "===> make..."
   cmdmake = ("make -j %d VERBOSE=1" % (mpi_procs) )

   os.system(cmdmake+" 2>stderr 1>stdout")
   stdout=open("stdout",'r').readlines()
   stderr=open("stderr",'r').readlines()
   
   logf = open(log_path, 'a') #append

   success= False
   logf.write('\n ====> MAKE: STDOUT FILE:\n \n')
   for line in stdout :
      logf.write(line)
      if "SUCCESS: FLUXO BUILD COMPLETE!" in line :
         success= True

   logf.write('\n ====> MAKE: STDERR FILE:\n \n')
   for line in stderr :
      logf.write(line)

   logf.close()
   if (ntail > 0) :
      nlines=len(stdout)
      i=0
      for line in stdout :
         i=i+1
         if ( i > (nlines -ntail) ) :
            sys.stdout.write("%s" % line)

   if (not success) :
      print "===============================================  "
      print " !!!! PROBLEM WITH BUILD, NOT FINISHED CORRECTLY!!!!! "
      err.extend(["caseID=%6d ,project= %s <=Make" % (case,project)])
      print "===============================================  "
      print "  "
      os.chdir(cwd)
      return success #=False
   else :
      print "===============================================  "
      print " Build finished sucessfully."
      print "===============================================  "
      print "  "


   os.chdir(cwd)
   doruntest=False
   if( run_test ) :
      if( len(run_test) == 3 ) :  # run tests
         doruntest = True
   if(doruntest) :
      #RUN TEST
      print "  "
      print "===> run test %s/%s ..." % (run_test[0],run_test[1])
      #change to test directory
      os.chdir(run_test[0])
      if (not os.path.isfile('./'+run_test[1])) :
         print "!!!!!!!!  PARAMTERFILE %s/%s does not exist!" %(run_test[0],run_test[1])
         os.chdir(cwd)
         return False

      projectname = read_prm(run_test[1],'ProjectName')
      projectnamex=("%s_%d_%s" % (projectname,case,project))
      success = False
      try :
         [L2,Linf,PID] = execute(cwd+"/"+builddir+"/bin/fluxo", run_test[1], projectnamex,\
                               [get_last_L2_error, get_last_Linf_error, get_cpu_per_dof],\
                               log = True, ntail = ntail ,\
                                     mpi_procs = mpi_procs )
         if(Linf) :
            print "   ... check Linf %s < %s ?" % (Linf[0],run_test[2])
            if(float(Linf[0]) < float(run_test[2])) :
               success = True
      except :
         success = False
      if (not success) :
         print "===============================================  "
         print " !!!! PROBLEM WITH RUN, NOT FINISHED CORRECTLY!!!!! "
         err.extend(["caseID=%6d ,project= %s <=Run" % (case,project)])
         print "===============================================  "
         print "  "
         os.chdir(cwd)
         return success #=False
      else :
         print "===============================================  "
         print " Run finished sucessfully."
         print "===============================================  "
         print "  "
   else :
      print " (no test specified...) "
        

   # success still True, only delete sucessfull build if keepdir=0
   os.chdir(cwd)
   if(keepdir ==0) :
      os.system('rm -rf '+builddir) #destroy directory
   return success

     
     
  
########################################################################################################
# parse a str with comma separated ranges: 1,5,10-12
def parse_range(astr):
    result = set()
    for part in astr.split(','):
       x = part.split('-')
       result.update(range(int(x[0]), int(x[-1]) + 1))
    return sorted(result)

########################################################################################################

# MAIN PROGRAM

########################################################################################################
parser = argparse.ArgumentParser(description='Tool to build/compile fluxo in different configurations',\
                                 formatter_class=argparse.RawTextHelpFormatter)
parser.add_argument('-p','--procs',type=int,default=1, help='    number of processors used for the make (DEFAULT=1)')

parser.add_argument('-withmpi', type=int, default=1,   help="1 : DEFAULT ,compile with mpi\n"
                                                            "0 : compile without mpi  ")

parser.add_argument('-buildhdf5', type=int, default=1, help="1 : DEFAULT, build hdf5 locally,\n"
                                                            "0 : use external hdf5 (modules)")

parser.add_argument('-hostname', type=str, default="", help="    cmake hostname, only needed if compiling on a cluster" )

parser.add_argument('-case', type=str, default='0',    help="0 : DEFAULT, run all cases,\n" 
                                                            "101,102-104 : list of specific cases to run (without spaces!) ")

parser.add_argument('-runtests',type=int, default='0', help="1 : make a test run with the executable and a given parameterfile,\n" 
                                                            "0 : DEFAULT, do not run any tests. ")

parser.add_argument('-ntail', type=int, default=5,     help='    number of last line output of cmake/make screenlog (DEFAULT=5)')

parser.add_argument('-keepdir', type=int, default=1,   help="1 : DEFAULT, keep all build directories,\n" 
                                                            "0 : delete sucessfull build directories")


args = parser.parse_args()

cases = parse_range(args.case)

if(args.withmpi == 0) :
   MPIOPT="OFF"
else :
   MPIOPT="ON"
globopts=["FLUXO_BUILD_MPI"        ,MPIOPT 
         ]


if(args.buildhdf5 == 0) :
   HDF5OPT="OFF"
else :
   HDF5OPT="ON"

globopts.extend([
          "FLUXO_BUILD_HDF5"       ,HDF5OPT
         ])

if(len(args.hostname) > 1 ) :
   globopts.extend([
                    "CMAKE_HOSTNAME"       ,args.hostname 
                   ])

builderr= []

dbg  = False #debug, no builds
stat = True

#============================================================================
#============================================================================
#first group, LINADV, 100 < caseID <200
#============================================================================
#============================================================================
caseID=100 
baseopts=[
           "FLUXO_EQNSYSNAME"       ,"linearscalaradvection"
          ,"FLUXO_TESTCASE"         ,"default"
         ]

TEST=[]
if(args.runtests == 1) :
   # relative path from tests folder, parameterfile,Linf[0]<crit for success
   TEST.extend(["freestream","parameter_freestream_linadv.ini", "1.0e-10" ])
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
caseID=caseID+1
if(cases[0] ==0 or (caseID in cases)) :
   pname="build_linadv_release_type1_GL"
   print "caseID: %d name: %s" % (caseID,pname)

   options=[]; options.extend(globopts) ; options.extend(baseopts)
   options.extend([
             "CMAKE_BUILD_TYPE"       ,"Release"
            ,"FLUXO_DISCTYPE"         ,"1"
            ,"FLUXO_DISC_NODETYPE"    ,"GAUSS-LOBATTO"
            ,"FLUXO_PARABOLIC"        ,"ON"
            ,"FLUXO_PARABOLIC_LIFTING","br1"
           ])
   
   if(not dbg ) : stat = test_fluxo(buildopts=options, case=caseID, project=pname,ntail = args.ntail,\
                          run_test=TEST , mpi_procs = args.procs , keepdir=args.keepdir, err=builderr )
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
caseID=caseID+1
if(cases[0] ==0 or (caseID in cases)) :
   pname="build_linadv_release_type2_GL"
   print "caseID: %d name: %s" % (caseID,pname)

   options=[]; options.extend(globopts) ; options.extend(baseopts)
   options.extend([
             "CMAKE_BUILD_TYPE"       ,"Release"
            ,"FLUXO_DISCTYPE"         ,"2"
            ,"FLUXO_DISC_NODETYPE"    ,"GAUSS-LOBATTO"
            ,"FLUXO_PARABOLIC"        ,"ON"
            ,"FLUXO_PARABOLIC_LIFTING","br1"
           ])
   
   if(not dbg ) : stat = test_fluxo(buildopts=options, case=caseID, project=pname, ntail = args.ntail ,\
                          run_test=TEST , mpi_procs = args.procs , keepdir=args.keepdir, err=builderr )
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
caseID=caseID+1
if(cases[0]==0 or (caseID in cases)) :
   pname="build_linadv_type2_br1"
   print "caseID: %d name: %s" % (caseID,pname)

   options=[]; options.extend(globopts) ; options.extend(baseopts)
   options.extend([
             "CMAKE_BUILD_TYPE"       ,"Debug"
            ,"FLUXO_DISCTYPE"         ,"2"
            ,"FLUXO_DISC_NODETYPE"    ,"GAUSS-LOBATTO"
            ,"FLUXO_PARABOLIC"        ,"ON"
            ,"FLUXO_PARABOLIC_LIFTING","br1"
           ])
   
   if(not dbg ) : stat = test_fluxo(buildopts=options, case=caseID, project=pname, ntail = args.ntail ,\
                          run_test=TEST , mpi_procs = args.procs , keepdir=args.keepdir, err=builderr )
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
caseID=caseID+1
if(cases[0]==0 or (caseID in cases)) :
   pname="build_linadv_type1_GL_br1"
   print "caseID: %d name: %s" % (caseID,pname)

   options=[]; options.extend(globopts) ; options.extend(baseopts)
   options.extend([
             "CMAKE_BUILD_TYPE"       ,"Debug"
            ,"FLUXO_DISCTYPE"         ,"1"
            ,"FLUXO_DISC_NODETYPE"    ,"GAUSS-LOBATTO"
            ,"FLUXO_PARABOLIC"        ,"ON"
            ,"FLUXO_PARABOLIC_LIFTING","br1"
           ])
   
   if(not dbg ) : stat = test_fluxo(buildopts=options, case=caseID, project=pname, ntail = args.ntail ,\
                          run_test=TEST , mpi_procs = args.procs , keepdir=args.keepdir, err=builderr )
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
caseID=caseID+1
if(cases[0]==0 or (caseID in cases)) :
   pname="build_linadv_type1_GL_br2"
   print "caseID: %d name: %s" % (caseID,pname)

   options=[]; options.extend(globopts) ; options.extend(baseopts)
   options.extend([
             "CMAKE_BUILD_TYPE"       ,"Debug"
            ,"FLUXO_DISCTYPE"         ,"1"
            ,"FLUXO_DISC_NODETYPE"    ,"GAUSS-LOBATTO"
            ,"FLUXO_PARABOLIC"        ,"ON"
            ,"FLUXO_PARABOLIC_LIFTING","br2"
           ])
   
   if(not dbg ) : stat = test_fluxo(buildopts=options, case=caseID, project=pname, ntail = args.ntail ,\
                           run_test=TEST , mpi_procs = args.procs , keepdir=args.keepdir, err=builderr )
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
caseID=caseID+1
if(cases[0]==0 or (caseID in cases)) :
   pname="build_linadv_type1_Gauss_br1"
   print "caseID: %d name: %s" % (caseID,pname)

   options=[]; options.extend(globopts) ; options.extend(baseopts)
   options.extend([
             "CMAKE_BUILD_TYPE"       ,"Debug"
            ,"FLUXO_DISCTYPE"         ,"1"
            ,"FLUXO_DISC_NODETYPE"    ,"GAUSS"
            ,"FLUXO_PARABOLIC"        ,"ON"
            ,"FLUXO_PARABOLIC_LIFTING","br1"
           ])
   
   if(not dbg ) : stat = test_fluxo(buildopts=options, case=caseID, project=pname, ntail = args.ntail ,\
                           run_test=TEST , mpi_procs = args.procs , keepdir=args.keepdir, err=builderr )
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
caseID=caseID+1
if(cases[0]==0 or (caseID in cases)) :
   pname="build_linadv_type1_Gauss_br2"
   print "caseID: %d name: %s" % (caseID,pname)

   options=[]; options.extend(globopts) ; options.extend(baseopts)
   options.extend([
             "CMAKE_BUILD_TYPE"       ,"Debug"
            ,"FLUXO_DISCTYPE"         ,"1"
            ,"FLUXO_DISC_NODETYPE"    ,"GAUSS"
            ,"FLUXO_PARABOLIC"        ,"ON"
            ,"FLUXO_PARABOLIC_LIFTING","br2"
           ])
   
   if(not dbg ) : stat = test_fluxo(buildopts=options, case=caseID, project=pname, ntail = args.ntail ,\
                          run_test=TEST , mpi_procs = args.procs , keepdir=args.keepdir, err=builderr )
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
caseID=caseID+1
if(cases[0]==0 or (caseID in cases)) :
   pname="build_linadv_type2_nopara_cart"
   print "caseID: %d name: %s" % (caseID,pname)

   options=[]; options.extend(globopts) ; options.extend(baseopts)
   options.extend([
             "CMAKE_BUILD_TYPE"       ,"Debug"
            ,"FLUXO_DISCTYPE"         ,"2"
            ,"FLUXO_DISC_CARTESIANFLUX","ON"           #=> runtest off
            ,"FLUXO_DISC_NODETYPE"    ,"GAUSS-LOBATTO"
            ,"FLUXO_PARABOLIC"        ,"OFF"
           ])
   
   if(not dbg ) : stat = test_fluxo(buildopts=options, case=caseID, project=pname, ntail = args.ntail ,\
                          run_test=None , mpi_procs = args.procs , keepdir=args.keepdir, err=builderr )
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~


#============================================================================
#============================================================================
#next group, MAXWELL, 200 < caseID <300
#============================================================================
#============================================================================
caseID=200 
baseopts=[
           "FLUXO_EQNSYSNAME"       ,"maxwell"
          ,"FLUXO_TESTCASE"         ,"default"
         ]
TEST=[]
#if(args.runtests == 1) :
#   # relative path from tests folder, parameterfile,Linf[0]<crit for success
#   TEST.extend(["freestream","parameter_freestream_maxwell.ini", "1.0e-10" ])
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
caseID=caseID+1
if(cases[0] ==0 or (caseID in cases)) :
   pname="build_maxwell_type1_GL"
   print "caseID: %d name: %s" % (caseID,pname)

   options=[]; options.extend(globopts) ; options.extend(baseopts)
   options.extend([
             "CMAKE_BUILD_TYPE"       ,"Debug"
            ,"FLUXO_DISCTYPE"         ,"1"
            ,"FLUXO_DISC_NODETYPE"    ,"GAUSS-LOBATTO"
            ,"FLUXO_PARABOLIC"        ,"OFF"
           ])
   
   if(not dbg ) : stat = test_fluxo(buildopts=options, case=caseID, project=pname, ntail = args.ntail ,\
                          run_test=TEST , mpi_procs = args.procs , keepdir=args.keepdir, err=builderr )
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
caseID=caseID+1
if(cases[0] ==0 or (caseID in cases)) :
   pname="build_maxwell_type2_GL"
   print "caseID: %d name: %s" % (caseID,pname)

   options=[]; options.extend(globopts) ; options.extend(baseopts)
   options.extend([
             "CMAKE_BUILD_TYPE"       ,"Debug"
            ,"FLUXO_DISCTYPE"         ,"2"
            ,"FLUXO_DISC_NODETYPE"    ,"GAUSS-LOBATTO"
            ,"FLUXO_PARABOLIC"        ,"OFF"
           ])
   
   if(not dbg ) : stat = test_fluxo(buildopts=options, case=caseID, project=pname, ntail = args.ntail ,\
                          run_test=TEST , mpi_procs = args.procs , keepdir=args.keepdir, err=builderr )
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
caseID=caseID+1
if(cases[0] ==0 or (caseID in cases)) :
   pname="build_maxwell_type2_GL_CART"
   print "caseID: %d name: %s" % (caseID,pname)

   options=[]; options.extend(globopts) ; options.extend(baseopts)
   options.extend([
             "CMAKE_BUILD_TYPE"       ,"Debug"
            ,"FLUXO_DISCTYPE"         ,"2"
            ,"FLUXO_DISC_NODETYPE"    ,"GAUSS-LOBATTO"
            ,"FLUXO_DISC_CARTESIANFLUX","ON"           #=> runtest off
            ,"FLUXO_PARABOLIC"        ,"OFF"
           ])
   
   if(not dbg ) : stat = test_fluxo(buildopts=options, case=caseID, project=pname, ntail = args.ntail ,\
                          run_test=None , mpi_procs = args.procs , keepdir=args.keepdir, err=builderr )
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~



#============================================================================
#============================================================================
#third group, MHD, 300 < caseID <400
#============================================================================
#============================================================================
caseID=300 
baseopts=[
           "FLUXO_EQNSYSNAME"       ,"mhd"
         ]
TEST=[]
if(args.runtests == 1) :
   # relative path from tests folder, parameterfile,Linf[0]<crit for success
   TEST.extend(["freestream","parameter_freestream_mhd.ini", "1.0e-10" ])
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
caseID=caseID+1
if(cases[0] ==0 or (caseID in cases)) :
   pname="build_mhd_release_type1_GL"
   print "caseID: %d name: %s" % (caseID,pname)

   options=[]; options.extend(globopts) ; options.extend(baseopts)
   options.extend([
             "CMAKE_BUILD_TYPE"       ,"Release"
            ,"FLUXO_DISCTYPE"         ,"1"
            ,"FLUXO_EQN_GLM"          ,"ON"
            ,"FLUXO_DISC_NODETYPE"    ,"GAUSS-LOBATTO"
            ,"FLUXO_PARABOLIC"        ,"ON"
            ,"FLUXO_PARABOLIC_LIFTING","br1"
            ,"FLUXO_TESTCASE"         ,"default"
           ])
   
   if(not dbg ) : stat = test_fluxo(buildopts=options, case=caseID, project=pname, ntail = args.ntail ,\
                          run_test=TEST , mpi_procs = args.procs , keepdir=args.keepdir, err=builderr )
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
caseID=caseID+1
if(cases[0] ==0 or (caseID in cases)) :
   pname="build_mhd_release_type2_GL"
   print "caseID: %d name: %s" % (caseID,pname)

   options=[]; options.extend(globopts) ; options.extend(baseopts)
   options.extend([
             "CMAKE_BUILD_TYPE"       ,"Release"
            ,"FLUXO_DISCTYPE"         ,"2"
            ,"FLUXO_EQN_GLM"          ,"ON"
            ,"FLUXO_DISC_NODETYPE"    ,"GAUSS-LOBATTO"
            ,"FLUXO_PARABOLIC"        ,"ON"
            ,"FLUXO_PARABOLIC_LIFTING","br2"
            ,"FLUXO_TESTCASE"         ,"default"
           ])
   
   if(not dbg ) : stat = test_fluxo(buildopts=options, case=caseID, project=pname, ntail = args.ntail ,\
                          run_test=TEST , mpi_procs = args.procs , keepdir=args.keepdir, err=builderr )
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
caseID=caseID+1
if(cases[0]==0 or (caseID in cases)) :
   pname="build_mhd_type2_br1"
   print "caseID: %d name: %s" % (caseID,pname)

   options=[]; options.extend(globopts) ; options.extend(baseopts)
   options.extend([
             "CMAKE_BUILD_TYPE"       ,"Debug"
            ,"FLUXO_DISCTYPE"         ,"2"
            ,"FLUXO_EQN_GLM"          ,"ON"
            ,"FLUXO_DISC_NODETYPE"    ,"GAUSS-LOBATTO"
            ,"FLUXO_PARABOLIC"        ,"ON"
            ,"FLUXO_PARABOLIC_LIFTING","br1"
            ,"FLUXO_TESTCASE"         ,"default"
           ])
   
   if(not dbg ) : stat = test_fluxo(buildopts=options, case=caseID, project=pname, ntail = args.ntail ,\
                          run_test=TEST , mpi_procs = args.procs , keepdir=args.keepdir, err=builderr )
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
caseID=caseID+1
if(cases[0]==0 or (caseID in cases)) :
   pname="build_mhd_type1_GL_br1"
   print "caseID: %d name: %s" % (caseID,pname)

   options=[]; options.extend(globopts) ; options.extend(baseopts)
   options.extend([
             "CMAKE_BUILD_TYPE"       ,"Debug"
            ,"FLUXO_DISCTYPE"         ,"1"
            ,"FLUXO_EQN_GLM"          ,"ON"
            ,"FLUXO_DISC_NODETYPE"    ,"GAUSS-LOBATTO"
            ,"FLUXO_PARABOLIC"        ,"ON"
            ,"FLUXO_PARABOLIC_LIFTING","br1"
            ,"FLUXO_TESTCASE"         ,"default"
           ])
   
   if(not dbg ) : stat = test_fluxo(buildopts=options, case=caseID, project=pname, ntail = args.ntail ,\
                          run_test=TEST , mpi_procs = args.procs , keepdir=args.keepdir, err=builderr )
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
caseID=caseID+1
if(cases[0]==0 or (caseID in cases)) :
   pname="build_mhd_type1_GL_br2_TCmhd"
   print "caseID: %d name: %s" % (caseID,pname)

   options=[]; options.extend(globopts) ; options.extend(baseopts)
   options.extend([
             "CMAKE_BUILD_TYPE"       ,"Debug"
            ,"FLUXO_DISCTYPE"         ,"1"
            ,"FLUXO_EQN_GLM"          ,"ON"
            ,"FLUXO_DISC_NODETYPE"    ,"GAUSS-LOBATTO"
            ,"FLUXO_PARABOLIC"        ,"ON"
            ,"FLUXO_PARABOLIC_LIFTING","br2"
            ,"FLUXO_TESTCASE"         ,"mhd_equilibrium"
           ])
   
   if(not dbg ) : stat = test_fluxo(buildopts=options, case=caseID, project=pname, ntail = args.ntail ,\
                          run_test=TEST , mpi_procs = args.procs , keepdir=args.keepdir, err=builderr )
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
caseID=caseID+1
if(cases[0]==0 or (caseID in cases)) :
   pname="build_mhd_type1_GL_br1_TCmhd_noGLM_anisodiff"
   print "caseID: %d name: %s" % (caseID,pname)

   options=[]; options.extend(globopts) ; options.extend(baseopts)
   options.extend([
             "CMAKE_BUILD_TYPE"       ,"Debug"
            ,"FLUXO_DISCTYPE"         ,"1"
            ,"FLUXO_EQN_GLM"          ,"OFF"
            ,"FLUXO_EQN_ANISO_HEAT"   ,"ON"
            ,"FLUXO_DISC_NODETYPE"    ,"GAUSS-LOBATTO"
            ,"FLUXO_PARABOLIC"        ,"ON"
            ,"FLUXO_PARABOLIC_LIFTING","br1"
            ,"FLUXO_TESTCASE"         ,"mhd_equilibrium"
           ])
   
   if(not dbg ) : stat = test_fluxo(buildopts=options, case=caseID, project=pname, ntail = args.ntail ,\
                          run_test=TEST , mpi_procs = args.procs , keepdir=args.keepdir, err=builderr )
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
caseID=caseID+1
if(cases[0]==0 or (caseID in cases)) :
   pname="build_mhd_type1_nopara"
   print "caseID: %d name: %s" % (caseID,pname)

   options=[]; options.extend(globopts) ; options.extend(baseopts)
   options.extend([
             "CMAKE_BUILD_TYPE"       ,"Debug"
            ,"FLUXO_DISCTYPE"         ,"1"
            ,"FLUXO_EQN_GLM"          ,"ON"
            ,"FLUXO_DISC_NODETYPE"    ,"GAUSS-LOBATTO"
            ,"FLUXO_PARABOLIC"        ,"OFF"
            ,"FLUXO_TESTCASE"         ,"default"
           ])
   
   if(not dbg ) : stat = test_fluxo(buildopts=options, case=caseID, project=pname, ntail = args.ntail ,\
                          run_test=TEST , mpi_procs = args.procs , keepdir=args.keepdir, err=builderr )
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
caseID=caseID+1
if(cases[0]==0 or (caseID in cases)) :
   pname="build_mhd_type2_nopara_noGLM_cart"
   print "caseID: %d name: %s" % (caseID,pname)

   options=[]; options.extend(globopts) ; options.extend(baseopts)
   options.extend([
             "CMAKE_BUILD_TYPE"       ,"Debug"
            ,"FLUXO_DISCTYPE"         ,"2"
            ,"FLUXO_EQN_GLM"          ,"OFF"
            ,"FLUXO_DISC_CARTESIANFLUX","ON"           #=> runtest off
            ,"FLUXO_DISC_NODETYPE"    ,"GAUSS-LOBATTO"
            ,"FLUXO_PARABOLIC"        ,"OFF"
            ,"FLUXO_TESTCASE"         ,"default"
           ])
   
   if(not dbg ) : stat = test_fluxo(buildopts=options, case=caseID, project=pname, ntail = args.ntail ,\
                          run_test=None , mpi_procs = args.procs , keepdir=args.keepdir, err=builderr )
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~


#============================================================================
#FINAL ERROR HANDLING:
#============================================================================
if(len(builderr) > 0 ) :
   print "!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!"
   print "!!!!!!!    WARNING, following builds failed:     !!!!!!!!"
   print " "
   for line in builderr :
      print "--> "+line
   print " "
   print "... see log.caseID_project files"
   print "   and dirx_caseID_project folders."
   print " rerun build script with -case option to run again spefic caseIDs."
   print "!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!"
else :
   print "/////////////////////////////////////////////////////////"
   print " "
   print " ==> ALL BUILDS WERE SUCCESSFULL!"
   print " "
   print "/////////////////////////////////////////////////////////"

