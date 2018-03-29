#!/usr/bin/python
# -*- coding: utf-8 -*-
#==================================================================================================================================
# Copyright (c) 2016 - 2017 Gregor Gassner
# Copyright (c) 2016 - 2017 Florian Hindenlang
# Copyright (c) 2016 - 2017 Andrew Winters
#
# This file is part of FLUXO (github.com/project-fluxo/fluxo). FLUXO is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by the Free Software Foundation, either version 3 
# of the License, or (at your option) any later version.
#
# FLUXO is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY; without even the implied warranty
# of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU General Public License v3.0 for more details.
#
# You should have received a copy of the GNU General Public License along with FLUXO. If not, see <http://www.gnu.org/licenses/>.
#==================================================================================================================================

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
                stage=0 , run_test=None, mpi_procs = 1 , err= None ):
   if ( len(buildopts) < 2 ) :
      print( "error, nobuild options given"  )
      err.extend(["caseID=%6d ,project= %s <== build opts." % (case,project)])
      return False
   # build directory
   builddir=("dirx_%d_%s" % (case,project))
   log_path=("../log_%d_%s.txt" % (case,project))

   cwd = os.getcwd()  #current working directory
   if (stage < 2) :
     os.system('rm -rf '+builddir) #destroy directory from previous run
     os.system('mkdir '+builddir )

     os.chdir(builddir)

     allopts=" "
     for i in range(int(int(len(buildopts))/2)) :
        allopts=allopts+" -D"+buildopts[2*i]+"="+buildopts[2*i+1]+" " 
     
     logf = open(log_path, 'w')
     print( "===> OPTIONS:" )
     print( allopts )
     
     logf.write("===> OPTIONS: \n %s \n" % (allopts) )
     
     cmdconfig = "cmake ../../."+allopts
     print( "  " )
     print( "===> configure..." )
     
     os.system(cmdconfig+" 2>std.err 1>std.out")
     stdout=open("std.out",'r').readlines()
     stderr=open("std.err",'r').readlines()
     
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
        for line in stdout[nlines-ntail:nlines] :
           sys.stdout.write("%s" % line)
     if(not success) :
        print( "================================================ " )
        print( " !!!! ERROR IN CMAKE, NOT FINISHED CORRECTLY!!!! " )
        print( "================================================ " )
        print( "  " )
        err.extend(["caseID=%6d ,project= %s <=Cmake" % (case,project)])
        os.chdir(cwd)
        return success #=False
     else : 
        print( "================================================ " )
        print( " Cmake finished successfully.                    " )                
        print( "================================================ " )
        print( "  " )
     #endif (not success)

     #MAKE
     print( "  " )
     print( "===> make..." )
     cmdmake = ("make -j %d VERBOSE=1" % (mpi_procs) )
     
     os.system(cmdmake+" 2>std.err 1>std.out")
     stdout=open("std.out",'r').readlines()
     stderr=open("std.err",'r').readlines()
     
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
        for line in stdout[nlines-ntail:nlines] :
           sys.stdout.write("%s" % line)
     
     if (not success) :
        print( "==================================================== " )
        print( "!!!! PROBLEM WITH BUILD, NOT FINISHED CORRECTLY!!!!! " )
        err.extend(["caseID=%6d ,project= %s <=Make" % (case,project)] )
        print( "==================================================== " )
        print( "  " )
        os.chdir(cwd)
        return success #=False
     else :
        print( "==================================================== " )
        print( " Build finished sucessfully.                         " )                  
        print( "==================================================== " )
        print( "  " )
     #endif (not success)
      
   #endif (stage <2)
   
   # RUN TEST: stage =1 and 2
   if (stage > 0 ) :
     os.chdir(cwd)
     
     doruntest=False
     if( run_test ) :
        if( len(run_test) == 3 ) :  # run tests
           doruntest = True
     if(doruntest) :
        if ( not os.path.exists(builddir)) :
          print( "error, build dir does not exist: %s " % builddir  )
          err.extend(["caseID=%6d ,project= %s <== no builddir" % (case,project)])
          return False
        #RUN TEST
        print( "  " )
        print( "===> run test %s/%s ..." % (run_test[0],run_test[1]) )
        #change to test directory
        os.chdir(run_test[0])
        if (not os.path.isfile('./'+run_test[1])) :
           print( "!!!!!!!!  PARAMTERFILE %s/%s does not exist!" %(run_test[0],run_test[1]) )
           os.chdir(cwd)
           return False
     
        projectname = read_prm(run_test[1],'ProjectName')
        projectnamex=("%d_%s_%s" % (case,project,projectname))
        success = False
        try :
           [L2,Linf,PID] = execute(cwd+"/"+builddir+"/bin/fluxo", run_test[1], projectnamex,\
                                 [get_last_L2_error, get_last_Linf_error, get_cpu_per_dof],\
                                 log = True, ntail = ntail ,\
                                       mpi_procs = mpi_procs )
           if(Linf) :
              print( "   ... check Linf %s < %s ?" % (Linf[0],run_test[2]) )
              if(float(Linf[0]) < float(run_test[2])) :
                 success = True
        except :
           success = False
        if (not success) :
           print( "================================================== " )
           print( "!!!! PROBLEM WITH RUN, NOT FINISHED CORRECTLY!!!!! " )
           err.extend(["caseID=%6d ,project= %s <=Run" % (case,project)])
           print( "================================================== " )
           print( "  " )
           os.chdir(cwd)
           return success #=False
        else :
           print( "==================================================  " )
           print( " Run finished sucessfully.                          " )
           print( "==================================================  " )
           print( "  " )
     else :
        print( " (no test specified...) " )
        success = True #ignore this test
     #endif (not success)
   #endif (stage >1)

   # success still True
   os.chdir(cwd)
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
parser.add_argument('-p','--procs',type=int,default=1, help='    number of processors used for the make / execute (DEFAULT=1)')

parser.add_argument('-withmpi', type=int, default=1,   help="1 : DEFAULT ,compile with mpi\n"
                                                            "0 : compile without mpi  ")

parser.add_argument('-buildhdf5', type=int, default=0, help="0 : DEFAULT, use external hdf5 (modules),\n"
                                                            "1 : build hdf5 locally")

parser.add_argument('-hostname', type=str, default="", help="    cmake hostname, only needed if compiling on a cluster" )

parser.add_argument('-case', type=str, default='0',    help="0 : DEFAULT, run all cases,\n" 
                                                            "101,102-104 : list of specific cases to run (without spaces!) ")

parser.add_argument('-stage',type=int, default='0',    help="0 : DEFAULT, only build code\n"
                                                            "1 : build code and run with executable and parameterfile\n" 
                                                            "2 : only run examples (checks if code builds exist) ")

parser.add_argument('-ntail', type=int, default=5,     help='    number of last line output of cmake/make screenlog (DEFAULT=5)')


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
          ,"_BUILD_FLUXO_POST"      ,"ON"
         ]

# relative path from tests folder, parameterfile,Linf[0]<crit for success
TEST=[]
TEST.extend(["freestream","parameter_freestream_linadv.ini", "1.0e-10" ])
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
caseID=caseID+1
if(cases[0] ==0 or (caseID in cases)) :
   pname="build_linadv_release_type1_GL"
   print( "caseID: %d name: %s" % (caseID,pname) )

   options=[]; options.extend(globopts) ; options.extend(baseopts)
   options.extend([
             "CMAKE_BUILD_TYPE"       ,"Release"
            ,"FLUXO_DISCTYPE"         ,"1"
            ,"FLUXO_DISC_NODETYPE"    ,"GAUSS-LOBATTO"
            ,"FLUXO_PARABOLIC"        ,"ON"
            ,"FLUXO_PARABOLIC_LIFTING","br1"
           ])
   
   if(not dbg ) : stat = test_fluxo(buildopts=options, case=caseID, project=pname,ntail = args.ntail,\
                          stage=args.stage , run_test=TEST , mpi_procs = args.procs , err=builderr )
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
caseID=caseID+1
if(cases[0] ==0 or (caseID in cases)) :
   pname="build_linadv_release_type2"
   print( "caseID: %d name: %s" % (caseID,pname) )

   options=[]; options.extend(globopts) ; options.extend(baseopts)
   options.extend([
             "CMAKE_BUILD_TYPE"       ,"Release"
            ,"FLUXO_DISCTYPE"         ,"2"
            ,"FLUXO_DISC_NODETYPE"    ,"GAUSS-LOBATTO"
            ,"FLUXO_PARABOLIC"        ,"ON"
            ,"FLUXO_PARABOLIC_LIFTING","br1"
           ])
   
   if(not dbg ) : stat = test_fluxo(buildopts=options, case=caseID, project=pname, ntail = args.ntail ,\
                          stage=args.stage , run_test=TEST , mpi_procs = args.procs , err=builderr )
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
caseID=caseID+1
if(cases[0]==0 or (caseID in cases)) :
   pname="build_linadv_type2_br1"
   print( "caseID: %d name: %s" % (caseID,pname) )

   options=[]; options.extend(globopts) ; options.extend(baseopts)
   options.extend([
             "CMAKE_BUILD_TYPE"       ,"Debug"
            ,"FLUXO_DISCTYPE"         ,"2"
            ,"FLUXO_DISC_NODETYPE"    ,"GAUSS-LOBATTO"
            ,"FLUXO_PARABOLIC"        ,"ON"
            ,"FLUXO_PARABOLIC_LIFTING","br1"
           ])
   
   if(not dbg ) : stat = test_fluxo(buildopts=options, case=caseID, project=pname, ntail = args.ntail ,\
                          stage=args.stage , run_test=TEST , mpi_procs = args.procs , err=builderr )
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
caseID=caseID+1
if(cases[0]==0 or (caseID in cases)) :
   pname="build_linadv_type1_GL_br1"
   print( "caseID: %d name: %s" % (caseID,pname) )

   options=[]; options.extend(globopts) ; options.extend(baseopts)
   options.extend([
             "CMAKE_BUILD_TYPE"       ,"Debug"
            ,"FLUXO_DISCTYPE"         ,"1"
            ,"FLUXO_DISC_NODETYPE"    ,"GAUSS-LOBATTO"
            ,"FLUXO_PARABOLIC"        ,"ON"
            ,"FLUXO_PARABOLIC_LIFTING","br1"
           ])
   
   if(not dbg ) : stat = test_fluxo(buildopts=options, case=caseID, project=pname, ntail = args.ntail ,\
                          stage=args.stage , run_test=TEST , mpi_procs = args.procs , err=builderr )
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
caseID=caseID+1
if(cases[0]==0 or (caseID in cases)) :
   pname="build_linadv_type1_GL_br2"
   print( "caseID: %d name: %s" % (caseID,pname) )

   options=[]; options.extend(globopts) ; options.extend(baseopts)
   options.extend([
             "CMAKE_BUILD_TYPE"       ,"Debug"
            ,"FLUXO_DISCTYPE"         ,"1"
            ,"FLUXO_DISC_NODETYPE"    ,"GAUSS-LOBATTO"
            ,"FLUXO_PARABOLIC"        ,"ON"
            ,"FLUXO_PARABOLIC_LIFTING","br2"
           ])
   
   if(not dbg ) : stat = test_fluxo(buildopts=options, case=caseID, project=pname, ntail = args.ntail ,\
                           stage=args.stage , run_test=TEST , mpi_procs = args.procs , err=builderr )
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
caseID=caseID+1
if(cases[0]==0 or (caseID in cases)) :
   pname="build_linadv_type1_Gauss_br1"
   print( "caseID: %d name: %s" % (caseID,pname) )

   options=[]; options.extend(globopts) ; options.extend(baseopts)
   options.extend([
             "CMAKE_BUILD_TYPE"       ,"Debug"
            ,"FLUXO_DISCTYPE"         ,"1"
            ,"FLUXO_DISC_NODETYPE"    ,"GAUSS"
            ,"FLUXO_PARABOLIC"        ,"ON"
            ,"FLUXO_PARABOLIC_LIFTING","br1"
           ])
   
   if(not dbg ) : stat = test_fluxo(buildopts=options, case=caseID, project=pname, ntail = args.ntail ,\
                           stage=args.stage , run_test=TEST , mpi_procs = args.procs , err=builderr )
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
caseID=caseID+1
if(cases[0]==0 or (caseID in cases)) :
   pname="build_linadv_type1_Gauss_br2"
   print( "caseID: %d name: %s" % (caseID,pname) )

   options=[]; options.extend(globopts) ; options.extend(baseopts)
   options.extend([
             "CMAKE_BUILD_TYPE"       ,"Debug"
            ,"FLUXO_DISCTYPE"         ,"1"
            ,"FLUXO_DISC_NODETYPE"    ,"GAUSS"
            ,"FLUXO_PARABOLIC"        ,"ON"
            ,"FLUXO_PARABOLIC_LIFTING","br2"
           ])
   
   if(not dbg ) : stat = test_fluxo(buildopts=options, case=caseID, project=pname, ntail = args.ntail ,\
                          stage=args.stage , run_test=TEST , mpi_procs = args.procs , err=builderr )
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
caseID=caseID+1
if(cases[0]==0 or (caseID in cases)) :
   pname="build_linadv_type2_nopara"
   print( "caseID: %d name: %s" % (caseID,pname) )

   options=[]; options.extend(globopts) ; options.extend(baseopts)
   options.extend([
             "CMAKE_BUILD_TYPE"       ,"Debug"
            ,"FLUXO_DISCTYPE"         ,"2"
            ,"FLUXO_DISC_NODETYPE"    ,"GAUSS-LOBATTO"
            ,"FLUXO_PARABOLIC"        ,"OFF"
           ])
   
   if(not dbg ) : stat = test_fluxo(buildopts=options, case=caseID, project=pname, ntail = args.ntail ,\
                          stage=args.stage , run_test=TEST , mpi_procs = args.procs , err=builderr )
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
# relative path from tests folder, parameterfile,Linf[0]<crit for success
TEST=[]
# TEST.extend(["freestream","parameter_freestream_maxwell.ini", "1.0e-10" ])
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
caseID=caseID+1
if(cases[0] ==0 or (caseID in cases)) :
   pname="build_maxwell_type1_GL"
   print( "caseID: %d name: %s" % (caseID,pname) )

   options=[]; options.extend(globopts) ; options.extend(baseopts)
   options.extend([
             "CMAKE_BUILD_TYPE"       ,"Debug"
            ,"FLUXO_DISCTYPE"         ,"1"
            ,"FLUXO_DISC_NODETYPE"    ,"GAUSS-LOBATTO"
            ,"FLUXO_PARABOLIC"        ,"OFF"
           ])
   
   if(not dbg ) : stat = test_fluxo(buildopts=options, case=caseID, project=pname, ntail = args.ntail ,\
                          stage=args.stage , run_test=TEST , mpi_procs = args.procs , err=builderr )
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
caseID=caseID+1
if(cases[0] ==0 or (caseID in cases)) :
   pname="build_maxwell_type2"
   print( "caseID: %d name: %s" % (caseID,pname) )

   options=[]; options.extend(globopts) ; options.extend(baseopts)
   options.extend([
             "CMAKE_BUILD_TYPE"       ,"Debug"
            ,"FLUXO_DISCTYPE"         ,"2"
            ,"FLUXO_DISC_NODETYPE"    ,"GAUSS-LOBATTO"
            ,"FLUXO_PARABOLIC"        ,"OFF"
           ])
   
   if(not dbg ) : stat = test_fluxo(buildopts=options, case=caseID, project=pname, ntail = args.ntail ,\
                          stage=args.stage , run_test=TEST , mpi_procs = args.procs , err=builderr )
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
# relative path from tests folder, parameterfile,Linf[0]<crit for success
TEST=[]
TEST.extend(["freestream","parameter_freestream_mhd.ini", "1.0e-10" ])
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
caseID=caseID+1
if(cases[0] ==0 or (caseID in cases)) :
   pname="build_mhd_release_type1_GL_br1"
   print( "caseID: %d name: %s" % (caseID,pname) )

   options=[]; options.extend(globopts) ; options.extend(baseopts)
   options.extend([
             "CMAKE_BUILD_TYPE"       ,"Release"
            ,"FLUXO_DISCTYPE"         ,"1"
            ,"FLUXO_EQN_GLM"          ,"ON"
            ,"FLUXO_EQN_NONCONS"      ,"OFF"
            ,"FLUXO_DISC_NODETYPE"    ,"GAUSS-LOBATTO"
            ,"FLUXO_PARABOLIC"        ,"ON"
            ,"FLUXO_PARABOLIC_LIFTING","br1"
            ,"FLUXO_TESTCASE"         ,"default"
           ])
   
   if(not dbg ) : stat = test_fluxo(buildopts=options, case=caseID, project=pname, ntail = args.ntail ,\
                          stage=args.stage , run_test=TEST , mpi_procs = args.procs , err=builderr )
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
caseID=caseID+1
if(cases[0] ==0 or (caseID in cases)) :
   pname="build_mhd_release_type1_Gauss_br2"
   print( "caseID: %d name: %s" % (caseID,pname) )

   options=[]; options.extend(globopts) ; options.extend(baseopts)
   options.extend([
             "CMAKE_BUILD_TYPE"       ,"Release"
            ,"FLUXO_DISCTYPE"         ,"1"
            ,"FLUXO_EQN_GLM"          ,"ON"
            ,"FLUXO_EQN_NONCONS"      ,"OFF"
            ,"FLUXO_DISC_NODETYPE"    ,"GAUSS"
            ,"FLUXO_PARABOLIC"        ,"ON"
            ,"FLUXO_PARABOLIC_LIFTING","br2"
            ,"FLUXO_TESTCASE"         ,"default"
           ])
   
   if(not dbg ) : stat = test_fluxo(buildopts=options, case=caseID, project=pname, ntail = args.ntail ,\
                          stage=args.stage , run_test=TEST , mpi_procs = args.procs , err=builderr )
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
caseID=caseID+1
if(cases[0] ==0 or (caseID in cases)) :
   pname="build_mhd_release_type2_br2"
   print( "caseID: %d name: %s" % (caseID,pname) )

   options=[]; options.extend(globopts) ; options.extend(baseopts)
   options.extend([
             "CMAKE_BUILD_TYPE"       ,"Release"
            ,"FLUXO_DISCTYPE"         ,"2"
            ,"FLUXO_EQN_GLM"          ,"ON"
            ,"FLUXO_EQN_NONCONS"      ,"OFF"
            ,"FLUXO_DISC_NODETYPE"    ,"GAUSS-LOBATTO"
            ,"FLUXO_PARABOLIC"        ,"ON"
            ,"FLUXO_PARABOLIC_LIFTING","br2"
            ,"FLUXO_TESTCASE"         ,"default"
           ])
   
   if(not dbg ) : stat = test_fluxo(buildopts=options, case=caseID, project=pname, ntail = args.ntail ,\
                          stage=args.stage , run_test=TEST , mpi_procs = args.procs , err=builderr )
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
caseID=caseID+1
if(cases[0]==0 or (caseID in cases)) :
   pname="build_mhd_type2_br1"
   print( "caseID: %d name: %s" % (caseID,pname) )

   options=[]; options.extend(globopts) ; options.extend(baseopts)
   options.extend([
             "CMAKE_BUILD_TYPE"       ,"Debug"
            ,"FLUXO_DISCTYPE"         ,"2"
            ,"FLUXO_EQN_GLM"          ,"ON"
            ,"FLUXO_EQN_NONCONS"      ,"ON"
            ,"FLUXO_DISC_NODETYPE"    ,"GAUSS-LOBATTO"
            ,"FLUXO_PARABOLIC"        ,"ON"
            ,"FLUXO_PARABOLIC_LIFTING","br1"
            ,"FLUXO_TESTCASE"         ,"default"
           ])
   
   if(not dbg ) : stat = test_fluxo(buildopts=options, case=caseID, project=pname, ntail = args.ntail ,\
                          stage=args.stage , run_test=TEST , mpi_procs = args.procs , err=builderr )
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
caseID=caseID+1
if(cases[0]==0 or (caseID in cases)) :
   pname="build_mhd_type1_GL_br2"
   print( "caseID: %d name: %s" % (caseID,pname) )

   options=[]; options.extend(globopts) ; options.extend(baseopts)
   options.extend([
             "CMAKE_BUILD_TYPE"       ,"Debug"
            ,"FLUXO_DISCTYPE"         ,"1"
            ,"FLUXO_EQN_GLM"          ,"ON"
            ,"FLUXO_EQN_NONCONS"      ,"OFF"
            ,"FLUXO_DISC_NODETYPE"    ,"GAUSS-LOBATTO"
            ,"FLUXO_PARABOLIC"        ,"ON"
            ,"FLUXO_PARABOLIC_LIFTING","br2"
            ,"FLUXO_TESTCASE"         ,"default"
           ])
   
   if(not dbg ) : stat = test_fluxo(buildopts=options, case=caseID, project=pname, ntail = args.ntail ,\
                          stage=args.stage , run_test=TEST , mpi_procs = args.procs , err=builderr )
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
caseID=caseID+1
if(cases[0]==0 or (caseID in cases)) :
   pname="build_mhd_type1_Gauss_br1"
   print( "caseID: %d name: %s" % (caseID,pname) )

   options=[]; options.extend(globopts) ; options.extend(baseopts)
   options.extend([
             "CMAKE_BUILD_TYPE"       ,"Debug"
            ,"FLUXO_DISCTYPE"         ,"1"
            ,"FLUXO_EQN_GLM"          ,"ON"
            ,"FLUXO_EQN_NONCONS"      ,"OFF"
            ,"FLUXO_DISC_NODETYPE"    ,"GAUSS"
            ,"FLUXO_PARABOLIC"        ,"ON"
            ,"FLUXO_PARABOLIC_LIFTING","br1"
            ,"FLUXO_PARABOLIC_LIFTING_VAR","cons_var"
            ,"FLUXO_TESTCASE"         ,"default"
           ])
   
   if(not dbg ) : stat = test_fluxo(buildopts=options, case=caseID, project=pname, ntail = args.ntail ,\
                          stage=args.stage , run_test=TEST , mpi_procs = args.procs , err=builderr )
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
caseID=caseID+1
if(cases[0]==0 or (caseID in cases)) :
   pname="build_mhd_type1_GL_br2_TCmhd"
   print( "caseID: %d name: %s" % (caseID,pname) )

   options=[]; options.extend(globopts) ; options.extend(baseopts)
   options.extend([
             "CMAKE_BUILD_TYPE"       ,"Debug"
            ,"FLUXO_DISCTYPE"         ,"1"
            ,"FLUXO_EQN_GLM"          ,"ON"
            ,"FLUXO_EQN_NONCONS"      ,"OFF"
            ,"FLUXO_DISC_NODETYPE"    ,"GAUSS-LOBATTO"
            ,"FLUXO_PARABOLIC"        ,"ON"
            ,"FLUXO_PARABOLIC_LIFTING","br2"
            ,"FLUXO_PARABOLIC_LIFTING_VAR","prim_var"
            ,"FLUXO_TESTCASE"         ,"mhd_equilibrium"
           ])
   
   if(not dbg ) : stat = test_fluxo(buildopts=options, case=caseID, project=pname, ntail = args.ntail ,\
                          stage=args.stage , run_test=TEST , mpi_procs = args.procs , err=builderr )
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
caseID=caseID+1
if(cases[0]==0 or (caseID in cases)) :
   pname="build_mhd_type1_GL_br1_TCmhd_noGLM_anisodiff"
   print( "caseID: %d name: %s" % (caseID,pname) )

   options=[]; options.extend(globopts) ; options.extend(baseopts)
   options.extend([
             "CMAKE_BUILD_TYPE"       ,"Debug"
            ,"FLUXO_DISCTYPE"         ,"1"
            ,"FLUXO_EQN_GLM"          ,"OFF"
            ,"FLUXO_EQN_NONCONS"      ,"OFF"
            ,"FLUXO_EQN_ANISO_HEAT"   ,"ON"
            ,"FLUXO_DISC_NODETYPE"    ,"GAUSS-LOBATTO"
            ,"FLUXO_PARABOLIC"        ,"ON"
            ,"FLUXO_PARABOLIC_LIFTING","br1"
            ,"FLUXO_PARABOLIC_LIFTING_VAR","cons_var"
            ,"FLUXO_TESTCASE"         ,"mhd_equilibrium"
           ])
   
   if(not dbg ) : stat = test_fluxo(buildopts=options, case=caseID, project=pname, ntail = args.ntail ,\
                          stage=args.stage , run_test=TEST , mpi_procs = args.procs , err=builderr )
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
caseID=caseID+1
if(cases[0]==0 or (caseID in cases)) :
   pname="build_mhd_type2_GL_NONCONS_br1_entropyLiftVar"
   print( "caseID: %d name: %s" % (caseID,pname) )

   options=[]; options.extend(globopts) ; options.extend(baseopts)
   options.extend([
             "CMAKE_BUILD_TYPE"       ,"Debug"
            ,"FLUXO_DISCTYPE"         ,"2"
            ,"FLUXO_EQN_GLM"          ,"ON"
            ,"FLUXO_EQN_NONCONS"      ,"ON"
            ,"FLUXO_DISC_NODETYPE"    ,"GAUSS-LOBATTO"
            ,"FLUXO_PARABOLIC"        ,"ON"
            ,"FLUXO_PARABOLIC_LIFTING","br1"
            ,"FLUXO_PARABOLIC_LIFTING_VAR","entropy_var"
            ,"FLUXO_TESTCASE"         ,"default"
           ])
   
   if(not dbg ) : stat = test_fluxo(buildopts=options, case=caseID, project=pname, ntail = args.ntail ,\
                          stage=args.stage , run_test=TEST , mpi_procs = args.procs , err=builderr )
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
caseID=caseID+1
if(cases[0]==0 or (caseID in cases)) :
   pname="build_mhd_type2_NONCONS_noGLM_nopara"
   print( "caseID: %d name: %s" % (caseID,pname) )

   options=[]; options.extend(globopts) ; options.extend(baseopts)
   options.extend([
             "CMAKE_BUILD_TYPE"       ,"Debug"
            ,"FLUXO_DISCTYPE"         ,"2"
            ,"FLUXO_EQN_GLM"          ,"OFF"
            ,"FLUXO_EQN_NONCONS"      ,"ON"
            ,"FLUXO_DISC_NODETYPE"    ,"GAUSS-LOBATTO"
            ,"FLUXO_PARABOLIC"        ,"OFF"
            ,"FLUXO_TESTCASE"         ,"default"
           ])
   
   if(not dbg ) : stat = test_fluxo(buildopts=options, case=caseID, project=pname, ntail = args.ntail ,\
                          stage=args.stage , run_test=TEST , mpi_procs = args.procs , err=builderr )
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
caseID=caseID+1
if(cases[0]==0 or (caseID in cases)) :
   pname="build_mhd_type2_GL_GLM_nopara"
   print( "caseID: %d name: %s" % (caseID,pname) )

   options=[]; options.extend(globopts) ; options.extend(baseopts)
   options.extend([
             "CMAKE_BUILD_TYPE"       ,"Debug"
            ,"FLUXO_DISCTYPE"         ,"2"
            ,"FLUXO_EQN_GLM"          ,"ON"
            ,"FLUXO_EQN_NONCONS"      ,"OFF"
            ,"FLUXO_DISC_NODETYPE"    ,"GAUSS-LOBATTO"
            ,"FLUXO_PARABOLIC"        ,"OFF"
            ,"FLUXO_TESTCASE"         ,"default"
           ])
   
   if(not dbg ) : stat = test_fluxo(buildopts=options, case=caseID, project=pname, ntail = args.ntail ,\
                          stage=args.stage , run_test=TEST , mpi_procs = args.procs , err=builderr )
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
caseID=caseID+1
if(cases[0]==0 or (caseID in cases)) :
   pname="build_mhd_type2_br1_prim_var_noGLM"
   print( "caseID: %d name: %s" % (caseID,pname) )

   options=[]; options.extend(globopts) ; options.extend(baseopts)
   options.extend([
             "CMAKE_BUILD_TYPE"       ,"Debug"
            ,"FLUXO_DISCTYPE"         ,"2"
            ,"FLUXO_EQN_GLM"          ,"OFF"
            ,"FLUXO_EQN_NONCONS"      ,"ON"
            ,"FLUXO_DISC_NODETYPE"    ,"GAUSS-LOBATTO"
            ,"FLUXO_PARABOLIC"        ,"ON"
            ,"FLUXO_PARABOLIC_LIFTING","br1"
            ,"FLUXO_PARABOLIC_LIFTING_VAR","prim_var"
            ,"FLUXO_TESTCASE"         ,"default"
           ])
   
   if(not dbg ) : stat = test_fluxo(buildopts=options, case=caseID, project=pname, ntail = args.ntail ,\
                          stage=args.stage , run_test=TEST , mpi_procs = args.procs , err=builderr )
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~


#============================================================================
#============================================================================
#fourth group, Navierstokes, 400 < caseID <500
#============================================================================
#============================================================================
caseID=400 
baseopts=[
           "FLUXO_EQNSYSNAME"       ,"navierstokes"
         ]
# relative path from tests folder, parameterfile,Linf[0]<crit for success
TEST=[]
TEST.extend(["freestream","parameter_freestream_navierstokes.ini", "1.0e-10" ])
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
caseID=caseID+1
if(cases[0] ==0 or (caseID in cases)) :
   pname="build_navierstokes_release_type1_br1_GL"
   print( "caseID: %d name: %s" % (caseID,pname) )

   options=[]; options.extend(globopts) ; options.extend(baseopts)
   options.extend([
             "CMAKE_BUILD_TYPE"       ,"Release"
            ,"FLUXO_DISCTYPE"         ,"1"
            ,"FLUXO_DISC_NODETYPE"    ,"GAUSS-LOBATTO"
            ,"FLUXO_PARABOLIC"        ,"ON"
            ,"FLUXO_PARABOLIC_LIFTING","br1"
            ,"FLUXO_PARABOLIC_LIFTING_VAR","cons_var"
            ,"FLUXO_TESTCASE"         ,"default"
           ])
   
   if(not dbg ) : stat = test_fluxo(buildopts=options, case=caseID, project=pname, ntail = args.ntail ,\
                          stage=args.stage , run_test=TEST , mpi_procs = args.procs , err=builderr )
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
caseID=caseID+1
if(cases[0] ==0 or (caseID in cases)) :
   pname="build_navierstokes_release_type2_br2"
   print( "caseID: %d name: %s" % (caseID,pname) )

   options=[]; options.extend(globopts) ; options.extend(baseopts)
   options.extend([
             "CMAKE_BUILD_TYPE"       ,"Release"
            ,"FLUXO_DISCTYPE"         ,"2"
            ,"FLUXO_DISC_NODETYPE"    ,"GAUSS-LOBATTO"
            ,"FLUXO_PARABOLIC"        ,"ON"
            ,"FLUXO_PARABOLIC_LIFTING","br2"
            ,"FLUXO_PARABOLIC_LIFTING_VAR","entropy_var"
            ,"FLUXO_TESTCASE"         ,"default"
           ])
   
   if(not dbg ) : stat = test_fluxo(buildopts=options, case=caseID, project=pname, ntail = args.ntail ,\
                          stage=args.stage , run_test=TEST , mpi_procs = args.procs , err=builderr )
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
caseID=caseID+1
if(cases[0]==0 or (caseID in cases)) :
   pname="build_navierstokes_type1_br1"
   print( "caseID: %d name: %s" % (caseID,pname) )

   options=[]; options.extend(globopts) ; options.extend(baseopts)
   options.extend([
             "CMAKE_BUILD_TYPE"       ,"Debug"
            ,"FLUXO_DISCTYPE"         ,"1"
            ,"FLUXO_DISC_NODETYPE"    ,"GAUSS-LOBATTO"
            ,"FLUXO_PARABOLIC"        ,"ON"
            ,"FLUXO_PARABOLIC_LIFTING","br1"
            ,"FLUXO_PARABOLIC_LIFTING_VAR","entropy_var"
            ,"FLUXO_TESTCASE"         ,"default"
           ])
   
   if(not dbg ) : stat = test_fluxo(buildopts=options, case=caseID, project=pname, ntail = args.ntail ,\
                          stage=args.stage , run_test=TEST , mpi_procs = args.procs , err=builderr )
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
caseID=caseID+1
if(cases[0]==0 or (caseID in cases)) :
   pname="build_navierstokes_type2_br2"
   print( "caseID: %d name: %s" % (caseID,pname) )

   options=[]; options.extend(globopts) ; options.extend(baseopts)
   options.extend([
             "CMAKE_BUILD_TYPE"       ,"Debug"
            ,"FLUXO_DISCTYPE"         ,"2"
            ,"FLUXO_DISC_NODETYPE"    ,"GAUSS-LOBATTO"
            ,"FLUXO_PARABOLIC"        ,"ON"
            ,"FLUXO_PARABOLIC_LIFTING","br1"
            ,"FLUXO_PARABOLIC_LIFTING_VAR","prim_var"
           ])
   
   if(not dbg ) : stat = test_fluxo(buildopts=options, case=caseID, project=pname, ntail = args.ntail ,\
                          stage=args.stage , run_test=TEST , mpi_procs = args.procs , err=builderr )
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
caseID=caseID+1
if(cases[0]==0 or (caseID in cases)) :
   pname="build_navierstokes_type1_Gauss_br1_viscsuth"
   print( "caseID: %d name: %s" % (caseID,pname) )

   options=[]; options.extend(globopts) ; options.extend(baseopts)
   options.extend([
             "CMAKE_BUILD_TYPE"       ,"Debug"
            ,"FLUXO_DISCTYPE"         ,"1"
            ,"FLUXO_EQN_VISCOSITY"    ,"sutherland"
            ,"FLUXO_DISC_NODETYPE"    ,"GAUSS"
            ,"FLUXO_PARABOLIC"        ,"ON"
            ,"FLUXO_PARABOLIC_LIFTING","br1"
            ,"FLUXO_PARABOLIC_LIFTING_VAR","cons_var"
           ])
   
   if(not dbg ) : stat = test_fluxo(buildopts=options, case=caseID, project=pname, ntail = args.ntail ,\
                          stage=args.stage , run_test=TEST , mpi_procs = args.procs , err=builderr )
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
caseID=caseID+1
if(cases[0]==0 or (caseID in cases)) :
   pname="build_navierstokes_type1_GL_br2_viscpow"
   print( "caseID: %d name: %s" % (caseID,pname) )

   options=[]; options.extend(globopts) ; options.extend(baseopts)
   options.extend([
             "CMAKE_BUILD_TYPE"       ,"Debug"
            ,"FLUXO_DISCTYPE"         ,"1"
            ,"FLUXO_EQN_VISCOSITY"    ,"powerlaw"
            ,"FLUXO_DISC_NODETYPE"    ,"GAUSS-LOBATTO"
            ,"FLUXO_PARABOLIC"        ,"ON"
            ,"FLUXO_PARABOLIC_LIFTING","br2"
            ,"FLUXO_PARABOLIC_LIFTING_VAR","cons_var"
           ])
   
   if(not dbg ) : stat = test_fluxo(buildopts=options, case=caseID, project=pname, ntail = args.ntail ,\
                          stage=args.stage , run_test=TEST , mpi_procs = args.procs , err=builderr )
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
caseID=caseID+1
if(cases[0]==0 or (caseID in cases)) :
   pname="build_navierstokes_type1_GL_nopara"
   print( "caseID: %d name: %s" % (caseID,pname) )

   options=[]; options.extend(globopts) ; options.extend(baseopts)
   options.extend([
             "CMAKE_BUILD_TYPE"       ,"Debug"
            ,"FLUXO_DISCTYPE"         ,"1"
            ,"FLUXO_DISC_NODETYPE"    ,"GAUSS-LOBATTO"
            ,"FLUXO_PARABOLIC"        ,"OFF"
           ])
   
   if(not dbg ) : stat = test_fluxo(buildopts=options, case=caseID, project=pname, ntail = args.ntail ,\
                          stage=args.stage , run_test=TEST , mpi_procs = args.procs , err=builderr )
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
caseID=caseID+1
if(cases[0]==0 or (caseID in cases)) :
   pname="build_navierstokes_type2_nopara"
   print( "caseID: %d name: %s" % (caseID,pname) )

   options=[]; options.extend(globopts) ; options.extend(baseopts)
   options.extend([
             "CMAKE_BUILD_TYPE"       ,"Debug"
            ,"FLUXO_DISCTYPE"         ,"2"
            ,"FLUXO_DISC_NODETYPE"    ,"GAUSS-LOBATTO"
            ,"FLUXO_PARABOLIC"        ,"OFF"
            ,"FLUXO_TESTCASE"         ,"default"
           ])
   
   if(not dbg ) : stat = test_fluxo(buildopts=options, case=caseID, project=pname, ntail = args.ntail ,\
                          stage=args.stage, run_test=TEST , mpi_procs = args.procs , err=builderr )
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
caseID=caseID+1
if(cases[0]==0 or (caseID in cases)) :
   pname="build_navierstokes_type1_GL_nopara_TC_angmom"
   print( "caseID: %d name: %s" % (caseID,pname) )

   options=[]; options.extend(globopts) ; options.extend(baseopts)
   options.extend([
             "CMAKE_BUILD_TYPE"       ,"Debug"
            ,"FLUXO_DISCTYPE"         ,"1"
            ,"FLUXO_DISC_NODETYPE"    ,"GAUSS-LOBATTO"
            ,"FLUXO_PARABOLIC"        ,"OFF"
            ,"FLUXO_TESTCASE"         ,"ns_angularmomentum"
           ])
   
   if(not dbg ) : stat = test_fluxo(buildopts=options, case=caseID, project=pname, ntail = args.ntail ,\
                          stage=args.stage , run_test=TEST , mpi_procs = args.procs , err=builderr )

#============================================================================
#FINAL ERROR HANDLING:
#============================================================================
if(len(builderr) > 0 ) :
   print( "!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!" )
   print( "!!!!!!!    WARNING, following builds failed:     !!!!!!!!" )
   print( " " )
   for line in builderr :
      print( "--> "+line )
   print( " " )
   print( "... see log.caseID_project files" )
   print( "   and dirx_caseID_project folders." )
   print( " rerun build script with -case option to run again spefic caseIDs." )
   print( "!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!" )
   sys.exit(100)
else :
   print( "/////////////////////////////////////////////////////////" )
   print( " " )
   if ( args.stage == 0) :
     print( " ==> ALL BUILDS SUCCESSFULL!" )
   if ( args.stage == 1) :
     print( " ==> ALL BUILDS AND RUNS SUCCESSFULL!" )
   if ( args.stage == 2) :
     print( " ==> ALL RUNS SUCCESSFULL!" )
   print( " " )
   print( "/////////////////////////////////////////////////////////" )
   sys.exit(0)

