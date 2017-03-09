#!/usr/bin/python
# -*- coding: utf-8 -*-

import os
import sys
import math
import shutil
import tempfile
import time
import argparse


########################################################################################################

def buildfluxo( buildopts=None , case=0, project="test" , ntail = 0 , mpi_procs = 1, keepdir=0, 
                err= None ):
   if len(buildopts) < 2 :
      print "error, nobuild options given"
      err.extend(["caseID=%6d ,project= %s" % (case,project)])
      return False
   # build directory
   builddir=("dirx_%d_%s" % (case,project))
   log_path=("../log.%d_%s" % (case,project))

   os.system('rm -rf '+builddir) #destroy directory from previous run
   os.system('mkdir '+builddir )
   os.chdir(builddir)

   allopts=" "
   for i in range((len(buildopts)/2)) :
      allopts=allopts+" -D"+buildopts[2*i]+"="+buildopts[2*i+1]+" " 

   print "===> OPTIONS:"
   print allopts

   logf = open(log_path, 'w')
   logf.write("===> OPTIONS: \n %s \n" % (allopts) )
 
   cmdconfig = "cmake ../../. "+allopts
   print "  "
   print "===> configure..."

   os.system(cmdconfig+" 2>stderr 1>stdout")
   stdout=open("stdout",'r').readlines()
   stderr=open("stderr",'r').readlines()

   cmake_err= False
   for line in stderr :
      if "ERROR" in line :
         cmake_err= True
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
   if(cmake_err) :
      print "===============================================  "
      print " !!!! ERROR IN CMAKE, NOT FINISHED CORRECTLY!!!!! "
      print "===============================================  "
      print "  "
      err.extend(["caseID=%6d ,project= %s" % (case,project)])
      return (not cmake_err)
   #MAKE
   print "  "
   print "===> make..."
   cmdmake = ("make -j %d VERBOSE=1" % (mpi_procs) )

   os.system(cmdmake+" 2>stderr 1>stdout")
   stdout=open("stdout",'r').readlines()
   stderr=open("stderr",'r').readlines()
   
   logf = open(log_path, 'a') #append
   build_success= False

   logf.write('MAKE: STDOUT FILE:\n')
   for line in stdout :
      logf.write(line)
      if "SUCCESS: FLUXO BUILD COMPLETE!" in line :
         build_success= True

   logf.write('MAKE: STDERR FILE:\n')
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

   os.chdir('../')
   print "===============================================  "
   if (build_success) :
      if(keepdir ==0) :
         os.system('rm -rf '+builddir) #destroy directory
      print " Build finished sucessfully."
   else :
      print " !!!! PROBLEM WITH BUILD, NOT FINISHED CORRECTLY!!!!! "
      err.extend(["caseID=%6d ,project= %s" % (case,project)])
   print "===============================================  "
   print "  "
   return build_success
  
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
parser.add_argument('-ntail', type=int, default=5,     help='    number of last line output of cmake/make screenlog (DEFAULT=5)')
parser.add_argument('-withmpi', type=int, default=1,   help="1 : DEFAULT ,compile with mpi\n"
                                                            "0 : compile without mpi  ")
parser.add_argument('-buildhdf5', type=int, default=1, help="1 : DEFAULT, build hdf5 locally,\n"
                                                            "0 : use external hdf5 (modules)")
parser.add_argument('-hostname', type=str, default="", help="    cmake hostname, only needed if compiling on a cluster" )
parser.add_argument('-keepdir', type=int, default=1,   help="1 : DEFAULT, keep all build directories,\n" 
                                                            "0 : delete sucessfull build directories")
parser.add_argument('-case', type=str, default='0',    help="0 : DEFAULT, run all cases,\n" 
                                                            "101,102-104 : list of specific cases to run (without spaces!) ")

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
#######################################################################
#first group, 100 < caseID <200
caseID=100 
baseopts=[
           "FLUXO_EQNSYSNAME"       ,"linearscalaradvection"
          ,"FLUXO_TESTCASE"         ,"default"
         ]
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
   
   if(not dbg ) : stat = buildfluxo(buildopts=options, case=caseID, project=pname, ntail = args.ntail ,\
                                    mpi_procs = args.procs , keepdir=args.keepdir, err=builderr )
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
   
   if(not dbg ) : stat = buildfluxo(buildopts=options, case=caseID, project=pname, ntail = args.ntail ,\
                                    mpi_procs = args.procs , keepdir=args.keepdir, err=builderr )
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
   
   if(not dbg ) : stat = buildfluxo(buildopts=options, case=caseID, project=pname, ntail = args.ntail ,\
                                    mpi_procs = args.procs , keepdir=args.keepdir, err=builderr )
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
   
   if(not dbg ) : stat = buildfluxo(buildopts=options, case=caseID, project=pname, ntail = args.ntail ,\
                                    mpi_procs = args.procs , keepdir=args.keepdir, err=builderr )
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
   
   if(not dbg ) : stat = buildfluxo(buildopts=options, case=caseID, project=pname, ntail = args.ntail ,\
                                    mpi_procs = args.procs , keepdir=args.keepdir, err=builderr )
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
   
   if(not dbg ) : stat = buildfluxo(buildopts=options, case=caseID, project=pname, ntail = args.ntail ,\
                                    mpi_procs = args.procs , keepdir=args.keepdir, err=builderr )
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
   
   if(not dbg ) : stat = buildfluxo(buildopts=options, case=caseID, project=pname, ntail = args.ntail ,\
                                    mpi_procs = args.procs , keepdir=args.keepdir, err=builderr )
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
caseID=caseID+1
if(cases[0]==0 or (caseID in cases)) :
   pname="build_linadv_type2_nopara_cart"
   print "caseID: %d name: %s" % (caseID,pname)

   options=[]; options.extend(globopts) ; options.extend(baseopts)
   options.extend([
             "CMAKE_BUILD_TYPE"       ,"Debug"
            ,"FLUXO_DISCTYPE"         ,"2"
            ,"FLUXO_DISC_CARTESIANFLUX","ON"
            ,"FLUXO_DISC_NODETYPE"    ,"GAUSS-LOBATTO"
            ,"FLUXO_PARABOLIC"        ,"OFF"
           ])
   
   if(not dbg ) : stat = buildfluxo(buildopts=options, case=caseID, project=pname, ntail = args.ntail ,\
                                    mpi_procs = args.procs , keepdir=args.keepdir, err=builderr )

######################################################################

#FINAL ERROR HANDLING:
if(len(builderr) > 0 ) :
  print "!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!"
  print "!!!!!!!    WARNING, following builds failed:     !!!!!!!!"
  print " "
  for line in builderr :
     print "--> "+line
  print " "
  print "... see log.caseID_project files"
  print "   and dirx_caseID_project folders."
  print "!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!"
else :
  print "/////////////////////////////////////////////////////////"
  print " "
  print " ==> ALL BUILDS WERE SUCCESSFULL!"
  print " "
  print "/////////////////////////////////////////////////////////"

