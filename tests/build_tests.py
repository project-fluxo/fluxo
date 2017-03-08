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

def buildfluxo( buildopts=None , project="test" , ntail = 0 , mpi_procs = 1, keepdir=0):
   # build directory
   builddir='dir_'+project
   os.system('rm -rf '+builddir) #destroy directory from previous run
   os.system('mkdir '+builddir )
   os.chdir(builddir)

   allopts=" "
   for i in range((len(buildopts)/2)) :
      allopts=allopts+" -D"+buildopts[2*i]+"="+buildopts[2*i+1]+" " 

   print "===> OPTIONS:"
   print allopts

   log_path="../log."+project
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
      logf.write(line)
      if "ERROR" in line :
         cmake_err= True
         sys.stdout.write("%s" % line)

   for line in stdout :
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
      return (not cmake_err)
   #MAKE
   print "  "
   print "===> make..."
   cmdmake = ("make -j %d VERBOSE=1" % (mpi_procs) )
   os.system(cmdmake+" 2>stderr 1>stdout")

   stdout=open("stdout",'r').readlines()
   stderr=open("stderr",'r').readlines()
   
   logf = open(log_path, 'a')
   build_success= False
   for line in stdout :
      logf.write(line)
      if "SUCCESS: FLUXO BUILD COMPLETE!" in line :
         build_success= True
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
   print "===============================================  "
   print "  "
   return build_success
  

########################################################################################################

# MAIN PROGRAM

########################################################################################################
parser = argparse.ArgumentParser(description='Tool to build/compile fluxo in different configurations',\
                                 formatter_class=argparse.RawTextHelpFormatter)
parser.add_argument('-p','--procs', type=int, default=1, help='number of processors used for the make')
parser.add_argument('-ntail', type=int, default=5, help='number of last line output of cmake/make screenlog')
parser.add_argument('-withmpi', type=int, default=1, help="(1) compile with mpi (default)\n"
                                                          "(0) compile without mpi  ")
parser.add_argument('-buildhdf5', type=int, default=1, help="(1) build hdf5 locally (default),\n"
                                                            "(0) use external hdf5 (modules)")
parser.add_argument('-keepdir', type=int, default=1, help="(1) keep all build directories (default),\n" 
                                                          "(0) delete sucessfull build directories")

args = parser.parse_args()


if(args.withmpi == 0) :
  MPIOPT="OFF"
else :
  MPIOPT="ON"

if(args.buildhdf5 == 0) :
  HDF5OPT="OFF"
else :
  HDF5OPT="ON"


builderr= "_"

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

options=[ "CMAKE_BUILD_TYPE"       ,"Release"
         ,"FLUXO_BUILD_MPI"        ,MPIOPT
         ,"FLUXO_BUILD_HDF5"       ,HDF5OPT
         ,"FLUXO_EQNSYSNAME"       ,"linearscalaradvection"
         ,"FLUXO_DISCTYPE"         ,"1"
         ,"FLUXO_DISC_NODETYPE"    ,"GAUSS-LOBATTO"
         ,"FLUXO_PARABOLIC"        ,"ON"
         ,"FLUXO_PARABOLIC_LIFTING","br1"
         ,"FLUXO_TESTCASE"         ,"default" 
        ]

pname="build_test_1"
stat = buildfluxo(buildopts=options, project=pname,\
                              ntail = args.ntail ,\
                              mpi_procs = args.procs , keepdir=args.keepdir )
if(not stat) :
  builderr= builderr+" "+pname
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
 

options=[ "CMAKE_BUILD_TYPE"       ,"Debug"
         ,"FLUXO_BUILD_MPI"        ,MPIOPT
         ,"FLUXO_BUILD_HDF5"       ,HDF5OPT
         ,"FLUXO_EQNSYSNAME"       ,"linearscalaradvection"
         ,"FLUXO_DISCTYPE"         ,"1"
         ,"FLUXO_DISC_NODETYPE"    ,"GAUSS-LOBATTO"
         ,"FLUXO_PARABOLIC"        ,"ON"
         ,"FLUXO_PARABOLIC_LIFTING","br1"
         ,"FLUXO_TESTCASE"         ,"default" 
        ]

pname="build_test_2"
stat = buildfluxo(buildopts=options, project=pname,\
                              ntail = args.ntail ,\
                              mpi_procs = args.procs , keepdir=args.keepdir )
if(not stat) :
  builderr= builderr+" "+pname
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
 

options=[ "CMAKE_BUILD_TYPE"       ,"Debug"
         ,"FLUXO_BUILD_MPI"        ,MPIOPT
         ,"FLUXO_BUILD_HDF5"       ,HDF5OPT
         ,"FLUXO_EQNSYSNAME"       ,"linearscalaradvection"
         ,"FLUXO_DISCTYPE"         ,"1"
         ,"FLUXO_DISC_NODETYPE"    ,"GAUSS"
         ,"FLUXO_PARABOLIC"        ,"ON"
         ,"FLUXO_PARABOLIC_LIFTING","br2"
         ,"FLUXO_TESTCASE"         ,"default" 
        ]

pname="build_test_3"
stat = buildfluxo(buildopts=options, project=pname,\
                              ntail = args.ntail ,\
                              mpi_procs = args.procs , keepdir=args.keepdir )
if(not stat) :
  builderr= builderr+" "+pname
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
 

options=[ "CMAKE_BUILD_TYPE"       ,"Debug"
         ,"FLUXO_BUILD_MPI"        ,MPIOPT
         ,"FLUXO_BUILD_HDF5"       ,HDF5OPT
         ,"FLUXO_EQNSYSNAME"       ,"linearscalaradvection"
         ,"FLUXO_DISCTYPE"         ,"2"
         ,"FLUXO_DISC_NODETYPE"    ,"GAUSS-LOBATTO"
         ,"FLUXO_PARABOLIC"        ,"ON"
         ,"FLUXO_PARABOLIC_LIFTING","br1"
         ,"FLUXO_TESTCASE"         ,"default" 
        ]

pname="build_test_4"
stat = buildfluxo(buildopts=options, project=pname,\
                              ntail = args.ntail ,\
                              mpi_procs = args.procs , keepdir=args.keepdir )
if(not stat) :
  builderr= builderr+" "+pname
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~



if(len(builderr) > 1 ) :
  print "!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!"
  print "!!!!! WARNING, following builds failed !!!!!"
  print builderr
  print "see log.build_test files.."
  print "!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!"

