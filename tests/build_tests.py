#!/usr/bin/python
# -*- coding: utf-8 -*-

import os
import sys
import math
import shutil
import tempfile
import time
import argparse
import subprocess


########################################################################################################

def buildfluxo( buildopts=None , log = True , log_path="out" , ntail = 0 , mpi_procs = 1):
   allopts=" "
   for i in range((len(buildopts)/2)) :
      allopts=allopts+" -D"+buildopts[2*i]+"="+buildopts[2*i+1]+" " 

   print "===> OPTIONS:"
   print allopts
   if log :
      f = open(log_path, 'w')
      f.write("===> OPTIONS: \n %s \n" % (allopts) )

   cmdconfig = "cmake ../../. "+allopts
   print "  "
   print "===> configure..."
   pc = subprocess.Popen(cmdconfig, stdout = subprocess.PIPE, stderr = subprocess.PIPE, shell=True)
   lines = []
   while pc.poll() is None :
      line = pc.stdout.readlines()
      if line :
         lines.extend(line)
   if log :
      for line in lines :
         f.write(line)

   if (ntail > 0) :
      nlines=len(lines)
      i=0
      for line in lines :
         i=i+1
         if ( i > (nlines -ntail) ) :
            print line,
   #MAKE
   print "  "
   print "===> make..."
   cmdmake = ("make -j %d VERBOSE=1" % (mpi_procs) )
   pm = subprocess.Popen(cmdmake, stdout = subprocess.PIPE, stderr = subprocess.PIPE, shell=True)
   lines = []
   while pm.poll() is None :
      line = pm.stdout.readlines()
      if line :
         lines.extend(line)
   if log :
      for line in lines :
         f.write(line)
      f.close()
   if (ntail > 0) :
      nlines=len(lines)
      i=0
      for line in lines :
         i=i+1
         if ( i > (nlines -ntail) ) :
            print line,
   build_success= False
   for l in lines[-30:] :
      if "SUCCESS: FLUXO BUILD COMPLETE!" in l :
         build_success= True
      if(build_success) :
        break
   print "===============================================  "
   if (build_success) :
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
parser.add_argument('-eqnsys', type=int, default=0, help="equation system to be used: \n"
                                                         "(0) linearscalaradvection (default) ,\n"
                                                         "(1) maxwell\n" )

args = parser.parse_args()

os.system('mkdir -p dir_build_test')
os.chdir('dir_build_test')

if(args.withmpi == 0) :
  MPIOPT="OFF"
else :
  MPIOPT="ON"

if(args.buildhdf5 == 0) :
  HDF5OPT="OFF"
else :
  HDF5OPT="ON"

if(args.eqnsys == 0) :
  EQNSYS="linearscalaradvection"
elif(args.eqnsys == 1) :
  EQNSYS="maxwell"
else :
  EQNSYS="linearscalaradvection"

builderr= "_"

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
outfile="../log.build_test_1"

options=[ "CMAKE_BUILD_TYPE"       ,"Release"
         ,"FLUXO_BUILD_MPI"        ,MPIOPT
         ,"FLUXO_BUILD_HDF5"       ,HDF5OPT
         ,"FLUXO_EQNSYSNAME"       ,EQNSYS
         ,"FLUXO_DISCTYPE"         ,"1"
         ,"FLUXO_DISC_NODETYPE"    ,"GAUSS-LOBATTO"
         ,"FLUXO_PARABOLIC"        ,"ON"
         ,"FLUXO_PARABOLIC_LIFTING","br1"
         ,"FLUXO_TESTCASE"         ,"default" 
        ]

stat = buildfluxo(buildopts=options,log= True, log_path=outfile,\
                              ntail = args.ntail ,\
                              mpi_procs = args.procs )
if(not stat) :
  builderr= builderr+" 1"
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
 
outfile="../log.build_test_2"

options=[ "CMAKE_BUILD_TYPE"       ,"Debug"
         ,"FLUXO_BUILD_MPI"        ,MPIOPT
         ,"FLUXO_BUILD_HDF5"       ,HDF5OPT
         ,"FLUXO_EQNSYSNAME"       ,EQNSYS
         ,"FLUXO_DISCTYPE"         ,"1"
         ,"FLUXO_DISC_NODETYPE"    ,"GAUSS-LOBATTO"
         ,"FLUXO_PARABOLIC"        ,"ON"
         ,"FLUXO_PARABOLIC_LIFTING","br1"
         ,"FLUXO_TESTCASE"         ,"default" 
        ]

stat = buildfluxo(buildopts=options,log= True, log_path=outfile,\
                              ntail = args.ntail ,\
                              mpi_procs = args.procs )
if(not stat) :
  builderr= builderr+" 2"
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
 
outfile="../log.build_test_3"

options=[ "CMAKE_BUILD_TYPE"       ,"Debug"
         ,"FLUXO_BUILD_MPI"        ,MPIOPT
         ,"FLUXO_BUILD_HDF5"       ,HDF5OPT
         ,"FLUXO_EQNSYSNAME"       ,EQNSYS
         ,"FLUXO_DISCTYPE"         ,"1"
         ,"FLUXO_DISC_NODETYPE"    ,"GAUSS"
         ,"FLUXO_PARABOLIC"        ,"ON"
         ,"FLUXO_PARABOLIC_LIFTING","br2"
         ,"FLUXO_TESTCASE"         ,"default" 
        ]

stat = buildfluxo(buildopts=options,log= True, log_path=outfile,\
                              ntail = args.ntail ,\
                              mpi_procs = args.procs )
if(not stat) :
  builderr= builderr+" 3"
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
 
outfile="../log.build_test_4"

options=[ "CMAKE_BUILD_TYPE"       ,"Debug"
         ,"FLUXO_BUILD_MPI"        ,MPIOPT
         ,"FLUXO_BUILD_HDF5"       ,HDF5OPT
         ,"FLUXO_EQNSYSNAME"       ,EQNSYS
         ,"FLUXO_DISCTYPE"         ,"2"
         ,"FLUXO_DISC_NODETYPE"    ,"GAUSS-LOBATTO"
         ,"FLUXO_PARABOLIC"        ,"ON"
         ,"FLUXO_PARABOLIC_LIFTING","br1"
         ,"FLUXO_TESTCASE"         ,"default" 
        ]

stat = buildfluxo(buildopts=options,log= True, log_path=outfile,\
                              ntail = args.ntail ,\
                              mpi_procs = args.procs )
if(not stat) :
  builderr= builderr+" 4"
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

outfile="../log.build_test_5"

options=[ "CMAKE_BUILD_TYPE"       ,"Debug"
         ,"FLUXO_BUILD_MPI"        ,MPIOPT
         ,"FLUXO_BUILD_HDF5"       ,HDF5OPT
         ,"FLUXO_EQNSYSNAME"       ,EQNSYS
         ,"FLUXO_DISCTYPE"         ,"2"
         ,"FLUXO_DISC_NODETYPE"    ,"GAUSS"
         ,"FLUXO_PARABOLIC"        ,"OFF"
         ,"FLUXO_TESTCASE"         ,"default" 
        ]

stat = buildfluxo(buildopts=options,log= True, log_path=outfile,\
                              ntail = args.ntail ,\
                              mpi_procs = args.procs )
if(not stat) :
  builderr= builderr+" 4"
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

if(len(builderr) > 1 ) :
  print "!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!"
  print "!!!!! WARNING, following builds failed !!!!!"
  print builderr
  print "see log.build_test files.."
  print "!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!"

