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

def copy2temporary(tmp_dir, f) :
    name = os.path.join(tmp_dir, os.path.basename(f))
    shutil.copy(f, name)
    return name

########################################################################################################
def execute(exec_path, prm_path, projectname, log = True, ntail = 0 ,\
      mpi_procs = 1) :
   if mpi_procs == 1 :
      cmd = []
   else :
      cmd = ["mpirun", "-np", "%d" % mpi_procs]
   cmd.append(exec_path)
   cmd.append(prm_path)
   #print cmd
   p = subprocess.Popen(cmd, stdout = subprocess.PIPE, stderr = subprocess.PIPE)
   lines = []
   while p.poll() is None :
      line = p.stdout.readlines()
      if line :
         lines.extend(line)
   if p.wait() != 0 :
       for line in lines :
           print line,
       print "!!!!!!!!!!!!!!!!!!!!!!!!!!!"
       print "!!     PROGRAM crashed!    !!"
       print "!!!!!!!!!!!!!!!!!!!!!!!!!!!"
       return None

   #for line in lines :
      #print line,
   if log :
      #log_path = os.path.splitext(os.path.basename(prm_path))[0] + ".log"
      log_path = "log."+projectname
      plines = open(prm_path, 'r').readlines()
      f = open(log_path, 'a')
      for pline in plines :
         f.write(pline)
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

########################################################################################################

def modify_prm(path, properties) :
   lines = open(path, 'r').readlines()
   # iterate over all lines of parameter file
   for i in range(len(lines)) :
      line = lines[i]
      # split line at '='. Before is the property
      tmp = line.split("=")
      if len(tmp) < 2: continue
      prop = tmp[0]
      # iterate over all properties, that must be changed and check if
      # one matches the property of the current line
      for key, value in properties.items() :
         if key == prop.strip() :
            # change property
            tmp = tmp[1].split("!")
            if len(tmp) > 1 :
               lines[i] = "%s= %s !%s" % (prop, str(value), tmp[1])
            else :
               lines[i] = "%s= %s\n" % (prop, str(value))
   # write parameter file
   f = open(path, 'w')
   for line in lines : 
      f.write(line)
   f.close()

########################################################################################################

def read_prm(path,param) :
   lines = open(path, 'r').readlines()
   # iterate over all lines of parameter file
   for line in lines :
      # split line at '='. Before is the property
      tmp = line.split("=")
      if len(tmp) < 2: continue
      if tmp[0].strip() == param :
          # split tmp at '!'. Before is the  value
          tmp = tmp[1].split("!")
          if len(tmp) >= 1 :
              return tmp[0].strip()
          else :
              return None

########################################################################################################

# MAIN PROGRAM

########################################################################################################
parser = argparse.ArgumentParser(description='Tool to generate a series of meshes')
#parser.add_argument('-p','--procs', type=int, default=1, help='number of processors used for the run')
parser.add_argument('-ntail', type=int, default=10, help='number of last line output of screenlog')
parser.add_argument('exe', help='path to executable')
parser.add_argument('prm',   help='path to parameter file')

args = parser.parse_args()

if not os.path.exists(args.prm) :
    print "parameter-file '%s' not found" % args.prm 
    sys.exit(1)


# copy executable and parameter file to a temporary directory
tmp_dir = tempfile.mkdtemp()
args.exe = copy2temporary(tmp_dir, args.exe)
args.prm = copy2temporary(tmp_dir, args.prm)
args.ntail = 10
args.procs = 1


# this generates 3 meshes
nElemsX = ['01','02','04' ]
nElemsY = ['01','02','04' ]
nElemsZ = ['01','02','04' ]

projectname = read_prm(args.prm,'ProjectName')

# loop over meshes
for i in range(0,len(nElemsX)) :

    projectnameX = projectname+'_'+nElemsX[i]+'_'+nElemsY[i]+'_'+nElemsZ[i] 
    modify_prm(args.prm, {'ProjectName' : projectnameX})
    print "               "
    print "%03.0i === > ProjectName: %s" % (i,projectnameX)
    print "               "
    # modify parameters by replacing string
    #    args.prm = [w.replace('NEX',nElemsX[i] ) for w in args.prm] 
    modify_prm(args.prm, {'DEFVAR=(INT):ne_x' : nElemsX[i]})
    modify_prm(args.prm, {'DEFVAR=(INT):ne_y' : nElemsY[i]})
    modify_prm(args.prm, {'DEFVAR=(INT):ne_z' : nElemsZ[i]})


    # execute hopr 
    start_time = time.time()
    try :
        execute(args.exe, args.prm, projectnameX, log = True, ntail = args.ntail ,\
                mpi_procs = args.procs)
    except :
        shutil.rmtree(tmp_dir)
        exit(1)
    end_time = time.time()


    #print end_time - start_time
    sys.stdout.flush()


 
shutil.rmtree(tmp_dir)
