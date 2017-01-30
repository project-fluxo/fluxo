#!/usr/bin/python
# -*- coding: utf-8 -*-

import os
import shutil
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
   if log :
      #log_path = os.path.splitext(os.path.basename(prm_path))[0] + ".log"
      log_path = "log."+projectname
      plines = open(prm_path, 'r').readlines()
      f = open(log_path, 'w')
      for pline in plines :
         f.write(pline)
      f.close()
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
      f = open(log_path, 'a')
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
   # iterate over all properties, that must be changed and check if
   # one matches the property of the current line
   for key, value in properties.items() :
      found = False
      # iterate over all lines of parameter file
      for i in range(len(lines)) :
         line = lines[i]
         # split line at '='. Before is the property
         if "DEFVAR=" in line.strip() :
            tmp = line.replace("DEFVAR=", "DEFVAR" ,1)
            tmp = tmp.split("=")
            tmp[0] = tmp[0].replace("DEFVAR", "DEFVAR=" ,1)
         else :
            tmp = line.split("=")
         if len(tmp) < 2: continue
         prop = tmp[0]
         if key.strip() == prop.strip() :
            # change property
            tmp = tmp[1].split("!")
            if len(tmp) > 1 :
               lines[i] = "%s= %s !%s" % (prop, str(value), tmp[1])
            else :
               lines[i] = "%s= %s\n" % (prop, str(value))
            found = True
      if (not found) :
         print " parameter not found: %s " % key
   # write parameter file
   f = open(path, 'w')
   for line in lines : 
      f.write(line)
   f.close()

########################################################################################################

def read_prm(path,param) :
   lines = open(path, 'r').readlines()
   found = False
   # iterate over all lines of parameter file
   for line in lines :
      # split line at '='. Before is the property
      if "DEFVAR=" in line.strip() :
         tmp = line.replace("DEFVAR=", 'DEFVAR' ,1)
         tmp = tmp.split("=")
         tmp[0] = tmp[0].replace("DEFVAR", 'DEFVAR=' ,1)
      else :
         tmp = line.split("=")
      if len(tmp) < 2: continue
      if tmp[0].strip() == param :
          # split tmp at '!'. Before is the  value
          tmp = tmp[1].split("!")
          if len(tmp) >= 1 :
             found = True
             return tmp[0].strip()
          else :
             return None
   if (not found) :
      print " parameter not found: %s " % param
