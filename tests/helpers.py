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
import shutil
import subprocess

########################################################################################################

def copy2temporary(tmp_dir, f) :
    name = os.path.join(tmp_dir, os.path.basename(f))
    shutil.copy(f, name)
    return name

########################################################################################################
def execute(exec_path, prm_path, projectname, analyze_fcts=None, log = True, ntail = 0 ,\
      mpi_procs = 1, L2 = 1., Linf=1.,PID=0. ) :
   if mpi_procs == 1 :
      cmd = " "
   else :
      cmd = ("mpirun -np %d " % mpi_procs)

   cmd=cmd+" "+exec_path.strip()+" "+prm_path.strip()
   if log :
      #log_path = os.path.splitext(os.path.basename(prm_path))[0] + ".log"
      log_path = "log_"+projectname+".txt"
      plines = open(prm_path, 'r').readlines()
      logf = open(log_path, 'w')
      for pline in plines :
         logf.write(pline)
      logf.close()
   cwd = os.getcwd()  #current working director
   os.system(cmd+" 2>std.err 1>std.out")
   stdout=open("std.out",'r').readlines()
   stderr=open("std.err",'r').readlines()
   if log :
      logf = open(log_path, 'a')
      logf.write('EXECUTE COMMAND:\n')
      logf.write(cmd.strip() + '\n')
      logf.write('EXECUTE: STDOUT FILE:\n')
      for line in stdout :
         logf.write(line)
      logf.write('EXECUTE: STDERR FILE:\n')
      for line in stderr :
         logf.write(line)
      logf.close()
   nlines=len(stdout)
   finishedline = [ line for line in stdout[nlines-10:nlines] if ('finished' in line.lower())]
   if (len(finishedline) == 0) :
      sys.exit(99)  ## be sure that execution marked "crashed"
   if (ntail > 0) :
      for line in stdout[nlines-ntail:nlines] :
        print( line.split("\n")[0] )
   if analyze_fcts :
      results = []
      if type(analyze_fcts) != list : 
         return analyze_fcts(stdout)
      for analyze_fct in analyze_fcts :
         tmp=analyze_fct(stdout)
         results.append(tmp)

      return results

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
            #tmp = tmp[1].split("!")
            #if len(tmp) > 1 :
            #   lines[i] = "%s= %s !%s" % (prop, str(value), tmp[1])
            #else :
            lines[i] = "%s= %s\n" % (prop, str(value))
            found = True
      if (not found) :
         print(  " parameter not found: %s " % key )
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
      print(  " parameter not found: %s " % param )

########################################################################################################

# write summary to screen
def write_summarytable(summary) :
    print(  "\n" )
    print(  "=" * 132 )
    print(  "  S U M M A R Y : " )
    sys.stdout.write("=" * 132)
    header=True
    for line in summary.split('\n') :
      if (header) :
        sys.stdout.write('\n'+line.replace(',',' '))
        sys.stdout.write("\n")
        sys.stdout.write("-" * 132)
        header=False
      else :
        sys.stdout.write('\n'+line.replace(',',' '))
    
    sys.stdout.write("\n")
    print(  "=" * 132 )

########################################################################################################

# extract the L2 error of the last timestep
def get_last_L2_error(lines) :
   for l in lines[-15:] :
      if "L_2" in l :
         tmp = l.split(":")[1]
         return [float(x) for x in tmp.split()]

########################################################################################################

# extract the L_inf error of the last timestep
def get_last_Linf_error(lines) :
   for l in lines[-15:] :
      if "L_inf" in l :
         tmp = l.split(":")[1]
         return [float(x) for x in tmp.split()]

########################################################################################################

def get_last_number(lines) :
   for line in reversed(lines) :
      tmp = line.split(' ')
      for t in reversed(tmp) :
         try :
            return float(t)
         except :
            pass

########################################################################################################

def get_cpu_per_dof(lines) :
   for line in reversed(lines) :
      if "CALCULATION TIME PER TSTEP/DOF: [" in line :
         return float(line.split("[")[1].split("sec")[0])
########################################################################################################

