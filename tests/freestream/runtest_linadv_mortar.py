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
import glob
import sys
import math
import shutil
import tempfile
import time
import argparse
import subprocess
import re

########################################################################################################
# use ../helpers.py
sys.path.append('../')

from helpers import copy2temporary, execute, modify_prm, read_prm
from helpers import get_last_L2_error, get_last_Linf_error, get_cpu_per_dof,write_summarytable

########################################################################################################

# MAIN PROGRAM

########################################################################################################
parser = argparse.ArgumentParser(description='Tool to run fluxo tests')
parser.add_argument('-p','--procs', type=int, default=1, help='number of processors used for the run')
parser.add_argument('-ntail', type=int, default=15, help='number of last line output of screenlog')
parser.add_argument('exe', help='path to executable')
parser.add_argument('prm',   help='path to parameter file')

args = parser.parse_args()

if not os.path.exists(args.prm) :
    print( "parameter-file '%s' not found" % args.prm )
    sys.exit(1)


# copy executable and parameter file to a temporary directory
tmp_dir = tempfile.mkdtemp()
args.exe = copy2temporary(tmp_dir, args.exe)
args.prm = copy2temporary(tmp_dir, args.prm)


#Mortar = ['proj','coll']
Mortar = ['coll']
nMortar=len(Mortar) 
# different polynomial degrees 
Degree = ['1','2','3','4','5','6']
#Degree = ['5','6']
nDegree=len(Degree) 
#Meshes = sorted(glob.glob('../meshes/*_mesh.h5'))
Meshes = [ 'cart_period_1level_1_2_2_elem_mesh.h5'
          ,'cart_period_1level_2_2_2_elem_mesh.h5'
          ,'cart_period_1level_1_2_1_elem_mesh.h5'
          ,'cart_period_1level_2_2_1_elem_mesh.h5'
          ,'cart_period_1level_1_1_2_elem_mesh.h5'
          ,'cart_period_1level_2_1_2_elem_mesh.h5'
          ,'cart_period_1level_1_4_4_elem_mesh.h5'
          ,'cart_period_1level_2_4_4_elem_mesh.h5'
          ,'cart_period_1level_1_4_1_elem_mesh.h5'
          ,'cart_period_1level_2_4_1_elem_mesh.h5'
          ,'cart_period_1level_1_1_4_elem_mesh.h5'
          ,'cart_period_1level_2_1_4_elem_mesh.h5'
          ,'cart_period_2level_1_4_4_elem_mesh.h5'
          ,'cart_period_2level_2_4_4_elem_mesh.h5'
          ,'cart_period_2level_1_4_1_elem_mesh.h5'
          ,'cart_period_2level_2_4_1_elem_mesh.h5'
          ,'cart_period_2level_1_1_4_elem_mesh.h5'
          ,'cart_period_2level_2_1_4_elem_mesh.h5'
         ]

nMeshes=len(Meshes) 
for m in range(0,nMeshes) :
    if( not os.path.exists(Meshes[m])) :
      print( " Mesh %4i %s  does not exist!!" % (m,Meshes[m]))
      sys.exit(1)
    else:
      print( " Mesh %4i %s  " % (m,Meshes[m]))

 
projectname = read_prm(args.prm,'ProjectName')

summaryfilename = '../summary_'+projectname+'.csv'

header=True
# loop over meshes
for r in range(0,nMortar) :
  for m in range(0,nMeshes) :
    for i in range(0,nDegree) :
      print( "               ")
      print( "Mortar: %s , Mesh: %s , Degree: %s  " % (Mortar[r],Meshes[m],Degree[i]))
      print( "               ")
      meshname = re.sub('\_mesh\.h5','',os.path.basename(Meshes[m]))

      projectnameX = projectname+'_mortar_'+Mortar[r]+'_Mesh_'+meshname+'_N_'+Degree[i]
      modify_prm(args.prm, {'ProjectName' : projectnameX})
      print( "               ")
      print( "%3i %3i === > ProjectName: %s" % (i,m,projectnameX))
      print( "               ")
      # modify parameters by replacing string
      #    args.prm = [w.replace('NEX',nElemsX[i] ) for w in args.prm] 
      modify_prm(args.prm, {'N' : Degree[i]})
      modify_prm(args.prm, {'MeshFile' : Meshes[m]})
      if (Mortar[r] == 'proj') :
        modify_prm(args.prm, {'whichMortar' : '0' })
      if (Mortar[r] == 'coll') :
        modify_prm(args.prm, {'whichMortar' : '1' })
      
      # execute fluxo
      start_time = time.time()
      try :
        [L2,Linf,PID] = execute(args.exe, args.prm, projectnameX,\
                                [get_last_L2_error, get_last_Linf_error, get_cpu_per_dof],\
                                log = True, ntail = args.ntail ,\
                                mpi_procs = args.procs )
      except :
      #    shutil.rmtree(tmp_dir)
          print( "\n\n\n    !!! SIMULATION CRASHED?! \n\n")
      #    exit(1)
      end_time = time.time()
      
      nVar=len(L2)
      if(header) :
        summaryheader=( "%5s " % "Mortar" )
        summaryheader=summaryheader+(", %-8s " % " Degree" )
        summaryheader=summaryheader+(", %20s " % "Meshname" )
        for ivar in range(0,nVar) :
          summaryheader=summaryheader+( ", %10s%2i%-9s " % ("   L2(",ivar+1,")") )
        for ivar in range(0,nVar) :
          summaryheader=summaryheader+( ", %10s%2i%-9s " % (" Linf(",ivar+1,")") )
        #summaryheader=summaryheader+( ", %8s   " % "PID" )
        summary=summaryheader
      
        sumfile = open(summaryfilename,'w') #overwrite old file
        sumfile.write(summaryheader)
        sumfile.close()
      
        header=False
      #endif header

    
      summaryline=( " %5s " % Mortar[r] )
      summaryline=summaryline+( ", %20s " % meshname )
      summaryline=summaryline+( ", %8s " % Degree[i] )
      for ivar in range(0,nVar) :
        summaryline=summaryline+(", %21.11e " % (L2[ivar]) )
      for ivar in range(0,nVar) :
        summaryline=summaryline+(", %21.11e " % (Linf[ivar]) )
      
      summary=summary+( "\n" )+summaryline
      
      sumfile = open(summaryfilename,'a')
      sumfile.write('\n'+summaryline)
      sumfile.close()
      
      sys.stdout.flush()
    #end for Mortar
  #end for Meshes
#end for Degree


sumfile.close()
write_summarytable(summary)
print( "table written to %s ..." % summaryfilename)
print( "=" * 132)
 
shutil.rmtree(tmp_dir)
