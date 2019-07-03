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


# different polynomial degrees 
Degree = ['4','5']
nDegree=len(Degree) 
#Meshes = sorted(glob.glob('../meshes/*_mesh.h5'))
Meshes = [ '../meshes/CartBoxPeriodic_02_02_02_mesh.h5'
          ,'../meshes/CartBoxPeriodic_04_04_04_mesh.h5'
#          ,'../meshes/CartBoxPeriodic_08_08_08_mesh.h5'
#          ,'../meshes/CartBoxPeriodic_16_16_16_mesh.h5'
#---
          ,'../meshes/ConformBoxTrilinear_02_02_02_mesh.h5'
          ,'../meshes/ConformBoxTrilinear_04_04_04_mesh.h5'
#          ,'../meshes/ConformBoxTrilinear_08_08_08_mesh.h5'
#          ,'../meshes/ConformBoxTrilinear_16_16_16_mesh.h5'
#---
          ,'../meshes/ConformBoxCurved_Ngeo3_02_02_02_mesh.h5'
          ,'../meshes/ConformBoxCurved_Ngeo3_04_04_04_mesh.h5'
#          ,'../meshes/ConformBoxCurved_Ngeo3_08_08_08_mesh.h5'
#          ,'../meshes/ConformBoxCurved_Ngeo3_16_16_16_mesh.h5'
#---
          ,'../meshes/DeformedBoxMortar_Ngeo_1_Level_01_mesh.h5'
          ,'../meshes/DeformedBoxMortar_Ngeo_1_Level_02_mesh.h5'
#          ,'../meshes/DeformedBoxMortar_Ngeo_1_Level_04_mesh.h5'
#          ,'../meshes/DeformedBoxMortar_Ngeo_1_Level_08_mesh.h5'
#---
          ,'../meshes/DeformedBoxMortar_Ngeo_2_Level_01_mesh.h5'
          ,'../meshes/DeformedBoxMortar_Ngeo_2_Level_02_mesh.h5'
#          ,'../meshes/DeformedBoxMortar_Ngeo_2_Level_04_mesh.h5'
#          ,'../meshes/DeformedBoxMortar_Ngeo_2_Level_08_mesh.h5'
#---
          ,'../meshes/DeformedBoxMortarPeriodic_Ngeo_1_Level_01_mesh.h5'
          ,'../meshes/DeformedBoxMortarPeriodic_Ngeo_1_Level_02_mesh.h5'
#          ,'../meshes/DeformedBoxMortarPeriodic_Ngeo_1_Level_04_mesh.h5'
#          ,'../meshes/DeformedBoxMortarPeriodic_Ngeo_1_Level_08_mesh.h5'
#---
          ,'../meshes/DeformedBoxMortarPeriodic_Ngeo_2_Level_01_mesh.h5'
          ,'../meshes/DeformedBoxMortarPeriodic_Ngeo_2_Level_02_mesh.h5'
#          ,'../meshes/DeformedBoxMortarPeriodic_Ngeo_2_Level_04_mesh.h5'
#          ,'../meshes/DeformedBoxMortarPeriodic_Ngeo_2_Level_08_mesh.h5'
         ]

nMeshes=len(Meshes) 
for m in range(0,nMeshes) :
    print( " Mesh %4i %s  " % (m,Meshes[m]))

if (nMeshes == 0 ) :
    print( " NO MESHES FOUND IN ../meshes  ")
    shutil.rmtree(tmp_dir)
    exit(1)
else :
    print( " %4i MESH(ES) FOUND IN ../meshes  " % nMeshes)
 
projectname = read_prm(args.prm,'ProjectName')

summaryfilename = '../summary_'+projectname+'.csv'

header=True
# loop over meshes
for i in range(0,nDegree) :
  for m in range(0,nMeshes) :
    print( "               ")
    print( "Degree: %s , Mesh: %s " % (Degree[i],Meshes[m]))
    print( "               ")
    meshname = re.sub('\_mesh\.h5','',os.path.basename(Meshes[m]))

    projectnameX = projectname+'_Degree_'+Degree[i]+'_Mesh_'+meshname 
    modify_prm(args.prm, {'ProjectName' : projectnameX})
    print( "               ")
    print( "%3i %3i === > ProjectName: %s" % (i,m,projectnameX))
    print( "               ")
    # modify parameters by replacing string
    #    args.prm = [w.replace('NEX',nElemsX[i] ) for w in args.prm] 
    modify_prm(args.prm, {'N' : Degree[i]})
    modify_prm(args.prm, {'MeshFile' : Meshes[m]})

    # execute fluxo
    start_time = time.time()
    try :
      [L2,Linf,PID] = execute(args.exe, args.prm, projectnameX,\
                              [get_last_L2_error, get_last_Linf_error, get_cpu_per_dof],\
                              log = True, ntail = args.ntail ,\
                              mpi_procs = args.procs )
    except :
        shutil.rmtree(tmp_dir)
        exit(1)
    end_time = time.time()

    nVar=len(L2)
    if(header) :
      summaryheader=("%-8s " % " Degree" )
      summaryheader=summaryheader+( ", %50s " % "Meshname" )
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

    
    summaryline=( "%8s " % Degree[i] )
    summaryline=summaryline+( ", %45s " % meshname )
    for ivar in range(0,nVar) :
      summaryline=summaryline+(", %21.11e " % (L2[ivar]) )
    for ivar in range(0,nVar) :
      summaryline=summaryline+(", %21.11e " % (Linf[ivar]) )

    summary=summary+( "\n" )+summaryline

    sumfile = open(summaryfilename,'a')
    sumfile.write('\n'+summaryline)
    sumfile.close()

    sys.stdout.flush()
  #end for Meshes
#end for Degree


sumfile.close()
write_summarytable(summary)
print( "table written to %s ..." % summaryfilename)
print( "=" * 132)
 
shutil.rmtree(tmp_dir)
