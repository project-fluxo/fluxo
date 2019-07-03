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
sys.path.append('../../../tests/')

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

case = ['EC','EC']


Degree= ['5','4']

riemann = ['11','11']

Meshes=[ '../ConformPeriodicDeform3D_Ngeo5_007_007_007_mesh.h5'
        ,'../ConformPeriodicDeform3D_Ngeo4_007_007_007_mesh.h5' ]

ncases=len(case)

CFL =    ['0.0125'] #['0.8' ,'0.4' , '0.2' ,'0.1' ,'0.05','0.025']
CFLstr = ['0012']   #['0800','0400', '0200','0100','0050','0025']
nCFL=len(CFL) 

 
projectname = read_prm(args.prm,'ProjectName')

summaryfilename = '../summary_'+projectname+'.csv'

header=True
# loop over cases and CFL
for i in range(0,ncases) :
  for m in range(0,nCFL) :
    print( "               ")
    print( "case: %s , CFL: %s " % (case[i],CFL[m]))
    print( "               ")
    projectnameX = projectname+'_'+case[i]+'_N'+Degree[i]+'_dt'+CFLstr[m]
    modify_prm(args.prm, {'ProjectName' : projectnameX})
    print( "               ")
    print( "%3i %3i === > ProjectName: %s" % (i,m,projectnameX))
    print( "               ")
    # modify parameters by replacing string
    #    args.prm = [w.replace('NEX',nElemsX[i] ) for w in args.prm] 
    modify_prm(args.prm, {'N' : Degree[i]})
    modify_prm(args.prm, {'Riemann' : riemann[i]})
    modify_prm(args.prm, {'MeshFile' : Meshes[i]})
    modify_prm(args.prm, {'CFLscale' : CFL[m]})

    # execute fluxo
    #start_time = time.time()
    try :
      [L2,Linf,PID] = execute(args.exe, args.prm, projectnameX,\
                              [get_last_L2_error, get_last_Linf_error, get_cpu_per_dof],\
                              log = True, ntail = args.ntail ,\
                              mpi_procs = args.procs )
    except :
      print( " crashed ")
    #    shutil.rmtree(tmp_dir)
    #    exit(1)
    #end_time = time.time()

  #end for CFL
#end for case


print( "=" * 132)
 
shutil.rmtree(tmp_dir)
