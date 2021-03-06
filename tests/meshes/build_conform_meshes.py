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
import subprocess

########################################################################################################
# use ../helpers.py
sys.path.append('../')

from helpers import copy2temporary, execute, modify_prm, read_prm

########################################################################################################

# MAIN PROGRAM

########################################################################################################
parser = argparse.ArgumentParser(description='Tool to generate a series of meshes')
#parser.add_argument('-p','--procs', type=int, default=1, help='number of processors used for the run')
parser.add_argument('-ntail', type=int, default=30, help='number of last line output of screenlog')
parser.add_argument('exe', help='path to executable')
parser.add_argument('prm',   help='path to parameter file')

args = parser.parse_args()

if not os.path.exists(args.prm) :
    print(  "parameter-file '%s' not found" % args.prm  )
    sys.exit(1)


# copy executable and parameter file to a temporary directory
tmp_dir = tempfile.mkdtemp()
args.exe = copy2temporary(tmp_dir, args.exe)
args.prm = copy2temporary(tmp_dir, args.prm)
args.ntail = 20
args.procs = 1

# this generates 4 3D refinements meshes
nElemsX = ['02','03','04','08' ]
nElemsY = ['02','03','04','08' ]
nElemsZ = ['02','03','04','08' ]

#this generates 5 2D refinements
#nElemsX = ['02','04','08','16' ,'32','64' ]
#nElemsY = ['02','04','08','16' ,'32','64' ]
#nElemsZ = ['01','01','01','01' ,'01','01' ]

projectname = read_prm(args.prm,'ProjectName')

# loop over meshes
for i in range(0,len(nElemsX)) :

    projectnameX = projectname+'_'+nElemsX[i]+'_'+nElemsY[i]+'_'+nElemsZ[i] 
    modify_prm(args.prm, {'ProjectName' : projectnameX})
    print(  "               " )
    print(  "%03.0i === > ProjectName: %s" % (i,projectnameX) )
    print(  "               " )
    # modify parameters by replacing string
    #    args.prm = [w.replace('NEX',nElemsX[i] ) for w in args.prm] 
    modify_prm(args.prm, {'DEFVAR=(INT):ne_x  ' : nElemsX[i]})
    modify_prm(args.prm, {'DEFVAR=(INT):ne_y  ' : nElemsY[i]})
    modify_prm(args.prm, {'DEFVAR=(INT):ne_z  ' : nElemsZ[i]})


    # execute hopr 
    start_time = time.time()
    #try :
    execute(args.exe, args.prm, projectnameX, log = True, ntail = args.ntail ,\
                mpi_procs = args.procs)
    #except :
    #    print(  " crashed... " )
    #    shutil.rmtree(tmp_dir)
    #    exit(1)
    end_time = time.time()


    #print(  end_time - start_time )
    sys.stdout.flush()


 
shutil.rmtree(tmp_dir)
