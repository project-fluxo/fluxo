!==================================================================================================================================
! Copyright (c) 2016 - 2017 Gregor Gassner
! Copyright (c) 2016 - 2017 Florian Hindenlang
! Copyright (c) 2016 - 2017 Andrew Winters
! Copyright (c) 2010 - 2016 Claus-Dieter Munz (github.com/flexi-framework/flexi)
!
! This file is part of FLUXO (github.com/project-fluxo/fluxo). FLUXO is free software: you can redistribute it and/or modify
! it under the terms of the GNU General Public License as published by the Free Software Foundation, either version 3 
! of the License, or (at your option) any later version.
!
! FLUXO is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY; without even the implied warranty
! of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU General Public License v3.0 for more details.
!
! You should have received a copy of the GNU General Public License along with FLUXO. If not, see <http://www.gnu.org/licenses/>.
!==================================================================================================================================
#include "defines.h"

!==================================================================================================================================
!> Control program of the Fluxo code. Initialization of the computation
!==================================================================================================================================
PROGRAM Fluxo
! MODULES
USE MOD_Globals
USE MOD_Commandline_Arguments
USE MOD_Restart,           ONLY:DefineParametersRestart,InitRestart,Restart,FinalizeRestart
USE MOD_Interpolation,     ONLY:DefineParametersInterpolation,InitInterpolation,FinalizeInterpolation
USE MOD_Mesh,              ONLY:DefineParametersMesh,InitMesh,FinalizeMesh
USE MOD_Mortar,            ONLY:InitMortar,FinalizeMortar
USE MOD_Equation,          ONLY:DefineParametersEquation,InitEquation,FinalizeEquation
USE MOD_IO_HDF5,           ONLY:DefineParametersIO_HDF5,InitIOHDF5
USE MOD_Output,            ONLY:DefineParametersOutput,InitOutput,FinalizeOutput
USE MOD_Analyze,           ONLY:DefineParametersAnalyze,InitAnalyze,FinalizeAnalyze
USE MOD_MPI,               ONLY:DefineParametersMPI,InitMPI
#if MPI
USE MOD_MPI,               ONLY:InitMPIvars,FinalizeMPI
#endif
USE MOD_ReadInTools,       ONLY:prms,IgnoredParameters,PrintDefaultParameterFile,FinalizeParameters
USE MOD_StringTools,       ONLY:STRICMP
USE MOD_TimeDisc,          ONLY:DefineParametersTimedisc,InitTimeDisc,FinalizeTimeDisc,TimeDisc
USE MOD_Testcase,          ONLY:DefineParametersTestcase,InitTestcase,FinalizeTestcase
USE MOD_GetBoundaryFlux,   ONLY:InitBC,FinalizeBC
USE MOD_DG,                ONLY:InitDG,FinalizeDG
#if PARABOLIC
USE MOD_Lifting,           ONLY:DefineParametersLifting,InitLifting,FinalizeLifting
#endif /*PARABOLIC*/
!IMPLICIT NONE
!!----------------------------------------------------------------------------------------------------------------------------------
!! LOCAL VARIABLES
REAL                    :: Time                              !< Used to measure simulation time
!!==================================================================================================================================
CALL InitMPI()
CALL ParseCommandlineArguments()
! Check if the number of arguments is correct
IF ((nArgs.LT.1).OR.(nArgs.GT.2)) THEN
  ! Print out error message containing valid syntax
  CALL CollectiveStop(__STAMP__,'ERROR - Invalid syntax. Please use: fluxo parameter.ini [restart.h5] or fluxo --help'// &
  '[option/section name] to print help for a single parameter, parameter sections or all parameters.')
END IF
CALL DefineParametersMPI()
CALL DefineParametersIO_HDF5()
CALL DefineParametersInterpolation()
CALL DefineParametersRestart()
CALL DefineParametersOutput()
CALL DefineParametersMesh()
CALL DefineParametersEquation()
CALL DefineParametersTestcase()
#if PARABOLIC
CALL DefineParametersLifting ()
#endif /*PARABOLIC*/
CALL DefineParametersTimedisc()
CALL DefineParametersAnalyze()
!
! check for command line argument --help or --markdown
IF (doPrintHelp.GT.0) THEN
  CALL PrintDefaultParameterFile(doPrintHelp.EQ.2, ParameterFile)
  STOP
END IF
CALL prms%read_options(ParameterFile)
!
CALL InitIOHDF5()
IF(MPIroot)THEN
  WRITE(UNIT_stdOut,'(132("="))')
  WRITE(UNIT_stdOut,'(A)') &
  "                  __________________   _______              _______     ______   ______      ______   __________________ "
  WRITE(UNIT_stdOut,'(A)') &
  "                 /                 /) /      /)            /      /)   /     /) /      |   _/     /) /                 /)"
  WRITE(UNIT_stdOut,'(A)') &
  "                /       __________// /      //            /      //   /     // /__     |_/´     _// /      _____      // "
  WRITE(UNIT_stdOut,'(A)') &
  "               /      /)__________) /      //            /      //   /     //  (__|          _/´_) /      /)___/     //  "
  WRITE(UNIT_stdOut,'(A)') &
  "              /      //___         /      //            /      //   /     //      |       _/´_/´  /      //   /     //   "
  WRITE(UNIT_stdOut,'(A)') &
  "             /           /)       /      //            /      //   /     //       |     /´ /´    /      //   /     //    "
  WRITE(UNIT_stdOut,'(A)') &
  "            /      _____//       /      //            /      //   /     //      _/´     |/´     /      //   /     //     "
  WRITE(UNIT_stdOut,'(A)') &
  "           /      /)____)       /      //            /      //   /     //    _/´        |      /      //   /     //      "
  WRITE(UNIT_stdOut,'(A)') &
  "          /      //            /      //_________   /      //___/     // __/´     _     |__   /      //___/     //       "
  WRITE(UNIT_stdOut,'(A)') &
  "         /      //            /                 /) /                 // /      _/´ |      /) /                 //        "
  WRITE(UNIT_stdOut,'(A)') &
  "        /______//            /_________________// /_________________// /_____/` _/´|_____// /_________________//         "
  WRITE(UNIT_stdOut,'(A)') &
  "        )______)             )_________________)  )_________________)  )_____)/´   )_____)  )_________________)          "
  WRITE(UNIT_stdOut,'(A)')
  WRITE(UNIT_stdOut,'(132("="))')
  !write CMAKE compilation options to screen. This file is build during cmake procedure to build/include
#include "configuration-cmake.f90"
  WRITE(UNIT_stdOut,'(132("="))')
END IF !MPIroot
! Measure init duration
StartTime=FLUXOTIME()
!
! Initialization
CALL InitInterpolation()
CALL InitMortar()
CALL InitRestart()
CALL InitOutput()
CALL InitMesh()
#if MPI
CALL InitMPIvars()
#endif
CALL InitEquation()
CALL InitBC()
CALL InitDG()
#if PARABOLIC
CALL InitLifting()
#endif /*PARABOLIC*/
CALL InitTimeDisc()
CALL Restart()
CALL InitAnalyze()
CALL InitTestcase()
! initialization finished
CALL IgnoredParameters()
!
! Measure init duration
Time=FLUXOTIME()
SWRITE(UNIT_stdOut,'(132("="))')
SWRITE(UNIT_stdOut,'(A,F8.2,A)') ' INITIALIZATION DONE! [',Time-StartTime,' sec ]'
SWRITE(UNIT_stdOut,'(132("="))')

! Run Simulation
CALL TimeDisc()

!Finalize
CALL FinalizeOutput()
CALL FinalizeAnalyze()
#if PARABOLIC
CALL FinalizeLifting()
#endif /*PARABOLIC*/
CALL FinalizeDG()
CALL FinalizeEquation()
CALL FinalizeBC()
CALL FinalizeInterpolation()
CALL FinalizeTimeDisc()
CALL FinalizeTestcase()
CALL FinalizeRestart()
CALL FinalizeMesh()
CALL FinalizeMortar()
! Measure simulation duration
Time=FLUXOTIME()
CALL FinalizeParameters()

#if MPI
CALL MPI_FINALIZE(iError)
IF(iError .NE. 0) STOP 'MPI finalize error'
CALL FinalizeMPI()
#endif
SWRITE(UNIT_stdOut,'(132("="))')
SWRITE(UNIT_stdOut,'(A,F8.2,A)') ' FLUXO FINISHED! [',Time-StartTime,' sec ]'
SWRITE(UNIT_stdOut,'(132("="))')
END PROGRAM Fluxo
