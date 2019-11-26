#!/bin/bash -l
# Standard output and error:
#SBATCH -o ./tjob.out.%j
#SBATCH -e ./tjob.err.%j
# Initial working directory:
#SBATCH -D ./
# Job Name:
#SBATCH -J convtest_MANSOL
# Queue (Partition):
#SBATCH --partition=express
# Number of nodes and MPI tasks per node:
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=32

# Wall clock limit:
#SBATCH --time=00:30:00

LOGFILE=log.$SLURM_JOB_ID.$SLURM_JOB_NAME

CASEDIR=/u/ipphinde/GIT/fluxo/tests/convergence

cd $CASEDIR

echo '~~~~~~~~~~~~~~~ BEGIN JOBINFO ~~~~~~~~~~~~~~~~~~~~~~~~' | awk '{printf "\n\n %s \n\n\n", $1}' > $LOGFILE
echo "jobid: " $SLURM_JOB_ID "jobname:" $SLURM_JOB_NAME >> $LOGFILE
echo '~~~~~~~~~~~~~~~ END JOBINFO ~~~~~~~~~~~~~~~~~~~~~~~~' | awk '{printf "\n\n %s \n\n\n", $1}' >> $LOGFILE


# Run the program:
#fluxes=("ECKEP" "LF" "HLL" "ECKEP-LF" "ECKEP-Roe" )
#flux_id=("9"    "1"  "22"  "16"       "44" )
#scheme="ECKEP"
#volflux="5"
fluxes=("LF" "HLL" "HLLC" "Roe"  )
flux_id=("1"  "22"    "2"   "3"     )
scheme="stDG"
volflux="0"
fluxo_exec="../../build/bin/fluxo"
nprocs="32"

ii=0
for flux in ${fluxes[@]}
do
  #######
  pname="NS_MANSOL2D_"${scheme}"_"${flux}
  cat parameter_convergence_NS_MANSOL2D.ini |sed -e 's/ProjectName/!ProjectName/g' |sed -e 's/Riemann/!Riemann/g' |sed -e 's/VolumeFlux/!VolumeFlux/g' > ${pname}.ini
  
  echo "ProjectName   = "${pname} >> ${pname}.ini
  echo "Riemann       = "${flux_id[$ii]} >> ${pname}.ini
  echo "VolumeFlux    = "${volflux} >> ${pname}.ini
  python runtest.py -p ${nprocs} ${fluxo_exec} ${pname}.ini >> $LOGFILE
  #tail -n 5 ${pname}.ini

  ii=$((ii+1))
done
