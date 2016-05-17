#PBS -l walltime=240:30:00,nodes=1:ppn=6

export CASE_DIR=/unicluster/home/matvey.kraposhin/Work/Unicfdlab/github-TrainingTracks/OpenFOAM/simpleFsi-OF3.0.0/cases/validation-laminar

cd $CASE_DIR
rm -rf log

decomposePar -cellDist | tee -a log

mpirun -np 6 -npernode 6 --bind-to core --report-bindings -machinefile $PBS_NODEFILE pimpleDyMFoam -parallel | tee -a log



