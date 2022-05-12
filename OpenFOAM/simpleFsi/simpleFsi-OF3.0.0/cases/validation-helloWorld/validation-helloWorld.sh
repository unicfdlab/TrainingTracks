#PBS -l walltime=240:00:00,nodes=1:ppn=2

export CASE_DIR=/unicluster/home/matvey.kraposhin/Work/Unicfdlab/github-TrainingTracks/OpenFOAM/simpleFsi-OF3.0.0/cases/validation-helloWorld

cd $CASE_DIR
rm -rf log

decomposePar -cellDist | tee -a log

mpirun -np 2 -npernode 2 --bind-to core --report-bindings -machinefile $PBS_NODEFILE pimpleDyMFoam -parallel | tee -a log



