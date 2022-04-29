#PBS -l walltime=240:00:00,nodes=1:ppn=2

export CASE_DIR=/unicluster/home/matvey.kraposhin/Work/Unicfdlab/github-TrainingTracks1/OpenFOAM/simpleFsi-OF4.1/cases/test-helloWorld

cd $CASE_DIR
rm -rf log

decomposePar -cellDist | tee -a log

mpirun -np 2 -npernode 2 --bind-to-core --report-bindings -machinefile $PBS_NODEFILE pimpleFoam -parallel | tee -a log



