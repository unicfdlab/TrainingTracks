#PBS -l walltime=240:30:00,nodes=1:ppn=6

export CASE_DIR=/unicluster/home/matvey.kraposhin/run/Unicfdlab/github-TrainingTracks1/OpenFOAM/simpleFsi-OF4.1/cases/test-laminar

cd $CASE_DIR
rm -rf log

decomposePar -cellDist | tee -a log

mpirun -np 6 -npernode 6 --bind-to-core --report-bindings -machinefile $PBS_NODEFILE pimpleFoam -parallel | tee -a log



