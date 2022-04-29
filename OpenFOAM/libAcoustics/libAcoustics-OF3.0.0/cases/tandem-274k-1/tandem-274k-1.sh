#PBS -l walltime=8192:00:00,nodes=1:ppn=12

cd /unicluster/home/matvey.kraposhin/Work/Unicfdlab/github-TrainingTracks/OpenFOAM/libAcoustics-OF3.0.0/cases/tandem-274k-1

rm -rf std.log
mpirun -np 12 -machinefile $PBS_NODEFILE pisoFoam -parallel | tee -a std.log



