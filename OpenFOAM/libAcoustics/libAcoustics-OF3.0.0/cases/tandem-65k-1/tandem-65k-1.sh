#PBS -l walltime=8192:00:00,nodes=1:ppn=4
#PBS -q workq

cd /unicluster/home/matvey.kraposhin/Work/Unicfdlab/github-TrainingTracks1/OpenFOAM/libAcoustics-OF3.0.0/cases/tandem-65k-1

#mpirun -np 12 -machinefile $PBS_NODEFILE pisoFoam -parallel | tee -a std.log
rm -rf std.log
decomposePar -force | tee -a std.log
mpirun -np 4 -machinefile $PBS_NODEFILE pisoFoam -parallel | tee -a std.log


