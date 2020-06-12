module load mpi/mpich/3.1.4/intel/mpivars
module load mpi/mvapich2/2.2a/gcc/mpivars
module load mpi/mvapich2/2.2a/intel/mpivars
module load mpi/openmpi/1.10.0/gcc/mpivars
module load mpi/openmpi/1.6.5/gcc/mpivars
module load mpi/openmpi/1.6.5/intel/mpivars
module load mpi/openmpi/1.8.4/gcc/mpivars
module load mpi/openmpi/1.8.4/intel/mpivars
module load compiler/cuda/8.0/compilervars
module load compiler/mpi/mpich/3.2/gnu
module load apps/lammps/gpu
nvcc -I/usr/mpi/gcc/openmpi-1.4.6/include -L/usr/mpi/gcc/openmpi-1.4.6/lib64 -lmpi main3.cu -o prog