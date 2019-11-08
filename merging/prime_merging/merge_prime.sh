export CCTBX_XFEL=~/ExaFEL/dials.viper
source ${CCTBX_XFEL}/build/setpaths.sh
export PATH=/usr/lib64/mpich/bin:$PATH
mpirun -n 64 prime.mpi_run ls49.phil
