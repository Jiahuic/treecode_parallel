This program is built on mvapich2/2.0-gcc-4.9.1, whose treecode part and direct sum part both are parallelized. Users can use different numbers of processors totest the preformance of parallelization.

To run the program:

* <p>On personal computers and/or workstations with the MPICH2 or OpenMPI library, we may use these libraries to lauch MPI-based parallel programs.

The calling syntax of <code>mpiexec</code> is:

<code>mpiexec [mpiexec_options] program_name [program_options]</code>

Example <code> mpiexec -n 4 ./ctreecode.exe</code></p>

* <p>On the SLURM queueing system, submit.job is a batch submission script for that.</p>



