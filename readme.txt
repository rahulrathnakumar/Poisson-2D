**************************************************************************
THIS FOLDER CONTAINS A PARALLEL IMPLEMENTATION OF THE NUMERICAL SOLUTION OF
THE POISSON EQUATION ON A 2D GRID USING JACOBI ITERATIONS. 
PROGRAMMER: RAHUL RATHNAKUMAR
**************************************************************************
Instructions to run: (Linux Environment Only)
1> Jacobi Iterations using MPI
Using Makefile:
1. Open Ifort using "module load intel-mpi/2018x"
2. Run "make -f MakefileMPI.mak"
Manual:
1. Open Ifort using "module load intel-mpi/2018x"
2. mpiifort -c jacobi.f90 jacobimpi.f90
3. mpiifort jacobi.o jacobimpi.o -o jacobimpi
4. mpirun -n 16 ./jacobimpi

2> Jacobi Iterations using Co-array Fortran
Using bash Script:
./coarray.sh
***************************************************************************