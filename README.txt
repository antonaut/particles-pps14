build instructions:

	make

run instructions for serial program:

	./linearserial -n <number of particles> -o <output file>

run instructions for pthreads:

	./linpthreads -n <number of particles> -p <number of threads> -o <output file>

run instructions for openMP:

	OMP_NUM_THREADS=<number of threads> ./linopenmp -n <number of particles> -o <output file>

run instructions for MPI:

	mpirun -np <number of processes> ./linmpi -n <number of particles> -o <output file>

visualize results

	cd vis
	make
	cd ../
	./vis/visualize

