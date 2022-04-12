#include <mpi.h>
#include <stdlib.h>
#include <stdio.h>

int pes;
int myid;

main(int argc, char **argv)
{
	int i;

	MPI_Init(&argc, &argv);
	MPI_Comm_size(MPI_COMM_WORLD, &pes);
	MPI_Comm_rank(MPI_COMM_WORLD, &myid);
	wd_micro_bench(4096);
	wd_micro_bench(65536);
	wd_micro_bench(1048576);

	MPI_Finalize();
}
