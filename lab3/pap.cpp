#include"C:\Program Files (x86)\Microsoft SDKs\MPI\Include\mpi.h"
#include <cstdio>
int world_size;

int periods[2] = { 0,0 };
int dims[2] = { 4,4 };
int reorder = 1;
int main(int argc, char** argv)
{
    MPI_Init(&argc, &argv);

    MPI_Comm_size(MPI_COMM_WORLD, &world_size);

    dims[0] = 4;
    dims[1] = 4;
    int rank;
    MPI_Comm comm2d;
    MPI_Cart_create(MPI_COMM_WORLD, 2, dims, periods, reorder, &comm2d);
    MPI_Comm_rank(comm2d, &rank);
    int cur_coords[2];
    MPI_Cart_get(comm2d, 2, dims, periods, cur_coords);
    printf("%d,  %d %d", rank, cur_coords[0], cur_coords[1]);
    MPI_Finalize();
    return 0;
}