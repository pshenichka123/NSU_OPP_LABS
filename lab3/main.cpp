#include"C:\Program Files (x86)\Microsoft SDKs\MPI\Include\mpi.h"
#define DIM 2
int world_size;

int n[3];
int px;
int py;
int* A_offsets;
int* B_offsets;
int* A_lens;
int* B_lens;
int dims[DIM];
int periods[DIM];
int reorder = 1;
int rank;
int cur_coords[DIM];

double* A_part;
double* B_part;
double* A;
double* B;

void count_ps()
{

}
void count_offsets_lens_A()
{
    A_offsets = new int[py];
    A_lens = new int[py];
    int perProcess = n[1] / py;
    int ost = n[1] % py;
    int start = 0;
    for (int i = 0; i < ost; i++) {
        A_offsets[i] = start;
        A_lens[i] = (perProcess + 1) * n[2];
        start += (perProcess + 1) * n[2];
    }
    for (int i = ost; i < py; i++) {
        A_offsets[i] = start;
        A_lens[i] = perProcess * n[2];
        start += perProcess * n[2];
    }

}



void count_offsets_lens_B()
{
    B_offsets = new int[px];
    B_lens = new int[px];
    int perProcess = n[3] / px;
    int ost = n[3] % px;
    int start = 0;
    for (int i = 0; i < ost; i++) {
        B_offsets[i] = start;
        B_lens[i] = (perProcess + 1) * n[2];
        start += (perProcess + 1) * n[2];
    }
    for (int i = ost; i < px; i++) {
        B_offsets[i] = start;
        B_lens[i] = perProcess * n[2];
        start += perProcess * n[2];
    }

}


void A_B_Init()
{

}



void Scatter_A(MPI_Comm comm2d)
{
    MPI_Group A_group;
    MPI_Comm A_comm;
    MPI_Group All;
    MPI_Comm_group(comm2d, &All);
    int coords[DIM] = { 0,0 };
    int ranks_for_A[py];
    for (int i = 0;i < py;i++)
    {
        coords[1] = i;
        MPI_Cart_rank(comm2d, coords, &(ranks_for_A[i]));
    }
    MPI_Comm_create(comm2d, A_group, &A_comm);
    MPI_Group_incl(All, py, ranks_for_A, &A_group);
    MPI_Scatterv(A, A_lens, A_offsets, MPI_DOUBLE, A_part, A_lens[0], MPI_DOUBLE, 0, comm2d);
    MPI_Bcast(A_part, A_lens[], )



}





int main(int argc, char** argv)
{
    MPI_Init(&argc, &argv);
    MPI_Comm_size(MPI_COMM_WORLD, &world_size);
    count_ps();
    dims[0] = px;
    dims[1] = py;
    count_offsets_lens();
    MPI_Comm comm2d;
    MPI_Cart_create(MPI_COMM_WORLD, DIM, dims, periods, reorder, &comm2d);
    MPI_Comm_rank(comm2d, &rank);
    MPI_Cart_get(comm2d, DIM, dims, periods, cur_coords);

    if (rank == 0) {
        A_B_Init();
    }
    Scatter_A(comm2d);
    Scatter_B();
    MPI_Send();









    MPI_Finalize();
}