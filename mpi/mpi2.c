
#include <math.h>
#include <stdio.h>
#include<stdlib.h>
#include"C:\Program Files (x86)\Microsoft SDKs\MPI\Include\mpi.h"
#include <stdbool.h>
static double eps=0.00001;
static double tay=0.00001;



void b_init(double *b, int len,int n) {
    for (int i = 0; i < len; i++) {
        b[i] = n + 1;
    }
}
void A_init(double *A, int n) {
    for (int i = 0; i < n; i++) {
        for (int j = 0; j < n; j++) {
            if (i == j) {
                A[i*n + j] = 2;
            }
            else {
                A[i*n + j] = 1;
            }
        }
    }
}


double count_occuracy(double* x, int x_len)
{
    double sum=0;
    for(int i=0;i<x_len;i++)
    {
        sum+=x[i]-1;
    }
    return sum;
}

double* mult_mat_on_vec(double* mat,int n, double* vector, int vec_len, int vec_offset_from_mat_string_beginning,int mat_string_num,double* ret_vec)
{

    for(int i=0;i<mat_string_num;i++)
    {

        for(int j=0;j<vec_len;j++)
        {
            ret_vec[i]+=mat[i*n+vec_offset_from_mat_string_beginning+j]*vector[j];
        }
    }
    return ret_vec;

}


void sub_vec_from_vec(double* a, double* b,double* ret_vec, int len)
{

    for(int i=0;i<len;i++)
    {
        ret_vec[i]=a[i]-b[i];
    }
}
void mult_vec_num(double* vec,double* ret_vec, int len, double mult)
{
    for(int i=0;i<len;i++)
    {
        ret_vec[i]=vec[i]*mult;
     //   printf("retvec%f  ",ret_vec[i]);

    }

}
double* count_new_x_part(double* A_part,double* b_part,double* x_part, int cur_proc_vector_offset,int cur_proc_element_in_vector,int n, int world_size,int world_ramk,int num_of_lines_in_matrix)
{
    double* x_part_new=calloc(cur_proc_element_in_vector,sizeof(double));
    int cur_iter_vector_len=cur_proc_element_in_vector;
    int cur_iter_vector_offset=cur_proc_vector_offset;
    for(int i=0;i<world_size;i++)
    {
    mult_mat_on_vec(A_part,n,x_part,cur_iter_vector_len,cur_iter_vector_offset,num_of_lines_in_matrix,x_part_new);
    MPI_Sendrecv_replace(x_part,cur_iter_vector_len,MPI_DOUBLE,(world_ramk + 1) % world_size,123, (world_ramk - 1 + world_size) % world_size,123,MPI_COMM_WORLD,MPI_STATUS_IGNORE);
    MPI_Sendrecv_replace(&cur_iter_vector_offset,1,MPI_INT, (world_ramk + 1) % world_size,423, (world_ramk - 1 + world_size) % world_size,423,MPI_COMM_WORLD,MPI_STATUS_IGNORE);
//    printf("offset:%d in world %d  ",cur_iter_vector_offset, world_ramk);
//    puts("iter");
}

    double* after_b_sub=calloc(cur_proc_element_in_vector,sizeof(double));
    sub_vec_from_vec(x_part_new,b_part,after_b_sub,cur_proc_element_in_vector);
//    for(int i=0;i<cur_proc_element_in_vector;i++)
//    {
      //  printf("afterb%f  ",after_b_sub[i]);
//    }
    double*Ax_b_tay=calloc(cur_proc_element_in_vector,sizeof(double));
    mult_vec_num(after_b_sub,Ax_b_tay,cur_proc_element_in_vector,tay);
//    printf("len=%d  ", cur_proc_element_in_vector);
//    for(int i=0;i<cur_proc_element_in_vector;i++)
//    {
     //   printf("aftertay%f  ",Ax_b_tay[i]);
//
//    }
    double* x_return=calloc(cur_proc_element_in_vector,sizeof(double));
    sub_vec_from_vec(x_part,Ax_b_tay,x_return,cur_proc_element_in_vector);
    free(x_part_new);
    free(after_b_sub);
    free(Ax_b_tay);
    return x_return;

}



double calc_top_module(double* x_part, double* A_part,  double* b_part, int cur_proc_vec_len,int world_rank, int world_size, int cur_proc_vector_offset, int n)
{
    double* x_part_new=calloc(cur_proc_vec_len,sizeof(double));
    int cur_iter_vector_len=cur_proc_vec_len;
    int cur_iter_vector_offset=cur_proc_vector_offset;
    for(int i=0;i<world_size;i++)
    {
    mult_mat_on_vec(A_part,n,x_part,cur_iter_vector_len,cur_iter_vector_offset,cur_iter_vector_len,x_part_new);
    MPI_Sendrecv_replace(x_part,cur_iter_vector_len,MPI_DOUBLE,(world_rank + 1) % world_size,223, (world_rank - 1 + world_size) % world_size,223,MPI_COMM_WORLD,MPI_STATUS_IGNORE);
    MPI_Sendrecv_replace(&cur_iter_vector_offset,1,MPI_INT, (world_rank + 1) % world_size,323, (world_rank - 1 + world_size) % world_size,323,MPI_COMM_WORLD,MPI_STATUS_IGNORE);
//    puts("iter");
    }

    double* after_b_sub=calloc(cur_proc_vec_len,sizeof(double));
    sub_vec_from_vec(x_part_new,b_part,after_b_sub,cur_proc_vec_len);
    double  sum=0;
    for(int i=0;i<cur_proc_vec_len;i++)
    {
        sum+=after_b_sub[i]*after_b_sub[i];
//        printf("  sum%f  \n",after_b_sub[i]);
    }
    free(x_part_new);
    free(after_b_sub);
    return sum;
}




int main(int argc, char* argv[])
{

    int n= atoi(argv[1]);
    MPI_Init(&argc, &argv);
    int world_rank;
    MPI_Comm_rank(MPI_COMM_WORLD, &world_rank);
    int world_size;
    MPI_Comm_size(MPI_COMM_WORLD, &world_size);

    int elements_in_vector_per_rank =n/world_size;
    int offsets_in_vector_per_rank = elements_in_vector_per_rank*world_rank;
    int number_of_lines_in_matrix=n/world_size;
    int num_elem_per_rank_matrix=number_of_lines_in_matrix*n;
    double* A=calloc(n*n,sizeof(double));
    A_init(A,n);

    double* A_part= calloc(num_elem_per_rank_matrix,sizeof(double));
    MPI_Scatter(
        A,
        num_elem_per_rank_matrix,
        MPI_DOUBLE,
        A_part,
        num_elem_per_rank_matrix,
        MPI_DOUBLE,
        0,
        MPI_COMM_WORLD
    );

    double* b_part=calloc(elements_in_vector_per_rank,sizeof(double));
    b_init(b_part,elements_in_vector_per_rank,n);

    double*x_part=calloc(elements_in_vector_per_rank,sizeof(double));

    bool is_good_ans=false;
    
    MPI_Barrier(MPI_COMM_WORLD);
    double start_time = MPI_Wtime();

    while(true )
    {

        double* x_part_new=count_new_x_part(
            A_part,
            b_part,
            x_part,
            offsets_in_vector_per_rank,
            elements_in_vector_per_rank,
            n,
            world_size,
            world_rank,
            number_of_lines_in_matrix
        );

      // for(int i=0;i<elements_in_vector_per_rank;i++){    printf("%f     ",x_part_new[i]); }

        free(x_part);
        x_part=x_part_new;
        double top_module=0;
        double top_module_part=calc_top_module(x_part,
            A_part,
            b_part,
            elements_in_vector_per_rank,
            world_rank,
            world_size,
            offsets_in_vector_per_rank,
            n);
//        printf("part=%f",top_module_part);
        MPI_Reduce(&top_module_part,&top_module,1,MPI_DOUBLE,MPI_SUM,0,MPI_COMM_WORLD);
        if(world_rank==0)
        {
            if(top_module/((n+1)*(n+1))<eps)
            {
                is_good_ans=true;


            }
//            printf("drob=%f\n",top_module/((n+1)*(n+1)));

        }
    //    printf( "ans is:%d       \n",is_good_ans);
        MPI_Bcast(&is_good_ans,1,MPI_C_BOOL,0,MPI_COMM_WORLD);
        if(is_good_ans)
        {
            break;
        }
      

    }

 MPI_Barrier(MPI_COMM_WORLD);
 double end_time = MPI_Wtime();

    double part_occur=count_occuracy(x_part,elements_in_vector_per_rank);
    double occur=0;
    MPI_Reduce(&part_occur,&occur,1,MPI_DOUBLE,MPI_SUM,0,MPI_COMM_WORLD);
    if(world_rank==0)
    {
        printf("occur=%f\n",occur);
        printf("time=%f\n",end_time-start_time);
    }



    free(A);
    free(A_part);
    free(b_part);
    free(x_part);


        MPI_Finalize();
        return 0;

}