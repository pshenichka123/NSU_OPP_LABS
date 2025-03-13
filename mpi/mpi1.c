


#include <math.h>
#include <stdio.h>
#include<stdlib.h>
#include"C:\Program Files (x86)\Microsoft SDKs\MPI\Include\mpi.h"
#include <stdbool.h>
static double eps=0.00001;
static double tay=0.00001;

void b_init(double *b, int n) {
    for (int i = 0; i < n; i++) {
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


double* Mult_mat_on_vec(double* mat, double* vec,int mat_strings_num,int vec_len,int return_vec_offset)
{   
    double* new_vec=(double*)calloc(vec_len,sizeof(double));
   

    for(int i=0;i<mat_strings_num;i++)
    {
        for(int j=0;j<vec_len;j++)
        {
            new_vec[return_vec_offset+ i]+=mat[i*vec_len+j]*vec[j];
        }
    }
    return new_vec;
}
void sub_vec_from_vec(double* vec1,int offset1, double* vec2,int offset2 , int num_elem_to_suc,double* ret_vec,int offset_ret_vec)
{
    for(int i=0;i<num_elem_to_suc;i++)
    {
        ret_vec[offset_ret_vec+i]=vec1[offset1+i]-vec2[offset2+i];
    }
    
}
void mult_vec_on_num(double* vec, int offset,int num_to_mul,double multiplyer,double* ret_vec)
{
    for(int i=0;i<num_to_mul;i++)
    {
        ret_vec[offset+i]=vec[i+offset]*multiplyer;
    }
}

double* count_x_new(double* A,int mat_col_size,double* x,double* b,int n, int vector_offset)
{



    double* after_mat_mult=Mult_mat_on_vec(A,x,mat_col_size,n,vector_offset);
   
    double* after_b_sub =(double*)calloc(n,sizeof(double));
        sub_vec_from_vec(after_mat_mult,vector_offset,b,vector_offset,mat_col_size,after_b_sub,vector_offset);



    double* after_b_sub_tay=(double*)calloc(n,sizeof(double));
    mult_vec_on_num(after_b_sub,vector_offset,mat_col_size,tay,after_b_sub_tay);
    double* x_new=(double*)calloc(n,sizeof(double));
   


    sub_vec_from_vec(x,vector_offset,after_b_sub_tay,vector_offset,mat_col_size,x_new,0);

    free(after_b_sub_tay);
    free(after_mat_mult);
    return x_new;

}


double* square_vec(double* vec, double* vec_squared ,int elem_num)
{
    for(int i=0;i<elem_num;i++)
    {
        vec_squared[i]=vec[i]*vec[i]; 
    }
    return vec_squared;
}


int *num_elem_per_rank_vector_init(int n, int world_size) {
    int *mas = (int*)calloc(world_size,sizeof(int));
    
    if(world_size>1)
    {
        for (int i = 0; i < world_size - 1; i++) {
            mas[i] = (n / world_size) ;
        }
        mas[world_size-1] = n % (world_size)+(n/world_size);
    }
    else{
        mas[0]=n;
    }
    return mas;
}
int* num_elem_per_rank_martrix_init(int n,int* num_elem_per_rank_vector,int world_size)
{
    int* mass=(int*)calloc(world_size,sizeof(int));

    for(int i = 0;i<world_size;i++)
    
    {
        mass[i]=num_elem_per_rank_vector[i]*n;

    }
return mass;


}
int* count_offset_for_vector(int n, int* num_elem_per_rank_vector,int world_size)
{
    int* mass=(int*)calloc(world_size,sizeof(int));
    mass[0]=0;
    for(int i=1;i<world_size;i++)
    {
        mass[i]=mass[i-1]+num_elem_per_rank_vector[i-1];
    }
    return mass;
}
int* count_offset_for_matrix(int n, int* num_elem_per_rank_matrix,int world_size)
{
    int* mass=(int*)calloc(world_size,sizeof(int));
    mass[0]=0;

    for(int i=1;i<world_size;i++)
    {
        mass[i]=mass[i-1]+num_elem_per_rank_matrix[i-1];
    }
    return mass;
}



double sum_vec_elem(double* vec, int size)
{
    double sum=0;
    for(int i=0;i<size;i++)
    {
        sum+=vec[i];
    }
    return sum;
}




double Part_of_norm_count(double* x,double* A_part, double* b, int vec_len, int mat_str_size )
{
double* after_mat_mul=Mult_mat_on_vec(A_part,x,mat_str_size,vec_len,0);



double* final=calloc(mat_str_size,sizeof(double));
sub_vec_from_vec(after_mat_mul,0,b,0,mat_str_size,final,0);
for(int i=0;i<mat_str_size;i++)
{
    final[i]*=final[i];
}
double sum=sum_vec_elem(final,mat_str_size);
free(after_mat_mul);

free(final);
return sum;
}


int main(int argc, char* argv[]) {

    MPI_Init(&argc, &argv);
    int world_rank;
    MPI_Comm_rank(MPI_COMM_WORLD, &world_rank);
   // printf("Wr=%d\n", world_rank);

    int world_size;
    MPI_Comm_size(MPI_COMM_WORLD, &world_size);
   // printf("WS=%d\n", world_size);

    double *x, *b, *A;
    int n = atoi(argv[1]);
   // printf("%d\n",n);
    int* num_elem_per_rank_vector = num_elem_per_rank_vector_init(n, world_size);
   // for(int i=0;i<world_size;i++){        printf("num_elem_vec=%d\n",num_elem_per_rank_vector[i]);}  
      int* num_elem_per_rank_matrix=num_elem_per_rank_martrix_init(n,num_elem_per_rank_vector,world_size);
    //  for(int i=0;i<world_size;i++){        printf("num_elem_mat=%d\n",num_elem_per_rank_matrix[i]);}   
       int* offset_for_vector=count_offset_for_vector(n,num_elem_per_rank_vector,world_size);
     //  for(int i=0;i<world_size;i++){        printf("offset_vec=%d\n",offset_for_vector[i]);}    
       int* offset_for_matrix= count_offset_for_matrix(n,num_elem_per_rank_matrix,world_size);
       //for(int i=0;i<world_size;i++){printf("offset_mat=%d\n",offset_for_matrix[i]);}

    x = (double*)calloc(n,sizeof(double));
    b = (double*)calloc(n,sizeof(double));
    A = (double*)calloc(n*n,sizeof(double));

    if(world_rank == 0) {
        A_init(A,n);
        
    }
    MPI_Barrier(MPI_COMM_WORLD);

    //for(int i=0;i<n;i++){for(int j=0;j<n;j++){printf("%f ",A[i*n+j]);}printf("\n");}
    b_init(b,n);
    //for(int i=0;i<n;i++){printf("b%f \n", b[i]);}
    double *A_part =(double*)calloc(num_elem_per_rank_matrix[world_rank],sizeof(double)); 

    MPI_Scatterv(A,num_elem_per_rank_matrix ,offset_for_matrix, MPI_DOUBLE,A_part,num_elem_per_rank_matrix[world_rank],MPI_DOUBLE,0,MPI_COMM_WORLD);
    MPI_Barrier(MPI_COMM_WORLD);

   // puts("skat");for(int i=0;i<num_elem_per_rank_vector[world_rank];i++){for(int j=0;j<n;j++){printf("%f ",A_part[i*n+j]);}printf("\n");}
    double* x_new;
    double* x_new_collect=(double*)calloc(n,sizeof(double));
    bool good_ans=false;
    int iter=0;
     while(iter<100){
        
    


        x_new=count_x_new(A_part,num_elem_per_rank_vector[world_rank],x,b,n, offset_for_vector[world_rank]);

        MPI_Barrier(MPI_COMM_WORLD);

        MPI_Gatherv(x_new,num_elem_per_rank_vector[world_rank],MPI_DOUBLE,x_new_collect,num_elem_per_rank_vector, offset_for_vector,MPI_DOUBLE,0,MPI_COMM_WORLD);
        MPI_Barrier(MPI_COMM_WORLD);
        free(x);
        free(x_new);
        x=x_new_collect;
        x_new_collect=(double*)calloc(n,sizeof(double));
        
        MPI_Bcast(x,n,MPI_DOUBLE,0,MPI_COMM_WORLD);
        MPI_Barrier(MPI_COMM_WORLD);

       //for(int i=0;i<num_elem_per_rank_vector[world_rank];i++){printf("new%f ", x_new[i]);}
       
        double part_of_sum=Part_of_norm_count(x,A_part,b,n,num_elem_per_rank_vector[world_rank]);
        double sum=0;
        MPI_Barrier(MPI_COMM_WORLD);
        MPI_Reduce(&part_of_sum,&sum,1,MPI_DOUBLE,MPI_SUM,0,MPI_COMM_WORLD );
        MPI_Barrier(MPI_COMM_WORLD);


        if(world_rank==0)
        {
            double b_norm=(n+1)*sqrt(n);
            if( sqrt(sum)/b_norm<eps )
            {
                    good_ans=true;
            }
    
        }
        MPI_Bcast(&good_ans,1,MPI_C_BOOL,0,MPI_COMM_WORLD);

        if(good_ans==true)
        {
            break;
        }
      
       
        //for(int i=0;i<n;i++){   printf("newxn_after_bcast: %f\n",x[i]); }
           
    }

    if (good_ans)
    {
     puts("win");
    }
    free(num_elem_per_rank_matrix);
    free(num_elem_per_rank_vector);
    free(A);
    free(A_part);
    free(b);
    free(x);
    free(x_new_collect);

    MPI_Finalize();
    return 0;
}

