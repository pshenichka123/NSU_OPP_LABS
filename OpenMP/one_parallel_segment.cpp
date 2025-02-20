#include <cmath>
#include <iostream>
#include <vector>
#include <cstdio>
#include <cstdlib>
#include <time.h>
static const double eps=0.00001;
static const double tay=0.00001;



void simple_init(double* A, double* b, int N, double* x)
{
for( int i=0;i<N;i++)
{
    for(int j=0;j<N;j++)
    {
        if(i==j)
        {
            A[i*N+j]=2;
        }
        else
        {
            A[i*N+j]=1;
        }
    }
}
for(int i=0;i<N;i++)
{
    b[i]=N+1;
    x[i]=0;
}


}

double* strong_init(double* A, double* b, int N)
{
    for( int i=0;i<N;i++)
    {
        for(int j=0;j<N;j++)
        {
            if(i==j)
            {
                A[i*N+j]=2;
            }
            else
            {
                A[i*N+j]=1;
            }
        }
    }

    double* u=(double*)calloc(N,sizeof(double));
    for( int i=0;i<N;i++)
    {
    u[i]=sin(2*M_PI*(double)i/(double)N);
    }
    for(int i=0;i<N;i++)
    {
        double multiply_current_result=0;
        for(int j=0;j<N;j++)
        {
            multiply_current_result+=A[i*N+j]*u[j];
        }
        b[i]=multiply_current_result;
    }
return u;

}



int main(void)
{

int N;
std::cout<<"mat size=";
std::cin>>N;

double* A=(double*)calloc(N*N,sizeof(double));
double* b=(double*)calloc(N,sizeof(double));
double* x=(double*)calloc(N,sizeof(double));

double* u;

//u=strong_init(A,b,N);
puts("aaa\n");

simple_init(A,b,N,x);
 clock_t start= clock();
while(!Good_Answer(A,b,x,N))
{
    x=Ansewer_iteration(x,A,b,N);
   
///puts("aaa\n");


}

clock_t end= clock();

//for(int i=0;i<N;i++){printf("%f ,  %f\n",x[i],0);}

puts("\n");
printf("time=%d\n", (end-start)/1000);

}
