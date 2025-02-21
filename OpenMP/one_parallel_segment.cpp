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
 bool is_good_answer=false;
 puts("after init\n");

 double* x_new;
 double* top;
double bot_module;
    double top_module;
     while (!is_good_answer)  
        {
        #pragma omp parallel
            {
          
            {
          #pragma omp single    
            x_new=(double*)calloc(N,sizeof(double));
            

        #pragma omp for
            for(int i=0;i<N;i++)
            {
                double correction=0;
                for(int j=0;j<N;j++)
                {
                    correction+=A[i*N+j]*x[j];
                }
                x_new[i]=x[i]-tay*(correction-b[i]);
            }
           
           #pragma omp single
            {
            free(x);
            x=x_new;
        
            top=(double*)calloc(N,sizeof(double));
            top_module=0;
            }

            #pragma omp for
            for(int i=0;i<N;i++)
            {
                top[i]=0;
                for(int j=0;j<N;j++)
                {
                    top[i]+=A[i*N+j]*x[j];
                }
                top[i]-=b[i];
            }
            
            #pragma omp for
            for(int i=0;i<N;i++)
            {
                top_module+=top[i]*top[i];
            }
            #pragma omp single
            {
            top_module=sqrt(top_module);
            free(top);
            bot_module=0;
            }
            
           #pragma omp for
            for(int i=0;i<N;i++)
            {
                bot_module+=b[i]*b[i];
            }
           
              {
                bot_module=sqrt(bot_module);
                if(top_module/bot_module<eps)
                {
                
                    is_good_answer=true;
                }
                else
                {
                    //printf("\n%f\n", top_module/bot_module);
                
                    is_good_answer=false;
                } 
              }            
            }
        }
        }
     
 
 
        puts("aaa\n");

 


clock_t end= clock();

    puts("\n");
    printf("time=%d\n", (end-start));
   // for(int i=0;i<N;i++){printf("%f ,  %f\n",x[i],0);}

    return 0;
}
