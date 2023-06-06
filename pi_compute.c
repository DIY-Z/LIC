#include "mpi.h"
#include <stdio.h>
/*  https://zhuanlan.zhihu.com/p/399150417  (这个实现方式不是send-receive，而是reduce)*/
double f(double);
double f(double x)
{
    return (4.0/(1.0+x*x));
}

int main(int argc, char *argv[])
{
    int myid, numprocs;

    MPI_Init(&argc, &argv);
    double pi;
    MPI_Comm_size(MPI_COMM_WORLD, &numprocs);
    MPI_Comm_rank(MPI_COMM_WORLD, &myid);
    printf("Process %d of %d\n", myid, numprocs);

    int n = 100;
    double h = 1.0 / (double)n;
    double sum = 0.0;
    for(int i = myid + 1; i <= n; i += numprocs)
    {
        double x = h * ((double)i - 0.5);
        sum += f(x);
    }
    double mypi = h * sum;
    MPI_Reduce(&mypi, &pi, 1, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);
    if (myid == 0)
    {
        printf("The result is %.10f.\n",pi);
    }    

    MPI_Finalize();
}