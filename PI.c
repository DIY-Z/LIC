#include<stdio.h>
#include<math.h>
#include "mpi.h"

double func(double xi)
{
    return (4.0 / (1.0 + xi * xi));
}

int main(int argc, char* argv[])
{
    int n, myid, numprocs, i;
    double pi, h, xi, res, startTime, endTime;
    pi = 0.0;
    res = 0.0;
    int root = 0;
    MPI_Init(&argc, &argv);
    MPI_Comm_rank(MPI_COMM_WORLD, &myid);
    MPI_Comm_size(MPI_COMM_WORLD, &numprocs);
    MPI_Status status;
    if(myid != root)
    {
        MPI_Bcast(&n, 1, MPI_INT, root, MPI_COMM_WORLD);
        h = 1.0 / (double)n;
        for(int i = myid; i <= n; i += (numprocs - 1))
        {
            xi = h * ((double)i - 0.5);
            res += func(xi);
            MPI_Send(&res, 1, MPI_DOUBLE, 0, 99, MPI_COMM_WORLD);
        }
        res = h * res;
    }
    else
    {
        printf("Enter N:");
        scanf("%d", &n);
        startTime = MPI_Wtime();
        MPI_Bcast(&n, 1, MPI_INT, root, MPI_COMM_WORLD);
        h = 1.0 / (double)n;
        for(i = 1; i < n; i += (numprocs - 1))
        {
            MPI_Recv(&res, 1, MPI_DOUBLE, i, 99, MPI_COMM_WORLD, &status);
            pi = pi + res;
        }
        endTime = MPI_Wtime();
        printf("\nPI is %f\nTime is : %f\n", pi, endTime - startTime);
    }
    MPI_Finalize();
    return 0;
}

/*这是一段存在bug的计算PI的代码，请改正它正确计算PI*/