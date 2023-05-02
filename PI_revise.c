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
    if(myid == root)
    {   // 在根进程中，从标准输入中读取n的值，并初始化startTime变量，以便记录程序执行时间。
        printf("Enter N:");
        scanf("%d", &n);
        startTime = MPI_Wtime();
    }
    MPI_Bcast(&n, 1, MPI_INT, root, MPI_COMM_WORLD);
    h = 1.0 / (double)n;
    for(i = myid; i < n; i += numprocs)  //在MPI_Bcast之后使用循环计算每个进程负责的区间，然后将结果累加到res中
    {
        xi = h * ((double)i + 0.5);
        res += func(xi);
    }
    res *= h;
    // 使用MPI_Reduce操作将各进程计算的结果求和，然后将结果发送回根进程。
    MPI_Reduce(&res, &pi, 1, MPI_DOUBLE, MPI_SUM, root, MPI_COMM_WORLD);
    if(myid == root)
    {
        endTime = MPI_Wtime();
        printf("\nPI is %f\nTime is : %f\n", pi, endTime - startTime);
    }
    MPI_Finalize();
    return 0;
}
