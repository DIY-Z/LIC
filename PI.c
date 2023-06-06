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
    if(myid != 0)
    {
        MPI_Bcast(&n, 1, MPI_INT, root, MPI_COMM_WORLD);  //从非根节点调用该函数，变量n会被赋值成从根节点接受到的数据
        h = 1.0 / (double)n;
        for(int i = myid; i <= n; i += (numprocs - 1))
        {
            xi = h * ((double)i - 0.5);
            res += func(xi);
            // MPI_Send(&res, 1, MPI_DOUBLE, 0, 99, MPI_COMM_WORLD);
        }
        res = h * res;
        MPI_Send(&res, 1, MPI_DOUBLE, 0, 99, MPI_COMM_WORLD);  //都算完了再发
    }
    else
    {
        printf("Enter N:");   //N相当于将0到1之间等分出来的矩形的个数(https://zhuanlan.zhihu.com/p/399150417)
        scanf("%d", &n);
        startTime = MPI_Wtime();
        MPI_Bcast(&n, 1, MPI_INT, root, MPI_COMM_WORLD); //从根节点调用该函数, 变量n里面的值会被发送到其它的节点上
        h = 1.0 / (double)n;
        // for(i = 1; i < n; i += (numprocs - 1))
        // {
        //     MPI_Recv(&res, 1, MPI_DOUBLE, i, 99, MPI_COMM_WORLD, &status);
        //     pi = pi + res;
        // }
        for(i = 1; i < numprocs; i++)
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

//主要修改了两处地方
//(1)首先将非根进程的那个循环里的MPI_Send移出来，移到循环外面，这样就变成了每个进程将自己的那部分算完了并求和完了，这时候再send给根进程
//(2)然后就是把根进程的那个循环去掉，不应该循环n，而是循环进程数，需要从每个进程那里receive, 然后再求和得到最后的结果
//(3)此外还有一些小错，之前就修改了，比如"某些变量事先没有定义就使用"等小问题