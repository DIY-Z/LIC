//
// Created by 洪锋 on 2022/1/27.
//
//  single process of LIC Visualization

#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include "mpi.h"
#define	 DISCRETE_FILTER_SIZE	1024
#define  LOWPASS_FILTR_LENGTH	8.00000f
#define	 LINE_SQUARE_CLIP_MAX	100000.0f
#define	 VECTOR_COMPONENT_MIN   0.050000f


void     ReadVector(int xres, int yres, float*  pVectr, int* pVectrFlag);
void	 NormalizVectrs(int  n_xres,  int     n_yres,  float*   pVectr);
void     GenBoxFiltrLUT(int  LUTsiz,  float*  p_LUT0,  float*   p_LUT1);
void     MakeWhiteNoise(int  n_xres,  int     n_yres,  unsigned char*  pNoise);
void	 FlowImagingLIC(int  n_xres,  int     n_yres,  float*   pVectr,   unsigned char*  pNoise,
                        unsigned char*  pImage,  float*  p_LUT0,  float*  p_LUT1,  float  krnlen, int start_row, int end_row);
void 	 WriteImage2PPM(int  n_xres,  int     n_yres,  int* pVectrFlag,   unsigned char*  pImage,     char*  f_name);

void saveArrayToFile(unsigned char arr[], int size, const char* prefix, int id) {
    char filename[100];
    snprintf(filename, sizeof(filename), "%s%d.txt", prefix, id);
    FILE* file = fopen(filename, "w");
    if (file == NULL) {
        printf("Failed to open file %s\n", filename);
        return;
    }

    // 写入表头
    fprintf(file, "Index\tValue\n");
    
    // 写入数据和索引
    for (int i = 0; i < size; i++) {
        fprintf(file, "%d\t%u\n", i, arr[i]);
    }

    fclose(file);
    printf("Array saved to file %s\n", filename);
}

int	main(int argc, char* argv[])
{
    clock_t start,end;
    start = clock();
    clock_t start_other, end_other;     //记录执行存储的起始时间和终止时间

    /*
        Note:
            根据所读取的文件中的第一列、第二列的数据来设置left, right, low, high这些值.
    */
    // 根据fielddata.txt发现，left和right应该是指第1列元素的最小值和最大值，low和high指第2列元素的最小值和最大值
    // int left = 110, right = 140, low =12, high = 32;
    float left = -180, right = 179.75, low = -80, high = 79.75;
    float res=0.25;
    int				n_xres = (int)((right-left)/res)+1;         //n_xres相当于经度个数，每个经度之间相隔0.25，经度范围[left,right]
    int				n_yres = (int)((high-low)/res)+1;           //n_yres相当于纬度个数，每个纬度之间相隔0.25，纬度范围[low, high]

    float*			pVectr = (float*         ) malloc( sizeof(float        ) * n_xres * n_yres * 2 );   //pVectr用于存后两列的值,U和V
    int*			pVectrFlag = (int*       ) malloc( sizeof(int          ) * n_xres * n_yres * 2 );   //pVectrFlag用于标记索引在pVectr是否是9999
    float*			p_LUT0 = (float*		 ) malloc( sizeof(float        ) * DISCRETE_FILTER_SIZE);
    float*			p_LUT1 = (float*		 ) malloc( sizeof(float        ) * DISCRETE_FILTER_SIZE);
    unsigned char*	pNoise = (unsigned char* ) malloc( sizeof(unsigned char) * n_xres * n_yres     );
    unsigned char*	pImage = (unsigned char* ) malloc( sizeof(unsigned char) * n_xres * n_yres     );

    ReadVector(n_xres, n_yres, pVectr, pVectrFlag);
    NormalizVectrs(n_xres, n_yres, pVectr);
    MakeWhiteNoise(n_xres, n_yres, pNoise);
    GenBoxFiltrLUT(DISCRETE_FILTER_SIZE, p_LUT0, p_LUT1);

    int myid, numprocs;
    MPI_Init(&argc, &argv);
    MPI_Comm_rank(MPI_COMM_WORLD, &myid);
    MPI_Comm_size(MPI_COMM_WORLD, &numprocs);
    int chunk_size = n_yres / numprocs;             //每个进程需要处理的个数

    // 计算每个进程的起始行和结束行
    int start_row = myid * chunk_size;
    int end_row = start_row + chunk_size;
    printf("\nmyid:%d,start_row:%d,end_row:%d\n",myid, start_row, end_row);
    // 分配每个进程需要的内存
    unsigned char*	pImage_local = (unsigned char* ) malloc( sizeof(unsigned char) * n_xres * n_yres);

    FlowImagingLIC(n_xres, n_yres, pVectr, pNoise, pImage_local, p_LUT0, p_LUT1, LOWPASS_FILTR_LENGTH, start_row, end_row);

    //汇总结果
    if(myid == 0){
        // 将根进程自己的计算结果拷贝到全局结果中
        for (int j = start_row; j < end_row; j++) {
            for (int i = 0; i < n_xres; i++) {
                pImage[j * n_xres + i] = pImage_local[(j - start_row) * n_xres + i];
            }
        }
        printf("myid(%d) copy itself over\n", myid);
        
        // 接收其他进程的计算结果
        for (int i = 1; i < numprocs; i++) {
            MPI_Recv(pImage + i * chunk_size * n_xres, chunk_size * n_xres, MPI_UNSIGNED_CHAR, i, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
            printf("receive over from rank(%d)\n", i);
        }

        const char* prefix = "array_";
        start_other = clock();
        saveArrayToFile(pImage, n_xres * n_yres, prefix, myid);
        end_other = clock();
    }else{
        // const char* prefix = "array_";
        // saveArrayToFile(pImage_local, n_xres * n_yres, prefix, myid);
        MPI_Send(pImage_local + start_row * n_xres, chunk_size * n_xres, MPI_UNSIGNED_CHAR, 0, 0, MPI_COMM_WORLD);
        printf("myid(%d) send over\n", myid);
    }

    WriteImage2PPM(n_xres, n_yres, pVectrFlag, pImage, "LIC-global.ppm");  //这里其实就是将fieldata.txt中的数据可视化成一张ppm图像

    free(pVectr);	pVectr = NULL;
    free(p_LUT0);	p_LUT0 = NULL;
    free(p_LUT1);	p_LUT1 = NULL;
    free(pNoise);	pNoise = NULL;
    free(pImage);	pImage = NULL;

    end = clock();
    if(myid == 0)
        printf("time=%f\n",(double)(end - start - (end_other - start_other))/1000000);

    //释放内存
    free(pImage_local); pImage_local = NULL;
    MPI_Finalize();

    return 0;
}


///		synthesize a saddle-shaped vector field     ///
//void	SyntheszSaddle(int  n_xres,  int  n_yres,  float*  pVectr)
//{
//    for(int  j = 0;  j < n_yres;  j ++)
//        for(int  i = 0;  i < n_xres;  i ++)
//        {
//            int	 index = (  (n_yres - 1 - j) * n_xres + i  )  <<  1;
//            pVectr[index    ] = - ( j / (n_yres - 1.0f) - 0.5f );
//            pVectr[index + 1] =     i / (n_xres - 1.0f) - 0.5f;
//        }
//}


///		read the vector field     ///
void ReadVector(int xres, int yres, float*  pVectr, int* pVectrFlag) {
    FILE *fp;
    if ((fp = fopen("./fielddataglobal.txt", "r")) == NULL) {
        printf("error in reading file !\n");
        exit(1);
    }
    float f1, f2, f3, f4;
    int index = 0;
    while (!feof(fp)) {
        if (fscanf(fp, "%f %f %f %f", &f1, &f2, &f3, &f4) == EOF){   //读到末尾就结束
            printf("Index:%d", index);
            break;
        }
        // printf( "%f %f %f %f \n", f1, f2, f3, f4);
        pVectr[index] = f3;
        if((int)(f3-9999)==0){
            pVectrFlag[index] = 1;
        }
        index++;

        pVectr[index] = f4;
        if((int)(f4-9999)==0){
            pVectrFlag[index] = 1;
        }
        index++;
    }
    fclose(fp);

}


///		normalize the vector field (no bug)    ///
void    NormalizVectrs(int  n_xres,  int  n_yres,  float*  pVectr)
{
    for(int	 j = 0;  j < n_yres;  j ++)
        for(int	 i = 0;  i < n_xres;  i ++)
        {
            int		index = (j * n_xres + i) << 1;      //左移一位，相当于乘2，目的是找到对应位置的索引
            float	vcMag = (float)(  sqrt( (double)(pVectr[index] * pVectr[index] + pVectr[index + 1] * pVectr[index + 1]) )  );

            float	scale = (vcMag == 0.0f) ? 0.0f : 1.0f / vcMag;
            pVectr[index    ] *= scale;
            pVectr[index + 1] *= scale;
        }
}


///		make white noise as the LIC input texture (no bug)    ///
void	MakeWhiteNoise(int  n_xres,  int  n_yres,  unsigned char*  pNoise)
{
    for(int  j = 0;   j < n_yres;  j ++)
        for(int  i = 0;   i < n_xres;  i ++)
        {
            int  r = rand();
            r = (  (r & 0xff) + ( (r & 0xff00) >> 8 )  ) & 0xff;     //0xff指255, 0xff00指65280， 65280 / 255 = 256，而256 = 16 * 16，正好如这个十六进制数所示
            pNoise[j * n_xres + i] = (unsigned char) r;
        }
}


///		generate box filter LUTs (no bug)    ///
void    GenBoxFiltrLUT(int  LUTsiz,  float*  p_LUT0,  float*  p_LUT1)
{
    for(int  i = 0;  i < LUTsiz;  i ++)  p_LUT0[i] = p_LUT1[i] = i;
}


///		write the LIC image to a PPM file     ///
void	WriteImage2PPM(int  n_xres,  int  n_yres,   int* pVectrFlag, unsigned char*  pImage,  char*  f_name)
{
    FILE*	o_file;
    if(   ( o_file = fopen(f_name, "w") )  ==  NULL   )
    {
        printf("Can't open output file\n");
        return;
    }

    fprintf(o_file, "P6\n%d %d\n255\n", n_xres, n_yres);

    for(int  j = n_yres-1;  j > -1 ;  j--)
        for(int  i = 0;  i < n_xres;  i ++)
        {
            unsigned  char	unchar = pImage[j * n_xres + i];
            ///leave the land pixel untouched
            if (pVectrFlag[(j * n_xres + i)*2]  == 1 ){
                unchar = (unsigned char) 255;
                //printf("%d %d \n", i, j);
            }
            //printf("%d %d %d\n", i, j, unchar);
            fprintf(o_file, "%c%c%c", unchar, unchar, unchar);
        }

    fclose (o_file);	o_file = NULL;
}


///		flow imaging (visualization) through Line Integral Convolution (no bug)     ///
void	FlowImagingLIC(int     n_xres,  int     n_yres,  float*  pVectr,  unsigned char*  pNoise,  unsigned char*  pImage,
                       float*  p_LUT0,  float*  p_LUT1,  float   krnlen, int start_row, int end_row)
{
    int		vec_id;						///ID in the VECtor buffer (for the input flow field)
    int		advDir;						///ADVection DIRection (0: positive;  1: negative)
    int		advcts;						///number of ADVeCTion stepS per direction (a step counter)
    int		ADVCTS = (int)(krnlen * 3);	///MAXIMUM number of advection steps per direction to break dead loops

    float	vctr_x;						///x-component  of the VeCToR at the forefront point
    float	vctr_y;						///y-component  of the VeCToR at the forefront point
    float	clp0_x;						///x-coordinate of CLiP point 0 (current)
    float	clp0_y;						///y-coordinate of CLiP point 0	(current)
    float	clp1_x;						///x-coordinate of CLiP point 1 (next   )
    float	clp1_y;						///y-coordinate of CLiP point 1 (next   )
    float	samp_x;						///x-coordinate of the SAMPle in the current pixel
    float	samp_y;						///y-coordinate of the SAMPle in the current pixel
    float	tmpLen;						///TeMPorary LENgth of a trial clipped-segment
    float	segLen;						///SEGment   LENgth
    float	curLen;						///CURrent   LENgth of the streamline
    float	prvLen;						///PReVious  LENgth of the streamline
    float	W_ACUM;						///ACcuMulated Weight from the seed to the current streamline forefront
    float	texVal;						///TEXture VALue
    float	smpWgt;						///WeiGhT of the current SaMPle
    float	t_acum[2];					///two ACcUMulated composite Textures for the two directions, perspectively
    float	w_acum[2];					///two ACcUMulated Weighting values   for the two directions, perspectively
    float*	wgtLUT = NULL;				///WeiGhT Look Up Table pointing to the target filter LUT
    float	len2ID = (DISCRETE_FILTER_SIZE - 1) / krnlen;	///map a curve LENgth TO an ID in the LUT

    ///for each pixel in the 2D output LIC image///
    for(int  j = start_row;	 j < end_row;  j ++)
        for(int  i = 0;	 i < n_xres;  i ++)
        {
            ///init the composite texture accumulators and the weight accumulators///
            t_acum[0] = t_acum[1] = w_acum[0] = w_acum[1] = 0.0f;

            ///for either advection direction///
            for(advDir = 0;  advDir < 2;  advDir ++)
            {
                ///init the step counter, curve-length measurer, and streamline seed///
                advcts = 0;
                curLen = 0.0f;
                clp0_x = i + 0.5f;
                clp0_y = j + 0.5f;

                ///access the target filter LUT///
                wgtLUT = (advDir == 0) ? p_LUT0 : p_LUT1;

                ///until the streamline is advected long enough or a tightly  spiralling center / focus is encountered///
                while( curLen < krnlen && advcts < ADVCTS )
                {
                    ///access the vector at the sample///
                    vec_id = ( (int)(clp0_y) * n_xres + (int)(clp0_x) )<<1;         //目的是找到目标位置的索引
                    vctr_x = pVectr[vec_id    ];
                    vctr_y = pVectr[vec_id + 1];

                    ///in case of a critical point///
                    if( vctr_x == 0.0f && vctr_y == 0.0f )
                    {
                        t_acum[advDir] = (advcts == 0) ? 0.0f : t_acum[advDir];		   ///this line is indeed unnecessary
                        w_acum[advDir] = (advcts == 0) ? 1.0f : w_acum[advDir];
                        break;
                    }

                    ///negate the vector for the backward-advection case///
                    vctr_x = (advDir == 0) ? vctr_x : -vctr_x;
                    vctr_y = (advDir == 0) ? vctr_y : -vctr_y;

                    ///clip the segment against the pixel boundaries --- find the shorter from the two clipped segments///
                    ///replace  all  if-statements  whenever  possible  as  they  might  affect the computational speed///
                    segLen = LINE_SQUARE_CLIP_MAX;
                    segLen = (vctr_x < -VECTOR_COMPONENT_MIN) ? ( (int)(     clp0_x         ) - clp0_x ) / vctr_x : segLen;
                    segLen = (vctr_x >  VECTOR_COMPONENT_MIN) ? ( (int)( (int)(clp0_x) + 1.5f ) - clp0_x ) / vctr_x : segLen;
                    segLen = (vctr_y < -VECTOR_COMPONENT_MIN) ? (((tmpLen = ((int)(clp0_y)-clp0_y)/vctr_y)<segLen)?tmpLen:segLen ): segLen;
                    segLen = (vctr_y >  VECTOR_COMPONENT_MIN) ? (((tmpLen = ((int)((int)(clp0_y) + 1.5f )-clp0_y)/vctr_y)<segLen ) ? tmpLen : segLen): segLen;

                    ///update the curve-length measurers///
                    prvLen = curLen;
                    curLen+= segLen;
                    segLen+= 0.0004f;

                    ///check if the filter has reached either end///
                    segLen = (curLen > krnlen) ? ( (curLen = krnlen) - prvLen ) : segLen;

                    ///obtain the next clip point///
                    clp1_x = clp0_x + vctr_x * segLen;
                    clp1_y = clp0_y + vctr_y * segLen;

                    ///obtain the middle point of the segment as the texture-contributing sample///
                    samp_x = (clp0_x + clp1_x) * 0.5f;
                    samp_y = (clp0_y + clp1_y) * 0.5f;

                    ///obtain the texture value of the sample///
                    texVal = pNoise[ (int)(samp_y) * n_xres + (int)(samp_x) ];

                    ///update the accumulated weight and the accumulated composite texture (texture x weight)///
                    W_ACUM = wgtLUT[ (int)(curLen * len2ID) ];
                    smpWgt = W_ACUM - w_acum[advDir];
                    w_acum[advDir]  = W_ACUM;
                    t_acum[advDir] += texVal * smpWgt;

                    ///update the step counter and the "current" clip point///
                    advcts ++;
                    clp0_x = clp1_x;
                    clp0_y = clp1_y;

                    ///check if the streamline has gone beyond the flow field///
                    if( clp0_x < 0.0f || clp0_x >= n_xres || clp0_y < 0.0f || clp0_y >= n_yres)  break;
                }
            }

            ///normalize the accumulated composite texture///
            texVal = (t_acum[0] + t_acum[1]) / (w_acum[0] + w_acum[1]);

            ///clamp the texture value against the displayable intensity range [0, 255]
            texVal = (texVal <   0.0f) ?   0.0f : texVal;
            texVal = (texVal > 255.0f) ? 255.0f : texVal;
            pImage[j * n_xres + i] = (unsigned char) texVal;
        }

}
