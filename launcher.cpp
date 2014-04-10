#include "build_value.h"
#include "io_utils.h"
#include "video_encode.h"
#include "data_alloc.h"
#include "algs_direct.h"
#include "algs_intra.h"
#include "algs_inter.h"
#include "algs_inter_intra.h"
#include "utils.h"
#include "uep_predict_utils.h"

#define DECODE __ALGO__

const double DERIVE_SNR[5][PXL]={{2.597844, 2.244942, 1.760056, 1.214787, 0.130591, 0.017250, 0.017250, 0.017250},
                                {2.382862, 1.963794, 1.582247, 1.234872, 0.765783, 0.043096, 0.013702, 0.013702},
                                {2.092984, 1.818853, 1.521252, 1.218653, 0.921529, 0.402081, 0.013660, 0.010884},
                                {1.933909, 1.607395, 1.418209, 1.193889, 0.953517, 0.706520, 0.176923, 0.009999},
                                {1.720262, 1.532340, 1.260163, 1.118219, 0.943704, 0.752738, 0.552381, 0.121300}};

const double COMP_SNR[4][PXL] = {{2.6493 , 2.233, 1.7531, 1.2123,  0.1237 , 0.01, 0.01, 0.0086} ,
                                 {2,2,1.75, 1.21 ,0.26, 0.26 , 0.26 ,0.26 },
                                 {1.5, 1.5, 1.25, 1.25, 0.75, 0.75, 0.5, 0.5},
                                 {1, 1, 1, 1, 1, 1, 1, 1}};
const double DERIVE_SNR_N1[PXL] = {3.015249, 2.409902, 1.813491, 0.728476, 0.032169, 0.000001, 0.000001, 0.000001 };

const double TEST[PXL] = {1.977668, 1.538618, 0.943222, 0.049886, 0.01, 0.01, 0.01, 0.01};
const double TEST_TEST[PXL]= {3.476999, 2.704871, 1.657564, 0.087482, 0.018022, 0.018022, 0.018022, 0.018022};

const double MINPOW_PSNR25[PXL] = {2.357380, 1.877585, 1.364773, 0.268487, 0.019715, 0.013750, 0.013750, 0.013750};
int main(int argc,char* argv[]){

    startRandom();

    const int h = __HEIGHT, w = __WIDTH, f = __FRAME ;
    const int lm = h*w;
    const int lu = lm+(G_L-1);
//    double snr = __SNR;
    const int len = __SNR_E-__SNR_S+1;
    double* snr = (double*)malloc(sizeof(double)*len);

    int ** G = getGenerator();
    double* PSNR = (double*) malloc(sizeof(double)*f);
    double *** Ly = new3d<double>(f,PXL,2*lu,0);
//    int *** x = new3d<int>(f,PXL,2*lu);
    int *** map_out = new3d<int>(f,PXL,lm);
    int*** Y = new3d<int>(f,h,w);
//   int*** _Y = new3d<int>(f,h,w);
//    int*** U = new3d<int>(f,h/2,w/2);
//    int*** V = new3d<int>(f,h/2,w/2);
    int*** _Y = NULL;
    int*** U = NULL;
    int*** V = NULL;
    double* weights = (double*)malloc(sizeof(double)*PXL);
//    printf("%lf\n",pow(10,snr/10));
    for(int i = 0 ; i < len ; ++i)
        snr[i] = __SNR_S+i;

    // read YUV
    if(argc == 1)
        yuv_random_read("foreman_cif.yuv",h,w,f,Y,U,V);
    else
        yuv_read(argv[1],h,w,f,Y,NULL,NULL);

    // encode
//    weight_predict_minMSE(weights,pow(10,snr/10));


    for(int i = 0 ; i < len; ++i){
        for(int j = 0 ; j < PXL ; ++j)
            weights[j] = MINPOW_PSNR25[j];//DERIVE_SNR[i][j];

        video_encode(Y,f,h,w,snr[i],G,Ly,map_out,weights);
    //    video_encode(Y,f,h,w,snr,G,x,map_out);
    //    generate_Ly(x,lu,f,snr,Ly,weights);

        // decode
        if(argc == 1)
            DECODE("foreman",Y,Ly,map_out,h,w,f,lu,G,PSNR,_Y);
        else if(argc > 2)
            DECODE(argv[2],Y,Ly,map_out,h,w,f,lu,G,PSNR,NULL);
        else
            DECODE("NONAME",Y,Ly,map_out,h,w,f,lu,G,PSNR,NULL);

        // compute Ia
    //    double** Ia = new2d<double>(f,PXL);
    //    computeIa(Ly,map_out,Y,f,h,w,Ia);


        // print result PSNR
#if __PROGRESS__
        double avg_psnr = 0.0;
        for(int j = 0 ; j < f ; ++j){
            printf("frame#%3d PSNR = %lf\n",j+1,PSNR[j]);
            avg_psnr += PSNR[j];
        }
        printf("AVERAGE PSNR = %lf\n",avg_psnr/f);
#endif

#if __OUTPUT_SEQ__
        yuv_write("hall",_Y,U,V,f,h,w);
#endif
        write_psnr_info(PSNR,snr[i]);
    }

    // free memory
    delete3d<double>(Ly);
//    delete3d<int>(x);
    delete3d<int>(map_out);
    delete2d<int>(G);
    deleteY(Y);
//    deleteY(_Y);
//    deleteY(U);
//    deleteY(V);
    free(snr);
    free(PSNR);
    free(weights);
}
