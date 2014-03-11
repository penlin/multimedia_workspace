#include "build_value.h"
#include "io_utils.h"
#include "video_encode.h"
#include "data_alloc.h"
#include "algs_direct.h"
#include "algs_intra.h"
#include "algs_inter.h"
#include "algs_inter_intra.h"
#include "utils.h"

#define DECODE __ALGO__

//const double LAST_GOOD[PXL] = {2.50,2.00,1.60,1.00,0.40,0.50,0.0,0.0};
//const double LAST_GOOD[PXL] = {2.50,2.10,1.70,1.10,0.25,0.35,0.0,0.0};
const double LAST_GOOD[PXL] = {2.52,2.10,1.69,1.10,0.26,0.33,0.01,0.00};

const double LAST_GOOD_SNR_2[PXL] = {2.22,1.75,1.59,1.15,0.86,0.38,0.06,0.00}; // 35 dB
const double LAST_GOOD_SNR_4[PXL] = {2.02,1.50,1.39,1.05,0.86,0.68,0.41,0.10}; // 46 dB
const double DERIVE_SNR_4[PXL] = {1.7441,1.5094, 1.2528,1.1347, 0.9412, 0.7502,0.5525, 0.1152}; // 47.8 dB
const double LAST_GOOD_SNR_6[PXL] = {1.97,1.50,1.09,0.90,0.86,0.68,0.56,0.45}; // 61 dB

int main(int argc,char* argv[]){

    startRandom();

    const int h = __HEIGHT, w = __WIDTH, f = __FRAME ;
    const int lm = h*w;
    const int lu = lm+(G_L-1);
    const int snr = __SNR;
    int ** G = getGenerator();
    double* PSNR = (double*) malloc(sizeof(double)*f);
    double *** Ly = new3d<double>(f,PXL,2*lu);
    int *** x = new3d<int>(f,PXL,2*lu);
    int *** map_out = new3d<int>(f,PXL,lm);
    int*** Y = new3d<int>(f,h,w);
    int*** _Y = new3d<int>(f,h,w);
    int*** U = new3d<int>(f,h/2,w/2);
    int*** V = new3d<int>(f,h/2,w/2);

    double* weights_prev = (double*)malloc(sizeof(double)*PXL);
    double* weights = (double*)malloc(sizeof(double)*PXL);
    for(int i = 0 ; i < PXL ; ++i)
        weights_prev[i] = DERIVE_SNR_4[i];

    // read YUV
    yuv_random_read("hall_cif.yuv",h,w,f,Y,NULL,NULL);

    // encode
    video_encode(Y,f,h,w,snr,G,x,map_out);

    double delta = 0.0;
    double performance[PXL] = {0.0};
    double lasthighest = -100.0;
    double highest = -1.0, lowest = 10000.0;
    const double tolerance = 0.07;
    int index_h = 0, index_l = 0;

    while(highest > lasthighest - tolerance){

        lasthighest = highest;
        highest = -1.0;
        lowest = 10000.0;
        index_h = index_l = -1;

        for(int i = 0 ; i < PXL ; ++i){

            // update weights
            for(int j = 0 ; j < PXL ; ++j){
                if(i==j)
                    weights[j] = weights_prev[j] + delta;
                else
                    weights[j] = weights_prev[j];
            }

            generate_Ly(x,lu,f,snr,Ly,weights);

            // decode
            DECODE("hall",Y,Ly,map_out,h,w,f,lu,G,PSNR,_Y);
            // print result PSNR
            performance[i] = 0.0;
            for(int frame = 0 ; frame < f ; ++frame){
                //printf("frame#%3d PSNR = %lf\n",i+1,PSNR[i]);
                performance[i] += PSNR[frame];
            }
            performance[i]/=f;
            printf("(%d)AVERAGE PSNR = %lf\n",i+1,performance[i]);

#if __OUTPUT_SEQ__
    yuv_write("hall",_Y,U,V,f,h,w);
#endif
            write_power_info(weights,performance[i]);

            if(performance[i]>highest){
                index_h = i;
                highest = performance[i];
            }

            if(performance[i]<lowest && weights_prev[i]>=delta-delta/100 ){
                index_l = i;
                lowest = performance[i];
            }

        }

        weights_prev[index_h] += delta;
        if(index_l>=0)
            weights_prev[index_l] -= delta;
        else
            break;


    }

    // free memory
    delete3d<double>(Ly);
    delete3d<int>(x);
    delete3d<int>(map_out);
    delete2d<int>(G);
    deleteY(Y);
    deleteY(_Y);
    deleteY(U);
    deleteY(V);

    free(PSNR);
    free(weights);
}
