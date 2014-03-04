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
const double LAST_GOOD[PXL] = {2.50,2.10,1.70,1.10,0.25,0.35,0.0,0.0};

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
        weights_prev[i] = LAST_GOOD[i];

    // read YUV
    yuv_random_read("hall_cif.yuv",h,w,f,Y,NULL,NULL);

    // encode
    video_encode(Y,f,h,w,snr,G,x,map_out);

    double delta = 0.05;
    double performance[PXL] = {0.0};
    double lasthighest = -100.0;
    double highest = -1.0, lowest = 10000.0;
    int index_h = 0, index_l = 0;

    while(highest > lasthighest){

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

            if(performance[i]<lowest && weights_prev[i]>delta/2 ){
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
