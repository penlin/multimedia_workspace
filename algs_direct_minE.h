#ifndef __ALGS_DIRECT_H
#define __ALGS_DIRECT_H
#include <omp.h>
#include <stdio.h>
#include <stdlib.h>
#include "frame.h"
#include "uep_predict_utils.h"

/**
*   @Penlin: algorithm for direct decode
*
*   @param: str     clip short name for record                          [char*]
*
*   @param: PSNR    computed PSNR between decoded frame and original frame [n_frame]
*   @param: imgr_out decoded clip                                       [n_frame*imgh*imgw]
*
*   @param: imgh, imgw, n_frame, lu, G
*
**/

void direct_minE_system(const char* str, FILE* fptr, const int &imgh, const int &imgw, const int &n_frame,int** G, double mse, double* PSNR ,int weight_type = 0, int*** img_out = NULL){


    double* weights = MALLOC(double,PXL);
    double snr = 0;

    if(weight_type){
        weight_predict_minPOWER(weights,mse);
        snr = 0.0;
        for(int j = 0 ; j < PXL ; ++j)
            snr += weights[j];
        snr /= PXL;
        for(int j = 0 ; j < PXL ; ++j)
            weights[j]/=snr;
        snr = 10*log10(snr);
        printf("Gamma_0 = %lf dB :[",snr);
        for(int j = 0 ; j < PXL ; ++j)
            printf("%lf, ",weights[j]);
        printf("]\n");
    }else{
        for(int j = 0 ; j < PXL ; ++j)
            weights[j] = 1;
        interpf(mse/21845,snr);
        snr = 10*log10(snr);
    }

    Frame* frame = new Frame(imgh,imgw,0,0);
    frame->encode_info(snr,G);

    // param initial
    const int lm = imgh*imgw;
    const int lu = lm + 2;



    double* Le1 = MALLOC(double,lu);
    for(int i=0; i<lu ;++i)
        Le1[i] = 0.5;

    // frame buffer
    Pixel** imgr = new2d<Pixel>(imgh,imgw);
    double** Ly = new2d<double>(PXL,2*lu,0);    //channel value
    int** map = new2d<int>(PXL,lm);      //interleaver map

    double** Lu_c = new2d<double>(PXL,lu,0);  //channel decoder output


    // decoding
    for(int f = 0 ; f < n_frame ; ++f){
        printf("Encoding frame#%d\n",f+1);

        frame->next(fptr,Ly,map,weights);

        printf("Decoding frame#%d\n",f+1);
        // BCJR decoding
#if __STATUS__
        printf("BCJR decoding ...%lf\n",getCurrentTime());
#endif
        #pragma omp parallel for
        for(int t_lvl = 0 ; t_lvl < PXL ; ++t_lvl)
            BCJR_decoding(lu, 1, Ly[t_lvl], Le1, NULL , Lu_c[t_lvl]);

        // deinterleave
#if __STATUS__
        printf("deinterleave ...%lf\n",getCurrentTime());
#endif

        // recover to image
        Lu2dec_img(Lu_c,lm,imgr,map);

        // compute PSNR
        PSNR[f] = frame->psnr(imgr);
#if __PSNR__
        printf("%s frame#%d PSNR = %lf\n",str,f+1,PSNR[f]);
#endif

        // imgr output
        if(img_out!=NULL)
            for(int i = 0, j = 0; i <imgh ; ++i)
                for(j=0 ; j < imgw ; ++j)
                    img_out[f][i][j] = imgr[i][j];

    }

    delete2d<double>(Ly);
    delete2d<int>(map);
    delete2d<double>(Lu_c);
    delete2d<Pixel>(imgr);
    DELETE(Le1);
    DELETE(weights);

    delete frame;
}


#endif // __ALGS_DIRECT_H