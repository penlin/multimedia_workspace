#ifndef __ALGS_DIRECT_H
#define __ALGS_DIRECT_H
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

void direct_system(const char* str, FILE* fptr, const size_t &imgh, const size_t &imgw, const size_t &n_frame,int** G,const double &snr, double* PSNR ,int weight_type = 0, int*** img_out = NULL){


    Frame* frame = new Frame(imgh,imgw,0,0);
    frame->encode_info(snr,G);

    // param initial
    const size_t lm = imgh*imgw;
    const size_t lu = lm + 2;

    double* weights = MALLOC(double,PXL);

    if(weight_type)
        weight_predict_minMSE(weights,pow(10,snr/10));
    else
        for(int j = 0 ; j < PXL ; ++j)
            weights[j] = 1;

    double* Le1 = MALLOC(double,lu);
    double* Le2 = MALLOC(double,lu);
    for(int i = 0 ; i < lu ; ++i)
        Le1[i] = Le2[i] = 0.5;

    // frame buffer
    PIXEL** imgr = new2d<PIXEL>(imgh,imgw);
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
        for(int t_lvl = 0 ; t_lvl < PXL ; ++t_lvl)
            BCJR_decoding( lu, 1, Ly[t_lvl], Le1, Le2, Lu_c[t_lvl]);

        // deinterleave
#if __STATUS__
        printf("deinterleave ...%lf\n",getCurrentTime());
#endif
//        deinterleave(Lu_c,map,lm);

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
    delete2d<PIXEL>(imgr);
    DELETE(Le1);
    DELETE(Le2);
    DELETE(weights);

    delete frame;
}


#endif // __ALGS_DIRECT_H
