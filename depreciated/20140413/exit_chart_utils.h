#ifndef __EXIT_CHART_UTILS_H
#define __EXIT_CHART_UTILS_H

#include <stdio.h>
#include <math.h>
#include "channel_code_utils.h"
#include "utils.h"

#define PXL 8

//I = 1- sum(log2(1 + exp(-L .* (2*x -1))),1)./size(x,1);
void mutual_information(double** Le, int*** imgr_bp, const int &lm, double* Ie_c){

    double tmp_sum = 0;

    for(int t_lvl = 0 , i = 0; t_lvl < PXL ; ++t_lvl){
        tmp_sum = 0;
        for(i=0; i < lm ; ++i){
            tmp_sum+=log2(1 + exp(-Le[t_lvl][i]*(2*imgr_bp[0][0][i+t_lvl*lm]-1)));

        }
        Ie_c[t_lvl] = 1 - tmp_sum/lm;
//        printf("temp_sum = %lf, Ie_c = %lf\n",tmp_sum,Ie_c[t_lvl]);
    }
}


void computeIa(double*** Ly, int*** map, int*** imgO, const int &n_frame, const int &imgh, const int &imgw , double** Ia){

    int*** imgr_bp = new3d<int>(PXL,imgh,imgw);
    const int lm = imgh*imgw;

    for(int i = 0 ; i < n_frame ; ++i){
        img2bp_frame(imgO[i],imgh,imgw,imgr_bp);
        deinterleave(Ly[i],map[i],lm);
        mutual_information(Ly[i],imgr_bp,lm,Ia[i]);
    }

    delete3d<int>(imgr_bp);
}


#endif // __EXIT_CHART_UTILS_H
