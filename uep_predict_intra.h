#ifndef __UEP_PREDICT_INTRA_H
#define __UEP_PREDICT_INTRA_H

#include "uep_predict_utils.h"


void weight_predict_Intra_minMSE(int8*** img_bp,const int &imgh, const int &imgw, double* weights, const double &gamma){

    double value = 0.0,beta;
    double eEn[PXL];
    int i,j,k,N=0;

    for(i = 0; i < PXL ; ++i){

        intra_beta_estimation(img_bp[i],beta,imgh,imgw);
        value= 0;
        for(j=1;j<imgh-1;++j){
            for(k=1;k<imgw-1;++k){
                N = img_bp[i][j-1][k]+img_bp[i][j+1][k]+img_bp[i][j][k-1]+img_bp[i][j][k+1];
                value+=(2*img_bp[i][j][k]-1)*(2*N-4);
            }
        }
        eEn[i] = exp(value*beta/((imgh-1)*(imgw-1)));
    }

    weight_predict_basic(weights,eEn,gamma);
}



double intra_psnr_est(int8*** img_bp,const int &imgh, const int &imgw, double* weights, const double &gamma){

    double fr[PXL];
    double frr[PXL]={0,0,0,0,0,0,0,0};
    double mse = 0, mse_ori = 0,psnr = 0, psnr_ori = 0;
    int N = 0;
    double value = 0.0;
    int i,j,k;
    double beta = 0.0;

    for(i = 0 ; i < PXL ; ++i){
        value = weights[i]*gamma;
        cutGAMMA(value);
        interp2(value,fr[i]);

        intra_beta_estimation(img_bp[i],beta,imgh,imgw);
        value=0;
        for(j=1;j<imgh-1;++j){
            for(k=1;k<imgw-1;++k){
                N = img_bp[i][j-1][k]+img_bp[i][j+1][k]+img_bp[i][j][k-1]+img_bp[i][j][k+1];
                value+=(2*img_bp[i][j][k]-1)*(2*N-4);
            }
        }
        value=value*beta*(1-2*fr[i])/(imgh*imgw);

        frr[i] = fr[i]/(fr[i]+(1-fr[i])*exp(value)); // N effecitve than n
        mse+=(frr[i]*ORDER2[i]);
        mse_ori+=(fr[i]*ORDER2[i]);
        //printf("%lf,",frr[i]);
    }
    psnr = 10*log10(65025/mse);
    psnr_ori = 10*log10(65025/mse_ori);
//    for(i = 0 ; i < PXL ; ++i)
//        printf("%lf,",fr[i]);
//    printf("%lf\n",psnr_ori);
    for(i = 0 ; i < PXL ; ++i)
        printf("%lf,",frr[i]);
    printf("%lf\n",psnr);
    return psnr;
}



#endif // __UEP_PREDICT_INTRA_H
