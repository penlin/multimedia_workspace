#ifndef __UEP_PREDICT_INTER_INTRA_H
#define __UEP_PREDICT_INTER_INTRA_H

#include "uep_predict_utils.h"

// for consective inter
void weight_predict_Inter_Intra_minMSE(Pixel** img, Pixel** img_prev, Pixel** img_nxt, const int &imgh, const int &imgw, double* weights, const double &gamma){

    double value1 = 0.0, value2 = 0.0, beta1, beta2;
    double eEn[PXL];
    int i,j,k,bit=1,mbSize=8,n_block=0;
    int8** img_bp = new2d<int8>(imgh,imgw,0);
    int8** img_bp_prev = new2d<int8>(imgh,imgw,0);
    int8** img_bp_nxt = new2d<int8>(imgh,imgw,0);
    int** mv_prev = new2d<int>(2,imgh*imgw/mbSize/mbSize,0);
    int** mv_nxt = new2d<int>(2,imgh*imgw/mbSize/mbSize,0);

    motionEstES(img,img_prev,imgh,imgw,mbSize,5,mv_prev);
    motionEstES(img,img_nxt,imgh,imgw,mbSize,5,mv_nxt);

    for(i = 0, bit = (1<<(PXL-1)); i < PXL ; ++i, bit>>=1){
        for(j=0;j<imgh;++j){
            for(k=0;k<imgw;++k){
                n_block = ((j/mbSize)*imgw/mbSize+k/mbSize);
                img_bp[j][k] = ((img[j][k] & bit)>0);
                img_bp_prev[j][k] = ((img_prev[j+mv_prev[0][n_block]][k+mv_prev[1][n_block]] & bit)>0);
                img_bp_nxt[j][k] = ((img_nxt[j+mv_nxt[0][n_block]][k+mv_nxt[1][n_block]] & bit)>0);
            }
        }
        inter2_beta_estimation(img_bp,img_bp_prev,img_bp_nxt,beta1,beta2,imgh,imgw);
        value1 = value2 = 0;
        for(j=0;j<imgh;++j)
            for(k=0;k<imgw;++k){
                value1+=(2*img_bp[j][k]-1)*(2*img_bp_prev[j][k]-1);
                value2+=(2*img_bp[j][k]-1)*(2*img_bp_nxt[j][k]-1);
            }

        eEn[i] = exp((value1*beta1+value2*beta2)/(imgh*imgw));
    }

    delete2d<int8>(img_bp);
    delete2d<int8>(img_bp_prev);
    delete2d<int8>(img_bp_nxt);
    delete2d<int>(mv_nxt);
    delete2d<int>(mv_prev);

    weight_predict_basic(weights,eEn,gamma);

}

// for inter pair
void weight_predict_Inter_Intra_minMSE(Pixel** img, Pixel** img_ref,const int &imgh, const int &imgw, double* weights, const double &gamma){

    double value = 0.0,beta;
    double eEn[PXL];
    int i,j,k,bit=1,mbSize=8,n_block=0;
    int8** img_bp = new2d<int8>(imgh,imgw,0);
    int8** img_bp_ref = new2d<int8>(imgh,imgw,0);
    int** mv = new2d<int>(2,imgh*imgw/mbSize/mbSize,0);

    motionEstES(img,img_ref,imgh,imgw,mbSize,5,mv);

    for(i = 0, bit = (1<<(PXL-1)); i < PXL ; ++i, bit>>=1){
        for(j=0;j<imgh;++j){
            for(k=0;k<imgw;++k){
                n_block = ((j/mbSize)*imgw/mbSize+k/mbSize);
                img_bp[j][k] = ((img[j][k] & bit)>0);
                img_bp_ref[j][k] = ((img_ref[j+mv[0][n_block]][k+mv[1][n_block]] & bit)>0);
            }
        }
        inter_beta_estimation(img_bp,beta,img_bp_ref,imgh,imgw);
        value= 0;
        for(j=0;j<imgh;++j)
            for(k=0;k<imgw;++k)
                value+=(2*img_bp[j][k]-1)*(2*img_bp_ref[j][k]-1);

        eEn[i] = exp(value*beta/(imgh*imgw));
    }

    delete2d<int8>(img_bp);
    delete2d<int8>(img_bp_ref);
    delete2d<int>(mv);

    weight_predict_basic(weights,eEn,gamma);

}

double inter_intra_basic_psnr_est(int8*** img_bp, int8*** img_bp_ref, int** mv1, const int &imgh, const int &imgw, double* weights, const double &gamma){

    double fr[PXL];
    double frr[PXL]={0,0,0,0,0,0,0,0};
    double mse = 0, mse_ori = 0,psnr = 0, psnr_ori = 0;
    double value = 0.0, En = 0.0, beta = 0.0 ;
    double suppress = 1.0;
    int i,j,k;
    int8** img_ref = new2d<int8>(imgh,imgw);

    for(i = 0 ; i < PXL ; ++i){

        value = weights[i]*gamma;
        cutGAMMA(value);
        interp2(value,fr[i]);
        suppress = 1-2*fr[i];
        motionComp(img_bp_ref[i],mv1,imgh,imgw,8,img_ref);

        En = 0.0;
        inter_beta_estimation(img_bp[i],beta,img_ref,imgh,imgw);
        for(j= 0 ; j < imgh ; ++j)
            for(k=0 ; k < imgw ; ++k)
                En+= (2*img_bp[i][j][k]-1)*(2*img_ref[j][k]-1);

        En = En*beta*suppress*suppress/(imgh*imgw);
        frr[i] = fr[i]/(fr[i]+(1-fr[i])*exp(En));

        mse+=(frr[i]*ORDER2[i]);
        mse_ori+=(fr[i]*ORDER2[i]);
//        printf("%lf,",frr[i]);
    }
    psnr = 10*log10(65025/mse);
    psnr_ori = 10*log10(65025/mse_ori);
//    for(i = 0 ; i < PXL ; ++i)
//        printf("%lf,",fr[i]);
//    printf("%lf\n",psnr_ori);

    for(i = 0 ; i < PXL ; ++i)
        printf("%lf,",frr[i]);
    printf("%lf\n",psnr);

    //printf("%lf,%lf\n",psnr_ori,psnr);

    delete2d<int8>(img_ref);
    return psnr;
}


double inter_intra_psnr_est(int8*** img_bp, int8*** img_bp_ref, int** mv1, const int &imgh, const int &imgw, double* weights, const double &gamma,double* frr_prev, int8*** img_bp_ref2 = NULL, int** mv2 = NULL){

    double fr[PXL];
    double frr[PXL]={0,0,0,0,0,0,0,0};
    double mse = 0, mse_ori = 0,psnr = 0, psnr_ori = 0;
    double value = 0.0, En = 0.0, En1 = 0.0, beta = 0.0, beta_prev = 0.0 ;
    double suppress = 1.0;
    int i,j,k;
    int8** img_ref = new2d<int8>(imgh,imgw);
    int8** img_ref2;

    if(img_bp_ref2!=NULL)
        img_ref2 = new2d<int8>(imgh,imgw);

    for(i = 0 ; i < PXL ; ++i){

        value = weights[i]*gamma;
        cutGAMMA(value);
        interp2(value,fr[i]);
        suppress = 1-2*fr[i];
        if(img_bp_ref2==NULL && frr_prev[i]>0){
            suppress = 1-2*frr_prev[i];
        }
        motionComp(img_bp_ref[i],mv1,imgh,imgw,8,img_ref);
        En = En1 = 0.0;
        if(img_bp_ref2==NULL){
            inter_beta_estimation(img_bp[i],beta,img_ref,imgh,imgw);
            for(j= 0 ; j < imgh ; ++j)
                for(k=0 ; k < imgw ; ++k)
                    En+= (2*img_bp[i][j][k]-1)*(2*img_ref[j][k]-1);

            En = En*beta*suppress/(imgh*imgw);
            frr[i] = fr[i]/(fr[i]+(1-fr[i])*exp(En));
        } else {
            motionComp(img_bp_ref2[i],mv2,imgh,imgw,8,img_ref2);
            inter2_beta_estimation(img_bp[i], img_ref2, img_ref, beta, beta_prev, imgh, imgw);
            for(j= 0 ; j < imgh ; ++j)
                for(k=0 ; k < imgw ; ++k){
                    En+= (2*img_bp[i][j][k]-1)*(2*img_ref[j][k]-1);
                    En1+= (2*img_bp[i][j][k]-1)*(2*img_ref2[j][k]-1);
                }

            En = En*beta_prev*(1-2*frr_prev[i])/(imgh*imgw);
            En1 = En1*beta*suppress/(imgh*imgw);
            frr[i] = fr[i]/(fr[i]+(1-fr[i])*exp(En+En1));
        }

        mse+=(frr[i]*ORDER2[i]);
        mse_ori+=(fr[i]*ORDER2[i]);
//        printf("%lf,",frr[i]);
        frr_prev[i] = frr[i];
    }
    psnr = 10*log10(65025/mse);
    psnr_ori = 10*log10(65025/mse_ori);
//    for(i = 0 ; i < PXL ; ++i)
//        printf("%lf,",fr[i]);
//    printf("%lf\n",psnr_ori);

    for(i = 0 ; i < PXL ; ++i)
        printf("%lf,",frr[i]);
    printf("%lf\n",psnr);

    //printf("%lf,%lf\n",psnr_ori,psnr);

    delete2d<int8>(img_ref);
    if(img_bp_ref2!=NULL){
        delete2d<int8>(img_ref2);
    }
    return psnr;
}


#endif // __UEP_PREDICT_INTER_INTRA_H
