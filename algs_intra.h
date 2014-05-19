#ifndef __ALGS_INTRA_H
#define __ALGS_INTRA_H
#include <omp.h>
#include <stdio.h>
#include <stdlib.h>
#include "frame.h"
#include "mrf_decoder_utils.h"
#include "uep_predict_intra.h"

/**
*   @Penlin: algorithm for intra decoding
*
*   @param: str     clip short name for record                          [char*]

*
*   @param: PSNR    computed PSNR between decoded frame and original frame [n_frame]
*   @param: imgr_out decoded clip                                       [n_frame*imgh*imgw]
*
*   @param: imgh, imgw, n_frame, lu, G
*
**/

void intra_system(const char* str, FILE* fptr, const int &imgh, const int &imgw, const int &n_frame,int** G, const double &snr, double* PSNR , int repredict = 0, int*** img_out = NULL){


    Frame* frame = new Frame(imgh,imgw,0,0);
    frame->encode_info(snr,G);
    const double EbN0 = pow(10,snr/10);
    // param initial
    const int lm = imgh*imgw;
    const int lu = lm + 2;

    double* ber = MALLOC(double, PXL);

    double* weights = MALLOC(double,PXL);
    for(int i = 0 ; i < PXL ; ++i)
        weights[i] = 1;

    // frame buffer
    Pixel** imgr = new2d<Pixel>(imgh,imgw);
    int8*** imgr_bp = new3d<int8>(PXL,imgh,imgw);
    double** Ly = new2d<double>(PXL,2*lu,0);    //channel value
    int** map = new2d<int>(PXL,lm);      //interleaver map

    // pstate, pout
    double** Le1 = new2d<double>(PXL,lu,0.5);//MALLOC(double,lu);
    double** Le2 = new2d<double>(PXL,lu,0.5);// MALLOC(double,lu);

    double* beta = MALLOC(double,PXL);

    double** Lu_c = new2d<double>(PXL,lu);  //channel decoder output
    double** Le_c = new2d<double>(PXL,lm); //channel extrinsic information
    double** Lu_s = new2d<double>(PXL,lm); //source decoder output
    double** Le_s = new2d<double>(PXL,lm); // source extrinsic information


    // decoding
    for(int f = 0 ; f < n_frame ; ++f){
#if __PROGRESS__
        printf("Encoding frame#%d\n",f+1);
#endif

        frame->read(fptr);
        if(repredict)
            weight_predict_Intra_minMSE(frame->img_bp,imgh,imgw,weights,EbN0);
        frame->encode(Ly,map,weights);

#if __PROGRESS__
        printf("Decoding frame#%d\n",f+1);
#endif
        // initialize channel/source value

        for(int i = 0 ; i < PXL ; ++i){

            for(int j = 0 ; j < lm ; ++j)
                Le_c[i][j] = Lu_s[i][j] = Le_s[i][j] = 0;

            for(int j = 0 ; j < lu ; ++j)
                Lu_c[i][j] = 0;
        }

        for(int iter = 0 ; iter < Niter ; ++iter){
#if __PROGRESS__
            printf("iter #%d\n",iter+1);
#endif

            // BCJR decoding
#if __STATUS__
            printf("BCJR decoding ...%lf\n",getCurrentTime());
#endif
            #pragma omp parallel for
            for(int t_lvl = 0 ; t_lvl < PXL ; ++t_lvl){
                for(int i = 0 ; i < lm ; ++i)
                    Le_s[t_lvl][i] = (Lu_s[t_lvl][map[t_lvl][i]]-Le_c[t_lvl][map[t_lvl][i]]);
                computeLe(Le_s[t_lvl],Le1[t_lvl],Le2[t_lvl],lm);
                BCJR_decoding(lu, 1, Ly[t_lvl], Le1[t_lvl], Le2[t_lvl], Lu_c[t_lvl]);

                // deinterleave
                for(int i = 0 ; i < lm ; ++i){
                    Le_c[t_lvl][map[t_lvl][i]] = (Lu_c[t_lvl][i] - Le_s[t_lvl][i]);
                    imgr_bp[t_lvl][0][map[t_lvl][i]] =((Le_c[t_lvl][map[t_lvl][i]]>=0)?1:0);
                }

                // MRF parameter estimation
                intra_beta_estimation(imgr_bp[t_lvl],beta[t_lvl],imgh,imgw);
            }

            // MRF decoding
#if __STATUS__
        printf("MRF decoding ... %lf\n",getCurrentTime());
#endif
            mrf_siso_intra(Le_c,beta,imgh,imgw,Lu_s,1);

#if __BETA__
        printf("frame#%d iter#%d,\n",f+1,iter+1);
        for(int i = 0 ; i < PXL ; ++i)
            printf("bp %d beta = %lf \n",i,beta[i]);
#endif

            // recover to image
            Lu2dec_img(Lu_c,lm,imgr,map);

            // compute PSNR
            double channel_psnr = frame->psnr(imgr);

            // recover to image
            Lu2dec_img(Lu_s,lm,imgr);
            // compute PSNR
            PSNR[f] = frame->psnr(imgr);

#if __PSNR__
            printf("%s frame#%d PSNR_iter%d = %lf (channel:%lf)\n",str,f+1,iter+1,PSNR[f],channel_psnr);
#endif
        }

#if !__PSNR__
        computeBER(imgr,frame->Y,lm,ber);
        printf("%lf\n",PSNR[f]);
#endif

        // imgr output
        if(img_out!=NULL)
            for(int i = 0, j = 0; i <imgh ; ++i)
                for(j=0 ; j < imgw ; ++j)
                    img_out[f][i][j] = imgr[i][j];

    }


    delete2d<double>(Ly);
    delete2d<int>(map);

    delete3d<int8>(imgr_bp);
    delete2d<Pixel>(imgr);

    DELETE(beta);
    DELETE(weights);
    DELETE(ber);

    delete2d<double>(Le1);
    delete2d<double>(Le2);
    delete2d<double>(Lu_c);
    delete2d<double>(Lu_s);
    delete2d<double>(Le_c);
    delete2d<double>(Le_s);

    delete frame;
}


#endif // __ALGS_INTRA_H
