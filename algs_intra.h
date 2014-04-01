#ifndef __ALGS_INTRA_H
#define __ALGS_INTRA_H

#include <stdio.h>
#include <stdlib.h>
#include "data_alloc.h"
#include "mrf_decoder_utils.h"
#include "image_process_utils.h"
#include "build_value.h"


/**
*   @Penlin: algorithm for intra decoding
*
*   @param: str     clip short name for record                          [char*]
*   @param: Y       original Y value of the clip , for compute PSNR     [n_frame*imgh*imgw]
*   @param: Ly_in   log likelihood value of bit planes                  [n_frame*PXL*(2*lu)]
*   @param: map_in  mapping for deinterleave                            [n_frame*PXL*lm]
*
*   @param: PSNR    computed PSNR between decoded frame and original frame [n_frame]
*   @param: imgr_out decoded clip                                       [n_frame*imgh*imgw]
*
*   @param: imgh, imgw, n_frame, lu, G
*
**/

void intra(const char* str, int*** Y, double*** Ly_in, int*** map_in, const int &imgh, const int &imgw, const int &n_frame, const int &lu, int** G, double* PSNR , int*** img_out = NULL){


    // param initial
    const int lm = imgh*imgw;
    double* Lu = (double*) malloc(sizeof(double)*lu);

    // frame buffer
    int** imgO;
    int** imgr = new2d<int>(imgh,imgw);
    int*** imgr_bp = new3d<int>(PXL,imgh,imgw);
    double** Ly;    //channel value
    int** map;      //interleaver map

    // pstate, pout
    const int Ns = pow(2,G_L-1);
    int ** pout = new2d<int>(Ns,4);
    int ** pstate = new2d<int>(Ns,2) ;
    double* Le1 = (double*)malloc(sizeof(double)*lu);
    double* Le2 = (double*)malloc(sizeof(double)*lu);

    for(int i = 0 ; i < lu ; ++i){
        Lu[i] = 0;
//        Le1[i] = Le2[i] = -log(2);
//        Le1[i] = Le2[i] = 0.5;
    }
    computeLe(Lu,Le1,Le2,lu);

    trellis(G,G_N,G_L,Ns,pout,pstate);


    double** beta = new2d<double>(n_frame,PXL);

    double** Lu_c = new2d<double>(PXL,lu);  //channel decoder output
    double** Le_c = new2d<double>(PXL,lm); //channel extrinsic information
    double** Lu_s = new2d<double>(PXL,lm); //source decoder output
    double** Le_s = new2d<double>(PXL,lm); // source extrinsic information


    // decoding
    for(int f = 0 ; f < n_frame ; ++f){
#if __PROGRESS__
        printf("Decoding frame#%d\n",f+1);
#endif

        // initialize channel/source value
        Ly = Ly_in[f];
        map = map_in[f];
        imgO = Y[f];

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

            for(int i = 0 ; i < PXL ; ++i)
                for(int j = 0 ; j < lm ; ++j)
                    //Le_s[i][j] = (Lu_s[i][j] - Le_c[i][j]);
                    Le_s[i][j] = (Lu_s[i][map[i][j]]-Le_c[i][map[i][j]]);

            // interleave
#if __STATUS__
        printf("interleave ...%lf\n",getCurrentTime());
#endif
            //interleave(Le_s,map,lm);

            // BCJR decoding
#if __STATUS__
        printf("BCJR decoding ...%lf\n",getCurrentTime());
#endif
            for(int t_lvl = 0,i = 0 ; t_lvl < PXL ; ++t_lvl){
                for(i = 0 ; i < lm ; ++i){
                    Lu[i] = Le_s[t_lvl][i];
//                    Le1[i] = -log(1+exp(Lu[i]));
//                    Le2[i] = (Le1[i] + Lu[i]);
//                    Le1[i] = 1/(1+exp(Lu[i]));      // P0
//                    Le2[i] = 1-Le1[i];              // P1
                }
                computeLe(Lu,Le1,Le2,lm);
                BCJR_decoding(Ns, lu, 1, Ly[t_lvl], Le1, Le2, pstate, pout, Lu_c[t_lvl]);
            }

            // deinterleave
#if __STATUS__
        printf("deinterleave ...%lf\n",getCurrentTime());
#endif
            for(int i = 0 ; i < PXL ; ++i)
                for(int j = 0 ; j < lm ; ++j)
                    Le_c[i][map[i][j]] = (Lu_c[i][j] - Le_s[i][j]);
            //deinterleave(Le_c,map,lm);

            // MRF parameter estimation
#if __STATUS__
        printf("MRF parameter estimation ...%lf\n",getCurrentTime());
#endif
            for(int i = 0 ; i <imgh ; ++i)
                for(int j = 0 ; j < imgw ; ++j)
                    for(int t_lvl=0 ; t_lvl < PXL ; ++t_lvl)
                        imgr_bp[t_lvl][i][j] = ((Le_c[t_lvl][j+i*imgw]>=0)?1:0);

            intra_beta_estimation(imgr_bp,beta[f],imgh,imgw);


            // MRF decoding
#if __STATUS__
        printf("MRF decoding ... %lf\n",getCurrentTime());
#endif
            mrf_siso_intra(Le_c,beta[f],imgh,imgw,Lu_s,1);
#if __BETA__
        printf("frame#%d iter#%d,\n",f+1,iter+1);
        for(int i = 0 ; i < PXL ; ++i)
            printf("bp %d beta = %lf \n",i,beta[f][i]);
#endif

            // recover to image
            for(int i = 0 ; i <imgh ; ++i){
                for(int j = 0 ; j < imgw ; ++j){
                    for(int t_lvl=0 ; t_lvl < PXL ; ++t_lvl)
                        imgr_bp[t_lvl][i][j] = ((Lu_c[t_lvl][j+i*imgw]>=0)?1:0);
                }
            }

            // construct image
            bin2dec_img(imgr_bp,imgh,imgw,imgr);

            // compute PSNR
            double channel_psnr = computePSNR(imgr,imgO,imgh,imgw);

            // recover to image
            for(int i = 0 ; i <imgh ; ++i){
                for(int j = 0 ; j < imgw ; ++j){
                    for(int t_lvl=0 ; t_lvl < PXL ; ++t_lvl)
                        imgr_bp[t_lvl][i][j] = ((Lu_s[t_lvl][j+i*imgw]>=0)?1:0);
                }
            }

            // construct image
            bin2dec_img(imgr_bp,imgh,imgw,imgr);

            // compute PSNR
            PSNR[f] = computePSNR(imgr,imgO,imgh,imgw);

            //printf("frame#%d iter#%d,avg beta=\n",f+1,iter+1);
            for(int i = 0 ; i < PXL ; ++i)
               printf("%lf,",beta[f][i]);
            printf("%lf\n",PSNR[f]-channel_psnr);
#if __PSNR__
            printf("%s frame#%d PSNR_iter%d = %lf\n",str,f+1,iter+1,PSNR[f]);
#endif
        }

        // imgr output
        if(img_out!=NULL)
            for(int i = 0, j = 0; i <imgh ; ++i)
                for(j=0 ; j < imgw ; ++j)
                    img_out[f][i][j] = imgr[i][j];

    }

    delete3d<int>(imgr_bp);
    delete2d<int>(imgr);
//    delete2d<int>(G);
    delete2d<int>(pstate);
    delete2d<int>(pout);

    delete2d<double>(beta);
    free(Lu);
    free(Le1);
    free(Le2);

    delete2d<double>(Lu_c);
    delete2d<double>(Lu_s);
    delete2d<double>(Le_c);
    delete2d<double>(Le_s);
}


#endif // __ALGS_INTRA_H
