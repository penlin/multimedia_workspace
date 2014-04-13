#ifndef __ALGS_INTRA_H
#define __ALGS_INTRA_H

#include <stdio.h>
#include <stdlib.h>
#include "frame.h"
#include "mrf_decoder_utils.h"


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

void intra_system(const char* str, FILE* fptr, const int &imgh, const int &imgw, const int &n_frame,int** G, int** pout, int** pstate,const double &snr, double* weights,double* PSNR , int*** img_out = NULL){


    Frame* frame = new Frame(imgh,imgw,0,0);
    frame->encode_info(snr,G);

    // param initial
    const int Ns = pow(2,G_L-1);
    const int lm = imgh*imgw;
    const int lu = lm + 2;

    // frame buffer
    int** imgr = new2d<int>(imgh,imgw);
    int*** imgr_bp = new3d<int>(PXL,imgh,imgw);
    double** Ly = new2d<double>(PXL,2*lu,0);    //channel value
    int** map = new2d<int>(PXL,lm);      //interleaver map

    // pstate, pout
    double* Lu = MALLOC(double,lu);
    double* Le1 = MALLOC(double,lu);
    double* Le2 = MALLOC(double,lu);
    for(int i = 0 ; i < lu ; ++i)
        Lu[i] = 0;

    computeLe(Lu,Le1,Le2,lu);

    double* beta = MALLOC(double,PXL);

    double** Lu_c = new2d<double>(PXL,lu);  //channel decoder output
    double** Le_c = new2d<double>(PXL,lm); //channel extrinsic information
    double** Lu_s = new2d<double>(PXL,lm); //source decoder output
    double** Le_s = new2d<double>(PXL,lm); // source extrinsic information


    // decoding
    for(int f = 0 ; f < n_frame ; ++f){
#if __PROGRESS__
        printf("Decoding frame#%d\n",f+1);
#endif

        frame->next(fptr,Ly,map,weights);
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

            for(int i = 0 ; i < PXL ; ++i)
                for(int j = 0 ; j < lm ; ++j)
                    Le_s[i][j] = (Lu_s[i][map[i][j]]-Le_c[i][map[i][j]]);


            // BCJR decoding
#if __STATUS__
            printf("BCJR decoding ...%lf\n",getCurrentTime());
#endif
            for(int t_lvl = 0,i = 0 ; t_lvl < PXL ; ++t_lvl){
                for(i = 0 ; i < lm ; ++i)
                    Lu[i] = Le_s[t_lvl][i];

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

            intra_beta_estimation(imgr_bp,beta,imgh,imgw);


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
            Lu2dec_img(Lu_c,imgh,imgw,imgr,map);
            // compute PSNR
            double channel_psnr = frame->psnr(imgr);

            // recover to image
            Lu2dec_img(Lu_s,imgh,imgw,imgr);
            // compute PSNR
            PSNR[f] = frame->psnr(imgr);

            //printf("frame#%d iter#%d,avg beta=\n",f+1,iter+1);
            for(int i = 0 ; i < PXL ; ++i)
               printf("%lf,",beta[i]);
            printf("%lf,%lf\n",channel_psnr,PSNR[f]-channel_psnr);
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


    delete2d<double>(Ly);
    delete2d<int>(map);

    delete3d<int>(imgr_bp);
    delete2d<int>(imgr);

    free(beta);
    free(Lu);
    free(Le1);
    free(Le2);

    delete2d<double>(Lu_c);
    delete2d<double>(Lu_s);
    delete2d<double>(Le_c);
    delete2d<double>(Le_s);

    delete frame;
}


#endif // __ALGS_INTRA_H
