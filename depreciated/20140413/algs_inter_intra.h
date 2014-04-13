#ifndef __ALGS_INTER_INTRA_H
#define __ALGS_INTER_INTRA_H

#include <stdio.h>
#include <stdlib.h>
#include "data_alloc.h"
#include "build_value.h"
#include "exit_chart_utils.h"
#include "io_utils.h"

/**
*   @Penlin: algorithm for inter+intra  decoding
*
*   @param: str     clip short name for record                          [char*]
*   @param: Y       original Y value of the clip , for compute PSNR     [imgh*imgw*n_frame]
*   @param: Ly_in   log likelihood value of bit planes                  [n_frame*PXL*(2*lu)]
*   @param: map_in  mapping for deinterleave                            [n_frame*PXL*lm]
*
*   @param: PSNR    computed PSNR between decoded frame and original frame [n_frame]
*   @param: imgr_out decoded clip                                       [n_frame*imgh*imgw]
*
*   @param: imgh, imgw, n_frame, lu, G
*
**/

void inter_intra(const char* str, int*** Y, double*** Ly_in, int*** map_in, const int &imgh, const int &imgw, const int &n_frame, const int &lu, int** G, double* PSNR , int*** img_out = NULL){

    // param initial
    const int lm = imgh*imgw;
    double* Lu = (double*) malloc(sizeof(double)*lu);

#if __EXIT_INFO__
    // EXIT
    int*** img_bp = new3d<int>(PXL,imgh,imgw);
    int*** img_bp_prev = new3d<int>(PXL,imgh,imgw);
#endif

    // frame buffer
    int** imgO;
    int** imgr = new2d<int>(imgh,imgw);
    int*** imgr_bp = new3d<int>(PXL,imgh,imgw);
    double** Ly;    //channel value
    int** map;      //interleaver map

    double** beta_s = new2d<double>(n_frame,PXL,0);
    double** beta_t = new2d<double>(n_frame,PXL,0);

    // pstate, pout
    const int Ns = pow(2,G_L-1);
    int ** pout = new2d<int>(Ns,4);
    int ** pstate = new2d<int>(Ns,2) ;
    double* Le1 = (double*)malloc(sizeof(double)*lu);
    double* Le2 = (double*)malloc(sizeof(double)*lu);

    for(int i = 0 ; i < lu ; ++i){
        Lu[i]  = 0;
//        Le1[i] = Le2[i] = -log(2);
//        Le1[i] = Le2[i] = 0.5;
    }
    computeLe(Lu,Le1,Le2,lu);

    trellis(G,G_N,G_L,Ns,pout,pstate);

    // frame buffer for previous frame
    int** imgO_prev;
    int** imgr_prev = new2d<int>(imgh,imgw);
    int*** imgr_bp_prev = new3d<int>(PXL,imgh,imgw);

    double** Ly_prev;    //channel value
    int** map_prev;      //interleaver map

    double** beta_s_prev = new2d<double>(n_frame,PXL,0);
    double** beta_t_prev = new2d<double>(n_frame,PXL,0);

    // frame buffer for previous 2 frame
    int*** imgr_bp_prev2 = new3d<int>(PXL,imgh,imgw);
    double** beta_t_prev2 = new2d<double>(n_frame,PXL,0);

    double** Le_s_prev2_inter = new2d<double>(PXL,lm,0); // source extrinsic information
    double** Le_s_prev2_intra = new2d<double>(PXL,lm,0); // source extrinsic information

    // ME
    const int mbSize = 8;
    const int me_range = 5;

    int** MV = new2d<int>(2,lm/(mbSize*mbSize),0);
    int** MV_prev = new2d<int>(2,lm/(mbSize*mbSize),0);
    double** Le_ref = new2d<double>(PXL,lm);

    // resetting buffer
    double** Lu_c = new2d<double>(PXL,lu);  //channel decoder output
    double** Le_c = new2d<double>(PXL,lm); //channel extrinsic information
    double** Lu_s = new2d<double>(PXL,lm); //source decoder output
    double** Le_s = new2d<double>(PXL,lm); // source extrinsic information
    double** Le_s_inter = new2d<double>(PXL,lm,0); // source extrinsic information
    double** Le_s_intra = new2d<double>(PXL,lm,0); // source extrinsic information

    double** Lu_c_prev = new2d<double>(PXL,lu);  //channel decoder output
    double** Le_c_prev = new2d<double>(PXL,lm); //channel extrinsic information
    double** Lu_s_prev = new2d<double>(PXL,lm); //source decoder output
    double** Le_s_prev = new2d<double>(PXL,lm); // source extrinsic information
    double** Le_s_prev_inter = new2d<double>(PXL,lm,0); // source extrinsic information
    double** Le_s_prev_intra = new2d<double>(PXL,lm,0); // source extrinsic information

    double*** imgr_soft_bp = new3d<double>(imgh,imgw,PXL);
    double*** imgr_soft_bp_prev = new3d<double>(imgh,imgw,PXL);


#if  __EXIT_INFO__
    // EXIT Ie_c, Ie_s
    double*** Ie_c = new3d<double>(n_frame,Niter,PXL);
    double*** Ie_c_prev = new3d<double>(n_frame,Niter,PXL);
    double*** Ie_s = new3d<double>(n_frame,Niter,PXL);
    double*** Ie_s_prev = new3d<double>(n_frame,Niter,PXL);
#endif
    // decoding [ start from the second frame ]
    for(int f = 1 ; f < n_frame ; ++f){
        printf("Decoding frame#%d\n",f+1);

        // initialize channel/source value

        Ly_prev = Ly_in[f-1];
        map_prev = map_in[f-1];
        imgO_prev = Y[f-1];

        Ly = Ly_in[f];
        map = map_in[f];
        imgO = Y[f];


#if __EXIT_INFO__
        img2bp_frame(imgO,imgh,imgw,img_bp);
        img2bp_frame(imgO_prev,imgh,imgw,img_bp_prev);
#endif

        for(int i = 0, j = 0  ; i < PXL ; ++i){
            for(j = 0 ; j < lm ; ++j){
                Le_c[i][j] = Lu_s[i][j] = Le_s[i][j] = Le_c_prev[i][j] = Le_s_prev[i][j] =  0;
                Lu_s_prev[i][j] = Le_s_prev2_intra[i][j] + Le_s_prev2_inter[i][j]*beta_t[f-1][i];
            }

            for(j = 0 ; j < lu ; ++j)
                Lu_c[i][j] = Lu_c_prev[i][j] = 0;
        }

//        motionEstES(imgO_prev,imgO,imgh,imgw,mbSize,me_range,MV);
//        motionEstES(imgO,imgO_prev,imgh,imgw,mbSize,me_range,MV_prev);

        // start iterative decode ~
        for(int iter = 0 ; iter < Niter ; ++iter){
            printf("iter #%d\n",iter+1);

            for(int i = 0 ; i < PXL ; ++i)
                for(int j = 0 ; j < lm ; ++j){
                    Le_s[i][j] = (Lu_s[i][j] - Le_c[i][j]);
                    Le_s_prev[i][j] = (Lu_s_prev[i][j] - Le_c_prev[i][j]);
                }

            // interleave
#if __STATUS__
            printf("interleave ...%lf\n",getCurrentTime());
#endif
            interleave(Le_s,map,lm);
            interleave(Le_s_prev,map_prev,lm);

            // BCJR decoding
#if __STATUS__
            printf("BCJR decoding ...%lf\n",getCurrentTime());
#endif
            for(int t_lvl = 0, i=0 ; t_lvl < PXL ; ++t_lvl){
                for(i = 0 ; i < lm ; ++i){
                    Lu[i] = Le_s[t_lvl][i];
//                    Le1[i] = -log(1+exp(Lu[i]));
//                    Le2[i] = (Le1[i] + Lu[i]);
                }

                computeLe(Lu,Le1,Le2,lm);

                BCJR_decoding(Ns, lu, 1, Ly[t_lvl], Le1, Le2, pstate, pout, Lu_c[t_lvl]);

                for(i = 0 ; i < lm ; ++i){
                    Lu[i] = Le_s_prev[t_lvl][i];
//                    Le1[i] = -log(1+exp(Lu[i]));
//                    Le2[i] = (Le1[i] + Lu[i]);
                }

                computeLe(Lu,Le1,Le2,lm);

                BCJR_decoding(Ns, lu, 1, Ly_prev[t_lvl], Le1, Le2, pstate, pout, Lu_c_prev[t_lvl]);
            }

            // sign detector for ME
#if __STATUS__
            printf("sign detecot for ME ...%lf\n",getCurrentTime());
#endif
            for(int i = 0, j=0, t_lvl=0,ii=0,jj=0 ; i <imgh ; ++i)
                for(j = 0 ; j < imgw ; ++j)
                    for(t_lvl=0 ; t_lvl < PXL ; ++t_lvl){
                        ii = map[t_lvl][j+i*imgw]/imgw;
                        jj = map[t_lvl][j+i*imgw]%imgw;
                        imgr_bp[t_lvl][ii][jj] = ((Lu_c[t_lvl][j+i*imgw]>=0)?1:0);

                        imgr_soft_bp[ii][jj][t_lvl] = exp(Lu_c[t_lvl][j+i*imgw]);

                        ii = map_prev[t_lvl][j+i*imgw]/imgw;
                        jj = map_prev[t_lvl][j+i*imgw]%imgw;
                        imgr_bp_prev[t_lvl][ii][jj] = ((Lu_c_prev[t_lvl][j+i*imgw]>=0)?1:0);

                        imgr_soft_bp_prev[ii][jj][t_lvl] = exp(Lu_c_prev[t_lvl][j+i*imgw]);
                    }

//            bin2dec_img(imgr_bp,imgh,imgw,imgr);
//            bin2dec_img(imgr_bp_prev,imgh,imgw,imgr_prev);

            // motion estimation
#if __STATUS__
            printf("Motion Estimation ...%lf\n",getCurrentTime());
#endif
//            motionEstES(imgr_prev,imgr,imgh,imgw,mbSize,me_range,MV);
//            motionEstES(imgr,imgr_prev,imgh,imgw,mbSize,me_range,MV_prev);

            motionEstES<double>(imgr_soft_bp_prev,imgr_soft_bp,imgh,imgw,mbSize,me_range,MV);
            motionEstES<double>(imgr_soft_bp,imgr_soft_bp_prev,imgh,imgw,mbSize,me_range,MV_prev);

            // generate extrinsic information
            for(int i = 0 ; i < PXL ; ++i)
                for(int j = 0 ; j < lm ; ++j){
                    Le_c[i][j] = (Lu_c[i][j] - Le_s[i][j]);
                    Le_c_prev[i][j] = (Lu_c_prev[i][j] - Le_s_prev[i][j]);
                }

            // deinterleave
#if __STATUS__
            printf("deinterleave ...%lf\n",getCurrentTime());
#endif
            deinterleave(Le_c,map,lm);
            deinterleave(Le_c_prev,map_prev,lm);

            for(int i = 0 ; i <imgh ; ++i)
                for(int j = 0 ; j < imgw ; ++j)
                    for(int t_lvl=0 ; t_lvl < PXL ; ++t_lvl){
                        imgr_bp[t_lvl][i][j] = ((Le_c[t_lvl][j+i*imgw]>=0)?1:0);
                        imgr_bp_prev[t_lvl][i][j] = ((Le_c_prev[t_lvl][j+i*imgw]>=0)?1:0);
                    }

#if __EXIT_INFO__
            // EXIT chart Ie_c
            mutual_information(Le_c,img_bp,lm,Ie_c[f][iter]);
            mutual_information(Le_c_prev,img_bp_prev,lm,Ie_c_prev[f][iter]);
#endif

            // MRF parameter estimation
#if __STATUS__
            printf("MRF parameter estimation ...%lf\n",getCurrentTime());
#endif
            // joint estimation
            intra_inter_beta_estimation(imgr_bp,imgr_bp_prev,MV_prev,beta_s[f],beta_t[f],imgh,imgw,mbSize);
            if(f==1)
                intra_inter_beta_estimation(imgr_bp_prev,imgr_bp,MV,beta_s_prev[f],beta_t_prev[f],imgh,imgw,mbSize);
            else
                intra_inter2_beta_estimation(imgr_bp,imgr_bp_prev,imgr_bp_prev2,MV,beta_s_prev[f],beta_t_prev[f],beta_t_prev2[f],imgh,imgw,mbSize);

#if __BETA__
            printf("frame#%d iter#%d,\n",f+1,iter+1);
            for(int i = 0 ; i < PXL ; ++i)
                printf("bp %d beta_s = %lf , beta_s_prev = %lf \n",i,beta_s[f][i],beta_s_prev[f][i]);
            for(int i = 0 ; i < PXL ; ++i)
                printf("bp %d beta_t = %lf , beta_t_prev = %lf , beta_t_prev2 = %lf\n",i,beta_t[f][i],beta_t_prev[f][i],beta_t_prev2[f][i]);
#endif

            if(iter == Niter-1){
                for(int i = 0 ;  i < PXL ; ++i)
                        motionComp(imgr_bp_prev[i],MV_prev,imgh,imgw,mbSize,imgr_bp_prev2[i]);
            }

            // MRF decoding
#if __STATUS__
            printf("MRF decoding ... %lf\n",getCurrentTime());
#endif
            motionComp(Le_c_prev,MV_prev,imgh,imgw,mbSize,Le_ref);
            mrf_siso_inter(Le_c,Le_ref,beta_t[f],imgh,imgw,Le_s_inter,0);
            mrf_siso_intra(Le_c,beta_s[f],imgh,imgw,Le_s_intra,0);

            motionComp(Le_c,MV,imgh,imgw,mbSize,Le_ref);
            mrf_siso_inter(Le_c_prev,Le_ref,beta_t_prev[f],imgh,imgw,Le_s_prev_inter,0);
            mrf_siso_intra(Le_c_prev,beta_s_prev[f],imgh,imgw,Le_s_prev_intra,0);

            for(int t_lvl = 0, i = 0 ; t_lvl < PXL ; ++t_lvl)
                for(i=0;i<lm;++i){
                    Lu_s[t_lvl][i] = (Le_s_inter[t_lvl][i] + Le_s_intra[t_lvl][i]);
                    Lu_s_prev[t_lvl][i] = (Le_s_prev_inter[t_lvl][i] + Le_s_prev_intra[t_lvl][i] + Le_s_prev2_inter[t_lvl][i]*beta_t_prev2[f][t_lvl]);
#if !__EXIT_INFO__
                    Lu_s[t_lvl][i]+=Le_c[t_lvl][i];
                    Lu_s_prev[t_lvl][i]+=Le_c_prev[t_lvl][i];
#endif
                }

#if __EXIT_INFO__
            // EXIT CHART Ie_s
            mutual_information(Lu_s,img_bp,lm,Ie_s[f][iter]);
            mutual_information(Lu_s_prev,img_bp_prev,lm,Ie_s_prev[f][iter]);

            for(int t_lvl = 0, i = 0 ; t_lvl < PXL ; ++t_lvl)
                for(i=0;i<lm;++i){
                    Lu_s[t_lvl][i]+=Le_c[t_lvl][i];
                    Lu_s_prev[t_lvl][i]+=Le_c_prev[t_lvl][i];
                }
#endif

            if(iter==Niter-1){
                for(int t_lvl = 0, i = 0 ; t_lvl < PXL ; ++t_lvl)
                    for(i=0;i<lm;++i){
                        Le_s_prev2_inter[t_lvl][i] = Le_s_inter[t_lvl][i]/beta_t[f][t_lvl];
                        Le_s_prev2_intra[t_lvl][i] = Le_s_intra[t_lvl][i];
                    }
            }


            // recover to image
            for(int i = 0 ; i <imgh ; ++i){
                for(int j = 0 ; j < imgw ; ++j){
                    for(int t_lvl=0 ; t_lvl < PXL ; ++t_lvl){
                        imgr_bp[t_lvl][i][j] = ((Lu_s[t_lvl][j+i*imgw]>=0)?1:0);
                        imgr_bp_prev[t_lvl][i][j] = ((Lu_s_prev[t_lvl][j+i*imgw]>=0)?1:0);
                    }
                }
            }

            // construct image
            bin2dec_img(imgr_bp,imgh,imgw,imgr);
            bin2dec_img(imgr_bp_prev,imgh,imgw,imgr_prev);

            // compute PSNR
            PSNR[f] = computePSNR(imgr,imgO,imgh,imgw);
            PSNR[f-1] = computePSNR(imgr_prev,imgO_prev,imgh,imgw);
#if __PSNR__
            printf("%s iter#%d frame#%d PSNR = %lf, frame#%d PSNR = %lf\n",str,iter+1,f,PSNR[f-1],f+1,PSNR[f]);
#endif

        }

        // imgr output
        if(img_out!=NULL)
            for(int i = 0, j = 0; i <imgh ; ++i)
                for(j=0 ; j < imgw ; ++j)
                    img_out[f-1][i][j] = imgr_prev[i][j];


    }

    // imgr output for last frame
    if(img_out!=NULL)
        for(int i = 0, j = 0; i <imgh ; ++i)
            for(j=0 ; j < imgw ; ++j)
                img_out[n_frame-1][i][j] = imgr[i][j];

#if __EXIT_INFO__
    write_exit_info(Ie_c,Ie_c_prev,Ie_s,Ie_s_prev);

    delete3d<double>(Ie_c);
    delete3d<double>(Ie_c_prev);
    delete3d<double>(Ie_s);
    delete3d<double>(Ie_s_prev);

#endif
    // free memory
    delete3d<int>(imgr_bp);
    delete3d<int>(imgr_bp_prev);
    delete3d<int>(imgr_bp_prev2);
    delete3d<double>(imgr_soft_bp);
    delete3d<double>(imgr_soft_bp_prev);

    delete2d<int>(imgr);
    delete2d<int>(imgr_prev);
    delete2d<int>(G);
    delete2d<int>(pstate);
    delete2d<int>(pout);
    delete2d<double>(beta_t);
    delete2d<double>(beta_s);
    delete2d<double>(beta_t_prev);
    delete2d<double>(beta_s_prev);
    delete2d<double>(beta_t_prev2);

    free(Lu);
    free(Le1);
    free(Le2);

    delete2d<double>(Lu_c);
    delete2d<double>(Lu_s);
    delete2d<double>(Le_c);
    delete2d<double>(Le_s);
    delete2d<double>(Le_s_inter);
    delete2d<double>(Le_s_intra);
    delete2d<double>(Lu_c_prev);
    delete2d<double>(Lu_s_prev);
    delete2d<double>(Le_c_prev);
    delete2d<double>(Le_s_prev);
    delete2d<double>(Le_s_prev_inter);
    delete2d<double>(Le_s_prev_intra);
    delete2d<double>(Le_s_prev2_inter);
    delete2d<double>(Le_s_prev2_intra);
    delete2d<double>(Le_ref);

    delete2d<int>(MV);
    delete2d<int>(MV_prev);
}


#endif // __ALGS_INTER_INTRA_H
