#ifndef __ALGS_INTER_INTRA_H
#define __ALGS_INTER_INTRA_H

#include <stdio.h>
#include <stdlib.h>
#include "frame.h"
#include "mrf_decoder_utils.h"
#include "uep_predict_inter_intra.h"

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

void inter_intra_system(const char* str, FILE* fptr, const int &imgh, const int &imgw, const int &n_frame,int** G,const double &snr, double* PSNR ,int weight_type = 0, int*** img_out = NULL){

    Frame frameMgr[2] = {Frame(imgh,imgw,0,0), Frame(imgh,imgw,0,1)};
    Frame* frame ;
    Frame* frame_prev;
    Frame* frame3;
    frameMgr[0].encode_info(snr,G);
    frameMgr[1].encode_info(snr,G);

    double* ber = MALLOC(double, PXL);

    // param initial
    const int lm = imgh*imgw;
    const int lu = lm + 2;

    const double EbN0 = pow(10,snr/10);
    double* weights = MALLOC(double,PXL);
    for(int i = 0 ; i < PXL; ++i)
        weights[i] = 1;

    // frame buffer
    Pixel** imgr = new2d<Pixel>(imgh,imgw);
    int8*** imgr_bp = new3d<int8>(PXL,imgh,imgw);

    double** LyMgr[2] = {new2d<double>(PXL,2*lu), new2d<double>(PXL,2*lu)};
    int** MapMgr[2] = {new2d<int>(PXL,lm), new2d<int>(PXL,lm)};

    double** Ly;    //channel value
    int** map;      //interleaver map

    double* beta_s = MALLOC(double,PXL);
    double* beta_t = MALLOC(double,PXL);

    // pstate, pout
    double* Lu = MALLOC(double,lu);
    double* Le1 = MALLOC(double,lu);
    double* Le2 = MALLOC(double,lu);

    for(int i = 0 ; i < lu ; ++i)
        Lu[i]  = 0;

    computeLe(Lu,Le1,Le2,lu);

    // frame buffer for previous frame
    Pixel** imgr_prev = new2d<Pixel>(imgh,imgw);
    int8*** imgr_bp_prev = new3d<int8>(PXL,imgh,imgw);

    double** Ly_prev;    //channel value
    int** map_prev;      //interleaver map

    double* beta_s_prev = MALLOC(double,PXL);
    double* beta_t_prev = MALLOC(double,PXL);

    // frame buffer for previous 2 frame
    int8*** imgr_bp_prev2 = new3d<int8>(PXL,imgh,imgw);
    double* beta_t_prev2 = MALLOC(double,PXL);

    double** Le_s_prev2_inter = new2d<double>(PXL,lm,0); // source extrinsic information
    double** Le_s_prev2_intra = new2d<double>(PXL,lm,0); // source extrinsic information

    // ME
    const int mbSize = 8;
    const int me_range = 5;

    int** MV = new2d<int>(2,lm/(mbSize*mbSize),0);
    int** MV_prev = new2d<int>(2,lm/(mbSize*mbSize),0);
    double** Le_ref = new2d<double>(PXL,lm);


    for(int i = 0 ; i < PXL ; ++i)
        beta_s[i] = beta_t[i] = beta_s_prev[i] = beta_t_prev[i] = beta_t_prev2[i] = 0;

    // assigning first frame
    if(weight_type){
        frame3 = new Frame(imgh,imgw,0,2);
        frame3->encode_info(snr,G);
        frameMgr[0].read(fptr);
        frameMgr[1].read(fptr);
        weight_predict_Inter_minMSE(frameMgr[0].Y,frameMgr[1].Y,imgh,imgw,weights,EbN0);
        frameMgr[0].encode(LyMgr[0],MapMgr[0],weights);

    }else
        frameMgr[0].next(fptr,LyMgr[0],MapMgr[0],weights);

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

    // decoding [ start from the second frame ]
    for(int f = 1 ; f < n_frame ; ++f){

    // initialize channel/source value
#if __PROGRESS__
        printf("Encoding frame#%d\n",f+1);
#endif
        Ly_prev = LyMgr[(f+1)%2];
        Ly = LyMgr[f%2];
        map_prev = MapMgr[(f+1)%2];
        map = MapMgr[f%2];

        if(weight_type){
            if(f==n_frame-1){
                switch(f%3){
                case 0:
                    weight_predict_Inter_Intra_minMSE(frameMgr[0].Y,frame3->Y,imgh,imgw,weights,EbN0);
                    frameMgr[0].encode(LyMgr[f%2],MapMgr[f%2],weights);
                    frame_prev = frame3;
                    frame = &frameMgr[0];
                    break;
                case 1:
                    weight_predict_Inter_Intra_minMSE(frameMgr[1].Y,frameMgr[0].Y,imgh,imgw,weights,EbN0);
                    frameMgr[1].encode(LyMgr[f%2],MapMgr[f%2],weights);
                    frame_prev = &frameMgr[0];
                    frame = &frameMgr[1];
                    break;
                case 2:
                    weight_predict_Inter_Intra_minMSE(frame3->Y,frameMgr[1].Y,imgh,imgw,weights,EbN0);
                    frame3->encode(LyMgr[f%2],MapMgr[f%2],weights);
                    frame_prev = &frameMgr[1];
                    frame = frame3;
                    break;
                }
            }else{
                switch(f%3){
                case 0:
                    frameMgr[1].read(fptr);
                    weight_predict_Inter_Intra_minMSE(frameMgr[0].Y,frame3->Y,frameMgr[1].Y,imgh,imgw,weights,EbN0);
                    frameMgr[0].encode(LyMgr[f%2],MapMgr[f%2],weights);
                    frame_prev = frame3;
                    frame = &frameMgr[0];
                    break;
                case 1:
                    frame3->read(fptr);
                    weight_predict_Inter_Intra_minMSE(frameMgr[1].Y,frameMgr[0].Y,frame3->Y,imgh,imgw,weights,EbN0);
                    frameMgr[1].encode(LyMgr[f%2],MapMgr[f%2],weights);
                    frame_prev = &frameMgr[0];
                    frame = &frameMgr[1];
                    break;
                case 2:
                    frameMgr[0].read(fptr);
                    weight_predict_Inter_Intra_minMSE(frame3->Y,frameMgr[1].Y,frameMgr[0].Y,imgh,imgw,weights,EbN0);
                    frame3->encode(LyMgr[f%2],MapMgr[f%2],weights);
                    frame_prev = &frameMgr[1];
                    frame = frame3;
                    break;
                }
            }
//            for(int i = 0 ; i < PXL ; ++i){
//                printf("%.3f, ",weights[i]);
//            }
//            printf("\n");

        }else{
            frameMgr[f%2].next(fptr,Ly,map,weights);
            frame_prev = &frameMgr[(f+1)%2];
            frame = &frameMgr[f%2];
        }

#if __PROGRESS__
        printf("Decoding frame#%d\n",f+1);
#endif

        for(int i = 0, j = 0  ; i < PXL ; ++i){
            for(j = 0 ; j < lm ; ++j){
                Le_c[i][j] = Lu_s[i][j] = Le_s[i][j] = Le_c_prev[i][j] = Le_s_prev[i][j] =  0;
                Lu_s_prev[i][j] = Le_s_prev2_intra[i][j] + Le_s_prev2_inter[i][j]*beta_t[i];
            }

            for(j = 0 ; j < lu ; ++j)
                Lu_c[i][j] = Lu_c_prev[i][j] = 0;
        }

        // start iterative decode ~
        for(int iter = 0 ; iter < Niter ; ++iter){
#if __PROGRESS__
            printf("iter #%d\n",iter+1);
#endif

#if __STATUS__
        printf("interleave ...%lf\n",getCurrentTime());
#endif
            // interleave
            for(int i = 0 ; i < PXL ; ++i)
                for(int j = 0 ; j < lm ; ++j){
                    Le_s[i][j] = (Lu_s[i][map[i][j]] - Le_c[i][map[i][j]]);
                    Le_s_prev[i][j] = (Lu_s_prev[i][map_prev[i][j]] - Le_c_prev[i][map_prev[i][j]]);
                }

            // BCJR decoding
#if __STATUS__
            printf("BCJR decoding ...%lf\n",getCurrentTime());
#endif
            for(int t_lvl = 0, i=0 ; t_lvl < PXL ; ++t_lvl){
                for(i = 0 ; i < lm ; ++i)
                    Lu[i] = Le_s[t_lvl][i];

                computeLe(Lu,Le1,Le2,lm);

                BCJR_decoding(lu, 1, Ly[t_lvl], Le1, Le2, Lu_c[t_lvl]);

                for(i = 0 ; i < lm ; ++i)
                    Lu[i] = Le_s_prev[t_lvl][i];

                computeLe(Lu,Le1,Le2,lm);
                BCJR_decoding(lu, 1, Ly_prev[t_lvl], Le1, Le2, Lu_c_prev[t_lvl]);
            }

            // sign detector for ME
#if __STATUS__
            printf("sign detecot for ME ...%lf\n",getCurrentTime());
#endif

            Lu2dec_img(Lu_c,lm,imgr,map);
            Lu2dec_img(Lu_c_prev,lm,imgr_prev,map_prev);

//            for(int i = 0, j=0, t_lvl=0,ii=0,jj=0 ; i <imgh ; ++i)
//                for(j = 0 ; j < imgw ; ++j)
//                    for(t_lvl=0 ; t_lvl < PXL ; ++t_lvl){
//                        ii = map[t_lvl][j+i*imgw]/imgw;
//                        jj = map[t_lvl][j+i*imgw]%imgw;
//                        imgr_bp[t_lvl][ii][jj] = ((Lu_c[t_lvl][j+i*imgw]>=0)?1:0);
//
////                        imgr_soft_bp[ii][jj][t_lvl] = exp(Lu_c[t_lvl][j+i*imgw]);
//
//                        ii = map_prev[t_lvl][j+i*imgw]/imgw;
//                        jj = map_prev[t_lvl][j+i*imgw]%imgw;
//                        imgr_bp_prev[t_lvl][ii][jj] = ((Lu_c_prev[t_lvl][j+i*imgw]>=0)?1:0);
//
////                        imgr_soft_bp_prev[ii][jj][t_lvl] = exp(Lu_c_prev[t_lvl][j+i*imgw]);
//                    }
//
//            bin2dec_img(imgr_bp,imgh,imgw,imgr);
//            bin2dec_img(imgr_bp_prev,imgh,imgw,imgr_prev);

            // motion estimation
#if __STATUS__
            printf("Motion Estimation ...%lf\n",getCurrentTime());
#endif
            motionEstES(imgr_prev,imgr,imgh,imgw,mbSize,me_range,MV);
            motionEstES(imgr,imgr_prev,imgh,imgw,mbSize,me_range,MV_prev);

//            motionEstES<double>(imgr_soft_bp_prev,imgr_soft_bp,imgh,imgw,mbSize,me_range,MV);
//            motionEstES<double>(imgr_soft_bp,imgr_soft_bp_prev,imgh,imgw,mbSize,me_range,MV_prev);

#if __STATUS__
            printf("deinterleave ...%lf\n",getCurrentTime());
#endif
            // generate extrinsic information
            for(int i = 0 ; i < PXL ; ++i)
                for(int j = 0 ; j < lm ; ++j){
                    Le_c[i][map[i][j]] = (Lu_c[i][j] - Le_s[i][j]);
                    Le_c_prev[i][map_prev[i][j]] = (Lu_c_prev[i][j] - Le_s_prev[i][j]);
                }

            for(int t_lvl=0, i=0 ; t_lvl < PXL ; ++t_lvl)
                for(i=0 ; i < lm ; ++i){
                    imgr_bp[t_lvl][0][i] = ((Le_c[t_lvl][i]>=0)?1:0);
                    imgr_bp_prev[t_lvl][0][i] = ((Le_c_prev[t_lvl][i]>=0)?1:0);
                }


            // MRF parameter estimation
#if __STATUS__
            printf("MRF parameter estimation ...%lf\n",getCurrentTime());
#endif
            // joint estimation
            intra_inter_beta_estimation(imgr_bp,imgr_bp_prev,MV_prev,beta_s,beta_t,imgh,imgw,mbSize);
            if(f==1)
                intra_inter_beta_estimation(imgr_bp_prev,imgr_bp,MV,beta_s_prev,beta_t_prev,imgh,imgw,mbSize);
            else
                intra_inter2_beta_estimation(imgr_bp,imgr_bp_prev,imgr_bp_prev2,MV,beta_s_prev,beta_t_prev,beta_t_prev2,imgh,imgw,mbSize);

#if __BETA__
            printf("frame#%d iter#%d,\n",f+1,iter+1);
            for(int i = 0 ; i < PXL ; ++i)
                printf("bp %d beta_s = %lf , beta_s_prev = %lf \n",i,beta_s[i],beta_s_prev[i]);
            for(int i = 0 ; i < PXL ; ++i)
                printf("bp %d beta_t = %lf , beta_t_prev = %lf , beta_t_prev2 = %lf\n",i,beta_t[i],beta_t_prev[i],beta_t_prev2[i]);
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
            mrf_siso_inter(Le_c,Le_ref,beta_t,lm,Le_s_inter,0);
            mrf_siso_intra(Le_c,beta_s,imgh,imgw,Le_s_intra,0);

            motionComp(Le_c,MV,imgh,imgw,mbSize,Le_ref);
            mrf_siso_inter(Le_c_prev,Le_ref,beta_t_prev,lm,Le_s_prev_inter,0);
            mrf_siso_intra(Le_c_prev,beta_s_prev,imgh,imgw,Le_s_prev_intra,0);

            for(int t_lvl = 0, i = 0 ; t_lvl < PXL ; ++t_lvl)
                for(i=0;i<lm;++i){
                    Lu_s[t_lvl][i] = (Le_c[t_lvl][i] + Le_s_inter[t_lvl][i] + Le_s_intra[t_lvl][i]);
                    Lu_s_prev[t_lvl][i] = (Le_c_prev[t_lvl][i] + Le_s_prev_inter[t_lvl][i] + Le_s_prev_intra[t_lvl][i] + Le_s_prev2_inter[t_lvl][i]*beta_t_prev2[t_lvl]);
                }

            if(iter==Niter-1){
                for(int t_lvl = 0, i = 0 ; t_lvl < PXL ; ++t_lvl)
                    for(i=0;i<lm;++i){
                        Le_s_prev2_inter[t_lvl][i] = Le_s_inter[t_lvl][i]/beta_t[t_lvl];
                        Le_s_prev2_intra[t_lvl][i] = Le_s_intra[t_lvl][i];
                    }
            }

            // recover to image
            Lu2dec_img(Lu_s,lm,imgr);
            Lu2dec_img(Lu_s_prev,lm,imgr_prev);

            // compute PSNR
            PSNR[f] = frame->psnr(imgr);
            PSNR[f-1] = frame_prev->psnr(imgr_prev);
#if __PSNR__
            printf("%s iter#%d frame#%d PSNR = %lf, frame#%d PSNR = %lf\n",str,iter+1,f,PSNR[f-1],f+1,PSNR[f]);
#endif

        }

        computeBER(imgr_prev,frame_prev->Y,lm,ber);
        printf("%lf\n",PSNR[f-1]);

        if(f==(n_frame-1)){
            computeBER(imgr,frame->Y,lm,ber);
            printf("%lf\n",PSNR[f]);
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


    // free memory
    delete3d<int8>(imgr_bp);
    delete3d<int8>(imgr_bp_prev);
    delete3d<int8>(imgr_bp_prev2);
//    delete3d<double>(imgr_soft_bp);
//    delete3d<double>(imgr_soft_bp_prev);

    delete2d<Pixel>(imgr);
    delete2d<Pixel>(imgr_prev);
    DELETE(beta_t);
    DELETE(beta_s);
    DELETE(beta_t_prev);
    DELETE(beta_s_prev);
    DELETE(beta_t_prev2);

    DELETE(Lu);
    DELETE(Le1);
    DELETE(Le2);
    DELETE(weights);
    DELETE(ber);

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

    if(weight_type)
        delete frame3;
}


#endif // __ALGS_INTER_INTRA_H
