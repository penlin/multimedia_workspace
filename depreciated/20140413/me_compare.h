#ifndef __ME_COMPARE_H
#define __ME_COMPARE_H

#include <stdio.h>
#include <stdlib.h>
#include "data_alloc.h"
#include "build_value.h"
#include "motion_estimation_utils.h"
#include "motion_estimation_soft_utils.h"


/**
*   @Penlin: comparison of the motion estimations
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

void compare(const char* str, int*** Y, double*** Ly_in, int*** map_in, const int &imgh, const int &imgw, const int &n_frame, const int &lu, int** G, double** PSNR){

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
        Lu[i]  = 0;
//        Le1[i] = Le2[i] = -log(2);
    }
    computeLe(Lu,Le1,Le2,lu);

    trellis(G,G_N,G_L,Ns,pout,pstate);


    // frame buffer for previous frame
    int** imgO_prev;
    int** imgr_prev = new2d<int>(imgh,imgw);
    int*** imgr_bp_prev = new3d<int>(PXL,imgh,imgw);

    double** Ly_prev;    //channel value
    int** map_prev;      //interleaver map

    // ME
    const int mbSize = 8;
    const int me_range = 5;

    int** MV_ori = new2d<int>(2,lm/(mbSize*mbSize),0);
    int** MV_ori_prev = new2d<int>(2,lm/(mbSize*mbSize),0);

    int** MV_hard = new2d<int>(2,lm/(mbSize*mbSize),0);
    int** MV_hard_prev = new2d<int>(2,lm/(mbSize*mbSize),0);

    int** MV = new2d<int>(2,lm/(mbSize*mbSize),0);
    int** MV_prev = new2d<int>(2,lm/(mbSize*mbSize),0);

    double** Le_ref = new2d<double>(PXL,lm);
    double** Le_hard = new2d<double>(PXL,lm);
    double** Le_soft = new2d<double>(PXL,lm);

    // assigning first frame
    Ly = Ly_in[0];
    map = map_in[0];
    imgO = Y[0];

    // resetting buffer
    double** Lu_c = new2d<double>(PXL,lu);  //channel decoder output

    double** Lu_c_prev = new2d<double>(PXL,lu);  //channel decoder output


    int*** imgr_bp_ref = new3d<int>(PXL,imgh,imgw);
    int*** imgr_bp_hard = new3d<int>(PXL,imgh,imgw);
    int*** imgr_bp_soft = new3d<int>(PXL,imgh,imgw);

    double*** imgr_soft_bp = new3d<double>(imgh,imgw,PXL);
    double*** imgr_soft_bp_prev = new3d<double>(imgh,imgw,PXL);

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

        motionEstES(imgO_prev,imgO,imgh,imgw,mbSize,me_range,MV_ori);
//        motionEstES(imgO,imgO_prev,imgh,imgw,mbSize,me_range,MV_ori_prev);
        motionComp(imgO,MV_ori,imgh,imgw,mbSize,imgr);
        PSNR[0][f] = computePSNR(imgr,imgO_prev,imgh,imgw);


        // BCJR decoding
        for(int t_lvl = 0, i=0 ; t_lvl < PXL ; ++t_lvl){
            BCJR_decoding(Ns, lu, 1, Ly[t_lvl], Le1, Le2, pstate, pout, Lu_c[t_lvl]);
            BCJR_decoding(Ns, lu, 1, Ly_prev[t_lvl], Le1, Le2, pstate, pout, Lu_c_prev[t_lvl]);
        }


        // sign detector for ME
//        for(int i = 0 , j = 0 ; i < imgh ; ++i)
//            for(j=0 ; j<imgw ; ++j)
//                imgr[i][j] = imgr_prev[i][j] = 0;


        for(int i = 0, j=0, t_lvl=0,ii=0,jj=0 ; i <imgh ; ++i)
            for(j = 0 ; j < imgw ; ++j)
                for(t_lvl=0 ; t_lvl < PXL ; ++t_lvl){
                    ii = map[t_lvl][j+i*imgw]/imgw;
                    jj = map[t_lvl][j+i*imgw]%imgw;
                    imgr_bp[t_lvl][ii][jj] = ((Lu_c[t_lvl][j+i*imgw]>=0)?1:0);

                    imgr_soft_bp[ii][jj][t_lvl] = Lu_c[t_lvl][j+i*imgw];
//                    imgr[ii][jj]+=llr_bp_to_img(Lu_c[t_lvl][j+i*imgw],t_lvl);

                    ii = map_prev[t_lvl][j+i*imgw]/imgw;
                    jj = map_prev[t_lvl][j+i*imgw]%imgw;
                    imgr_bp_prev[t_lvl][ii][jj] = ((Lu_c_prev[t_lvl][j+i*imgw]>=0)?1:0);

//                    imgr_prev[ii][jj]+=llr_bp_to_img(Lu_c_prev[t_lvl][j+i*imgw],t_lvl);
                    imgr_soft_bp_prev[ii][jj][t_lvl] = Lu_c_prev[t_lvl][j+i*imgw];
                }

//            bin2dec_img(imgr_bp,imgh,imgw,imgr);
//            bin2dec_img(imgr_bp_prev,imgh,imgw,imgr_prev);

        // motion estimation
//        motionEstES(imgr_prev,imgr,imgh,imgw,mbSize,me_range,MV);
//        motionEstES(imgr,imgr_prev,imgh,imgw,mbSize,me_range,MV_prev);

        bin2dec_img(imgr_bp,imgh,imgw,imgr);
        bin2dec_img(imgr_bp_prev,imgh,imgw,imgr_prev);

        motionEstES(imgr_prev,imgr,imgh,imgw,mbSize,me_range,MV_hard);

        motionEstES<double>(imgr_soft_bp_prev,imgr_soft_bp,imgh,imgw,mbSize,me_range,MV);

        // deinterleave
        deinterleave(Lu_c,map,lm);
//        deinterleave(Lu_c_prev,map_prev,lm);

//        motionComp(Lu_c_prev,MV_prev,imgh,imgw,mbSize,Le_ref);

/*
        motionComp(Lu_c,MV_ori,imgh,imgw,mbSize,Le_ref);
        motionComp(Lu_c,MV_hard,imgh,imgw,mbSize,Le_hard);
        motionComp(Lu_c,MV,imgh,imgw,mbSize,Le_soft);

        // recover to image
        for(int i = 0 ; i <imgh ; ++i){
            for(int j = 0 ; j < imgw ; ++j){
                for(int t_lvl=0 ; t_lvl < PXL ; ++t_lvl){
                    imgr_bp_ref[t_lvl][i][j] = ((Le_ref[t_lvl][j+i*imgw]>=0)?1:0);
                    imgr_bp_hard[t_lvl][i][j] = ((Le_hard[t_lvl][j+i*imgw]>=0)?1:0);
                    imgr_bp_soft[t_lvl][i][j] = ((Le_soft[t_lvl][j+i*imgw]>=0)?1:0);
                }
            }
        }

        // construct image
        bin2dec_img(imgr_bp_ref,imgh,imgw,imgr);
        bin2dec_img(imgr_bp_hard,imgh,imgw,imgr_prev);

        // compute PSNR
        PSNR[0][f] = computePSNR(imgr_prev,imgr,imgh,imgw);

        bin2dec_img(imgr_bp_soft,imgh,imgw,imgr_prev);
        PSNR[1][f] = computePSNR(imgr_prev,imgr,imgh,imgw);
//        PSNR[f-1] = computePSNR(imgr_prev,imgO_prev,imgh,imgw);

        bin2dec_img(imgr_bp_hard,imgh,imgw,imgr);
        PSNR[2][f] = computePSNR(imgr_prev,imgr,imgh,imgw);
*/

        motionComp(imgO,MV_ori,imgh,imgw,mbSize,imgr);
        motionComp(imgO,MV_hard,imgh,imgw,mbSize,imgr_prev);
        PSNR[1][f] = computePSNR(imgr_prev,imgr,imgh,imgw);

        motionComp(imgO,MV,imgh,imgw,mbSize,imgr_prev);
        PSNR[2][f] = computePSNR(imgr_prev,imgr,imgh,imgw);


    }


    // free memory
    delete3d<int>(imgr_bp_ref);
    delete3d<int>(imgr_bp_hard);
    delete3d<int>(imgr_bp_soft);
    delete3d<double>(imgr_soft_bp);
    delete3d<double>(imgr_soft_bp_prev);

    delete3d<int>(imgr_bp);
    delete3d<int>(imgr_bp_prev);
    delete2d<int>(imgr);
    delete2d<int>(imgr_prev);
    delete2d<int>(G);
    delete2d<int>(pstate);
    delete2d<int>(pout);

    free(Lu);
    free(Le1);
    free(Le2);

    delete2d<double>(Lu_c);
    delete2d<double>(Lu_c_prev);
    delete2d<double>(Le_ref);
    delete2d<double>(Le_hard);
    delete2d<double>(Le_soft);

    delete2d<int>(MV);
    delete2d<int>(MV_prev);
    delete2d<int>(MV_hard);
    delete2d<int>(MV_hard_prev);
    delete2d<int>(MV_ori);
    delete2d<int>(MV_ori_prev);

}


#endif // __ALGS_INTER_H
