#ifndef __ALGS_UEP_H
#define __ALGS_UEP_H
#include <omp.h>
#include <stdio.h>
#include <stdlib.h>
#include "frame.h"
#include "mrf_decoder_utils.h"
#include "uep_predict_intra.h"
#include "uep_predict_inter.h"
#include "uep_predict_inter_intra.h"

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

int flag_output = 1;

void yuv_write(FILE* fptr, Pixel** Y, const int& height, const int& width){
    if(fptr==NULL){
        printf("FILE can't be read");
        return;
    }
    int lm = height*width;
    unsigned char* buffer = (unsigned char*)malloc(sizeof(unsigned char)*lm);

    for(int i = 0, j = 0 ; i < height ; ++i){
        for(j = 0 ; j < width ; ++j)
            buffer[j+i*width] = Y[i][j];
    }
    fwrite(buffer,1,lm,fptr);

    DELETE(buffer);
}


void intra_system(const char* str, FILE* fptr, const int &imgh, const int &imgw, const int &n_frame,int** G, const double &snr, double* PSNR , int repredict = 0, int*** img_out = NULL){


    FILE* fout = NULL;
    if(flag_output){
        char buf[50];
        sprintf(buf,"%s%s_snr%d.yuv",__SEQ_DIR,str,(int)snr);
        fout = fopen(buf,"w+b");
    }
    Frame* frame = new Frame(imgh,imgw,0,0);
    frame->encode_info(snr,G);
    const double EbN0 = pow(10,snr/10);
    // param initial
    const int lm = imgh*imgw;

    double* weights = MALLOC(double,PXL);
    for(int i = 0 ; i < PXL ; ++i)
        weights[i] = 1;

    double sTime = 0.0;
    double timeUEP = 0.0;

    // decoding
    for(int f = 0 ; f < n_frame ; ++f){
#if __PROGRESS__
        printf("Processing frame#%d\n",f+1);
#endif

        frame->read(fptr);

        sTime = getCurrentTime();
        weight_predict_Intra_minMSE(frame->img_bp,imgh,imgw,weights,EbN0);
        timeUEP += getCurrentTime()-sTime;

        if(flag_output)
            yuv_write(fout,frame->Y,imgh,imgw);

    }
    printf("average time consumption: %f s\n",timeUEP/n_frame);
    DELETE(weights);
    delete frame;
    if(flag_output)
        fclose(fout);
}


void inter_system(const char* str, FILE* fptr, const int &imgh, const int &imgw, const int &n_frame,int** G, const double &snr, double* PSNR , int repredict = 0, int*** img_out = NULL){


    FILE* fout = NULL;
    if(flag_output){
        char buf[50];
        sprintf(buf,"%s%s_snr%d.yuv",__SEQ_DIR,str,(int)snr);
        fout = fopen(buf,"w+b");
    }
    Frame* frame = new Frame(imgh,imgw,0,0);
    frame->encode_info(snr,G);
    const double EbN0 = pow(10,snr/10);
    // param initial
    const int lm = imgh*imgw;

    double* weights = MALLOC(double,PXL);
    for(int i = 0 ; i < PXL ; ++i)
        weights[i] = 1;

    double sTime = 0.0;
    double timeUEP = 0.0;

    // decoding
    for(int f = 0 ; f < n_frame ; ++f){
#if __PROGRESS__
        printf("Processing frame#%d\n",f+1);
#endif

        frame->read(fptr);

        sTime = getCurrentTime();
        weight_predict_Intra_minMSE(frame->img_bp,imgh,imgw,weights,EbN0);
        timeUEP += getCurrentTime()-sTime;

        if(flag_output)
            yuv_write(fout,frame->Y,imgh,imgw);

    }
    printf("average time consumption: %f s\n",timeUEP/n_frame);
    DELETE(weights);
    delete frame;
    if(flag_output)
        fclose(fout);
}


#endif // __ALGS_UEP_H
