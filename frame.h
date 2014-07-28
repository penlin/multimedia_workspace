#ifndef __FRAME_H
#define __FRAME_H

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>
#include "rsc_encode.h"
#include "utils.h"
#include "data_alloc.h"
#include "channel_code_utils.h"
#include "image_process_utils.h"
#include "io_utils.h"

const int TYPE_Y = 0;
const int TYPE_YUV_420 = 1;
const int TYPE_YUV_444 = 2;

class Frame {

public:

    Frame(int h, int w, int t , int tag){
        height = h ;
        width = w;
        lm = h*w;
        lu = lm +2 ;
        type = t;
        frame_tag = tag;

        Y = new2d<Pixel>(h,w,0);
        img_bp = new3d<int8>(PXL,height,width,0);
        x = MALLOC(int8,2*lu);
        if(type == TYPE_YUV_420){
            U = new2d<Pixel>(h/2,w/2,0);
            V = new2d<Pixel>(h/2,w/2,0);
        }else if(type == TYPE_YUV_444){
            U = new2d<Pixel>(h,w,0);
            V = new2d<Pixel>(h,w,0);
        }else{
            U = V = NULL;
        }
    }

    ~Frame(){
        if(Y!=NULL)
            delete2d<Pixel>(Y);
        if(img_bp!=NULL)
            delete3d<int8>(img_bp);
        if(U!=NULL)
            delete2d<Pixel>(U);
        if(V!=NULL)
            delete2d<Pixel>(V);
        if(x!=NULL)
            DELETE(x);
    }

    void read(FILE* fptr, const int skip = 0 ){
        if(fptr==NULL){
            printf("FILE can't be read");
            return;
        }

//        unsigned char*  buffer = MALLOC(unsigned char,lm*3/2);//(unsigned char*)malloc(sizeof(unsigned char)*lm*3/2);
//        const int offset_u = lm, offset_v = lm*5/4;
        fseek(fptr,lm*3/2*skip,SEEK_CUR);
        fread(Y[0],1,lm,fptr);

//        fread(buffer,1,lm*3/2,fptr);
//
//        for(int i = 0 ; i < lm ; ++i)
//            Y[0][i] = (int) buffer[i];

        if(U!=NULL && V!=NULL && type == TYPE_YUV_420){
//            for(int i = 0 ; i < lm/4 ; ++i){
//                U[0][i] = (int) buffer[offset_u + i];
//                V[0][i] = (int) buffer[offset_v + i];
//            }
            fread(U[0],1,lm/4,fptr);
            fread(V[0],1,lm/4,fptr);
        }else
            fseek(fptr,lm/2,SEEK_CUR);

//        DELETE(buffer);
        img2bp_frame(Y,height,width,img_bp);
    }

    void write(FILE* fptr){
        if(fptr==NULL){
            printf("FILE can't be read");
            return;
        }
        // TODO:
        unsigned char* buffer = (unsigned char*)malloc(sizeof(unsigned char)*lm);

        for(int i = 0, j = 0 ; i < height ; ++i){
            for(j = 0 ; j < width ; ++j)
                buffer[j+i*width] = Y[i][j];
        }
        fwrite(buffer,1,lm,fptr);

        DELETE(buffer);
    }


    void encode_info(const double &EbN0dB, int** G_ptr){
        puncture = 0;                     // puncture or not
        rate = 1/(double)(2-puncture);       // code rate
        a = 1;                         // Fading amplitude. a=1 -> AWGN channel
        EbN0 = pow(10,EbN0dB/10);      // convert Eb/N0[dB] to normal number
        G = G_ptr;
    }

    void encode(double** Ly, int** map_out, double* weights){

        double L_c[PXL];                            // reliability value of the channel
        double sigma[PXL];                          // standard deviation of AWGN noise

        for(int i = 0 ; i < PXL ; ++i){
            L_c[i] = 4*a*weights[i]*EbN0*rate;
            sigma[i] = a*sqrt(2*rate*weights[i]*EbN0);
        }

        for(int t_lvl = 0 ; t_lvl < PXL ; ++t_lvl){
            // interleave
            random_sequence(0,lm-1,map_out[t_lvl]);
            rsc_encode(G,G_L,img_bp[t_lvl],map_out[t_lvl],lm,width,1,x);

            for(int i = 0 ; i < 2*lu ; ++i)
                Ly[t_lvl][i] = 0.5*L_c[t_lvl]*(2*x[i] - 1) + sigma[t_lvl]*gaussian_noise();  // add noise   + sigma*gaussian_noise()
        }
    }

    void next(FILE* fptr, double** Ly, int** map_out, double* weights ){
        read(fptr);
        encode(Ly,map_out,weights);
    }

    void copy(Frame* frame){
        for(int i = 0 ; i < lm ; ++i)
            Y[0][i] = frame->Y[0][i];
        for(int i = 0 ; i < PXL*lm ; ++i)
            img_bp[0][0][i] = frame->img_bp[0][0][i];

        if(type==TYPE_YUV_420){
            for(int i = 0 ; i < lm/4 ; ++i){
                U[0][i] = frame->U[0][i];
                V[0][i] = frame->V[0][i];
            }
        }
    }

    double psnr(Pixel** img){
        return computePSNR(img,Y,lm);
    }

    int getWidth(){return width;}
    int getHeight(){return height;}
    int getLm(){return lm;}
    int getTag(){return frame_tag;}

    Pixel** Y ;
    Pixel** U ;
    Pixel** V ;
    int8*** img_bp ;
/*
    void getImgBp(int*** bp){
        if(img_bp!=NULL){
            bp = img_bp;
            return;
        }

        img_bp = new3d<int>(PXL,height,width,0);
        img2bp_frame(Y,height,width,img_bp);
        bp = img_bp;
    }
*/

private:

    int width;
    int height;
    int lm;
    int type;
    int frame_tag;
    int lu;

    // for encode

    int8* x ;
    int** G ;
    double EbN0 ;
    double a;
    int puncture ;
    double rate ;       // code rate

};


#endif // __FRAME_H

