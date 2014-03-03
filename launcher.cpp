#include "build_value.h"
#include "io_utils.h"
#include "video_encode.h"
#include "data_alloc.h"
#include "algs_direct.h"
#include "algs_intra.h"
#include "algs_inter_weiting.h"
#include "algs_inter_intra.h"
#include "utils.h"

#define DECODE __ALGO__


void de_dpcm(int*** Y, int*** dpcm_Y, const int& h, const int&w, const int &f){

    int i,j,k;

    for(k=1;k<f;++k){

        for(i=0;i<h;++i)
            for(j=0;j<w;++j){
                dpcm_Y[k][i][j] = dpcm_Y[k][i][j];
                dpcm_Y[k][i][j] += dpcm_Y[k-1][i][j];
                dpcm_Y[k][i][j]%=256;
            }

        printf("PSNR[%d]=%f\n",k+1,computePSNR(Y[k],dpcm_Y[k],h,w));
    }


}


int main(int argc,char* argv[]){

    startRandom();

    const int h = __HEIGHT, w = __WIDTH, f = __FRAME ;
    const int lm = h*w;
    const int lu = lm+(G_L-1);
    int snr = __SNR;
    int ** G = getGenerator();
    double* PSNR = (double*) malloc(sizeof(double)*f);
    double *** Ly = new3d<double>(f,PXL,2*lu);
    int *** map_out = new3d<int>(f,PXL,lm);
    int*** Y = new3d<int>(f,h,w);
    int*** dpcm_Y = new3d<int>(f,h,w);
    int*** _Y = new3d<int>(f,h,w);
    int*** U = new3d<int>(f,h/2,w/2);
    int*** V = new3d<int>(f,h/2,w/2);

    // read YUV
    if(argc == 1)
//        Y = yuv_read("stefan_cif.yuv",h,w,f);
        yuv_read("hall_cif.yuv",h,w,f,Y,U,V);
    else
//        Y = yuv_read(argv[1],h,w,f);
        yuv_read(argv[1],h,w,f,Y,NULL,NULL);

//    dpcm(Y,dpcm_Y,h,w,f);

//    de_dpcm(Y,dpcm_Y,h,w,f);
    // encode
    video_encode(Y,f,h,w,snr,G,Ly,map_out);

    // decode
    if(argc == 1)
        DECODE("hall",Y,Ly,map_out,h,w,f,lu,G,PSNR,_Y);
    else if(argc > 2)
        DECODE(argv[2],Y,Ly,map_out,h,w,f,lu,G,PSNR,NULL);
    else
        DECODE("NONAME",Y,Ly,map_out,h,w,f,lu,G,PSNR,NULL);

    // compute Ia
//    double** Ia = new2d<double>(f,PXL);
//    computeIa(Ly,map_out,Y,f,h,w,Ia);


    // print result PSNR
    for(int i = 0 ; i < f ; ++i)
        printf("frame#%3d PSNR = %lf\n",i+1,PSNR[i]);

#if __OUTPUT_SEQ__
    yuv_write("hall",_Y,U,V,f,h,w);
#endif
    write_psnr_info(PSNR);

    // free memory
    delete3d<double>(Ly);
    delete3d<int>(map_out);
    delete2d<int>(G);
    deleteY(Y);
    deleteY(_Y);
    deleteY(dpcm_Y);
    deleteY(U);
    deleteY(V);

    free(PSNR);
}
