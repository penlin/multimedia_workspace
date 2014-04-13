#include "build_value.h"
#include "io_utils.h"
#include "video_encode.h"
#include "data_alloc.h"
#include "utils.h"
#include "me_compare.h"

int main(int argc,char* argv[]){

    startRandom();

    const int h = __HEIGHT, w = __WIDTH, f = 10 ;
    const int lm = h*w;
    const int lu = lm+(G_L-1);
    int snr = __SNR;
    int ** G = getGenerator();
//    double* PSNR = (double*) malloc(sizeof(double)*f);
    double** PSNR = new2d<double>(3,f,0);
    double *** Ly_in = new3d<double>(f,PXL,2*lu);
    int *** map_out = new3d<int>(f,PXL,lm);
    int*** Y = new3d<int>(f,h,w);

    // read YUV
    if(argc == 1)
//        Y = yuv_read("stefan_cif.yuv",h,w,f);
        yuv_read("stefan_cif.yuv",h,w,f,Y,NULL,NULL);
    else
//        Y = yuv_read(argv[1],h,w,f);
        yuv_read(argv[1],h,w,f,Y,NULL,NULL);

    // encode
    video_encode(Y,f,h,w,snr,G,Ly_in,map_out);

    // decode
    compare("stefan",Y,Ly_in,map_out,h,w,f,lu,G,PSNR);


    for(int i = 0 , j = 0 ; i < 3 ; ++i){
        for(j = 1 ; j < f ; ++j)
            printf("%3.3f ",PSNR[i][j]);
        printf("\n");
    }

    // free memory
    delete3d<double>(Ly_in);
    delete3d<int>(map_out);
    delete2d<int>(G);
    deleteY(Y);
    delete2d<double>(PSNR);

}
