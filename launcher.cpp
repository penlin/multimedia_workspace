#include "algs_direct.h"
#include "algs_intra.h"
#include "algs_inter.h"
#include "io_utils.h"
#include "uep_predict_utils.h"

#define DECODE __ALGO__


const double DERIVE_SNR[5][PXL]={{2.597844, 2.244942, 1.760056, 1.214787, 0.130591, 0.017250, 0.017250, 0.017250},
                                {2.382862, 1.963794, 1.582247, 1.234872, 0.765783, 0.043096, 0.013702, 0.013702},
                                {2.092984, 1.818853, 1.521252, 1.218653, 0.921529, 0.402081, 0.013660, 0.010884},
                                {1.933909, 1.607395, 1.418209, 1.193889, 0.953517, 0.706520, 0.176923, 0.009999},
                                {1.720262, 1.532340, 1.260163, 1.118219, 0.943704, 0.752738, 0.552381, 0.121300}};

const double MINPOW_PSNR25[PXL] = {2.357380, 1.877585, 1.364773, 0.268487, 0.019715, 0.013750, 0.013750, 0.013750};


int main(int argc,char* argv[]){

    startRandom();

    const int h = __HEIGHT, w = __WIDTH, f = __FRAME ;

    const int len = __SNR_E-__SNR_S+1;
    double snr[len];

    // pstate, pout
    int ** G = getGenerator();
    const int Ns = pow(2,G_L-1);
    int ** pout = new2d<int>(Ns,4);
    int ** pstate = new2d<int>(Ns,2) ;

    FILE* fptr = fopen(__SEQ__,"r+b");
    assert(fptr!=NULL);
    rewind(fptr);

    double* PSNR = (double*) malloc(sizeof(double)*f);
    double* weights = (double*)malloc(sizeof(double)*PXL);

    trellis(G,G_N,G_L,Ns,pout,pstate);

    for(int i = 0 ; i < len ; ++i)
        snr[i] = __SNR_S+i;

//    weight_predict_minMSE(weights,pow(10,snr/10));


    for(int i = 0 ; i < len; ++i){
        fseek(fptr,h*w*3/2*__SKIP,SEEK_SET);

        for(int j = 0 ; j < PXL ; ++j)
            weights[j] = DERIVE_SNR[i][j];

        // decode
//        direct_system(__TAG__,fptr,h,w,f,G,pout,pstate,snr[i],weights,PSNR,NULL);
//        intra_system(__TAG__,fptr,h,w,f,G,pout,pstate,snr[i],weights,PSNR,NULL);
//        inter_system(__TAG__,fptr,h,w,f,G,pout,pstate,snr[i],weights,PSNR,NULL);
        if(argc == 1)
            DECODE(__TAG__,fptr,h,w,f,G,pout,pstate,snr[i],weights,PSNR,NULL);
        else if(argc > 2)
            DECODE(argv[2],fptr,h,w,f,G,pout,pstate,snr[i],weights,PSNR,NULL);
        else
            DECODE(argv[1],fptr,h,w,f,G,pout,pstate,snr[i],weights,PSNR,NULL);

        // print result PSNR
#if __PROGRESS__
        double avg_psnr = 0.0;
        for(int j = 0 ; j < f ; ++j){
            printf("frame#%3d PSNR = %lf\n",j+1,PSNR[j]);
            avg_psnr += PSNR[j];
        }
        printf("AVERAGE PSNR = %lf\n",avg_psnr/f);
#endif

#if __OUTPUT_SEQ__
        //yuv_write(__TAG__,_Y,U,V,f,h,w);
#endif
        write_psnr_info(PSNR,snr[i]);
    }

    // free memory
    delete2d<int>(G);
    delete2d<int>(pstate);
    delete2d<int>(pout);
    free(PSNR);
    free(weights);
    fclose(fptr);
}
