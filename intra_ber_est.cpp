#include "frame.h"
#include "uep_predict_utils.h"
#include "mrf_decoder_utils.h"


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

    FILE* fptr = fopen(__SEQ__,"r+b");
    assert(fptr!=NULL);
    rewind(fptr);
    Frame* frame = new Frame(h,w,TYPE_Y,0);
    double* PSNR = MALLOC(double,f);
    double* weights = MALLOC(double,PXL);
//    double* beta = MALLOC(double,PXL);

    for(int i = 0 ; i < len ; ++i)
        snr[i] = __SNR_S+i;

    int weight_type = 0;
    // 0: EEP
    // 1: UEP
    if(argc > 1)
        weight_type = atoi(argv[1]);


    for(int i = 0 ; i < len; ++i){
        fseek(fptr,h*w*3/2*__SKIP,SEEK_SET);

        for(int j = 0 ; j < PXL ; ++j)
            weights[j] = (weight_type?DERIVE_SNR[i][j]:1);
        for(int j = 0 ; j < f ; ++ j){
            frame->read(fptr);
            //intra_beta_estimation(frame->img_bp,beta,h,w);
            PSNR[j] = intra_psnr_est(frame->img_bp,h,w,weights,pow(10,snr[i]/10));
        }

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
    }

    // free memory
    free(PSNR);
    free(weights);
    fclose(fptr);
}
