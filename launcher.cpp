#include "build_value.h"
#include "algs_direct.h"
#include "algs_intra.h"
#include "algs_inter.h"
#include "algs_inter_pair.h"
#include "io_utils.h"


#define DECODE __ALGO__


const double DERIVE_SNR[5][PXL]={{2.597844, 2.244942, 1.760056, 1.214787, 0.130591, 0.017250, 0.017250, 0.017250},
                                {2.382862, 1.963794, 1.582247, 1.234872, 0.765783, 0.043096, 0.013702, 0.013702},
                                {2.092984, 1.818853, 1.521252, 1.218653, 0.921529, 0.402081, 0.013660, 0.010884},
                                {1.933909, 1.607395, 1.418209, 1.193889, 0.953517, 0.706520, 0.176923, 0.009999},
                                {1.720262, 1.532340, 1.260163, 1.118219, 0.943704, 0.752738, 0.552381, 0.121300}};

const double MINPOW_PSNR25[PXL] = {2.357380, 1.877585, 1.364773, 0.268487, 0.019715, 0.013750, 0.013750, 0.013750};

const double DERIVE_INTRA_SNR[5][PXL] = {{1.542550, 1.738844, 1.455300, 1.401314, 1.168949, 0.641986, 0.041631, 0.010000},
                                        {1.468446, 1.559849, 1.369071, 1.310368, 1.211347, 0.925365, 0.143461, 0.011934},
                                        {1.316359, 1.459535, 1.241203, 1.226068, 1.163622, 0.972940, 0.578915, 0.041608},
                                        {1.225215, 1.260888, 1.192017, 1.174171, 1.100385, 0.963340, 0.748424, 0.335652},
                                        {1.178798, 1.228975, 1.131290, 1.106988, 1.014853, 0.959307, 0.779196, 0.600574}};

const char* FILENAME[4] = {__FOREMAN, __HALL, __STEFAN, __AKIYO};

int main(int argc,char* argv[]){

    startRandom();

    const int h = __HEIGHT, w = __WIDTH;
    int f = __FRAME ;
#ifdef MEM_MGR_H
    initMemMgr(h*w*PXL*sizeof(double)*5);
#endif

    const int len = __SNR_E-__SNR_S+1;
    double snr[len];

    // pstate, pout
    int ** G = getGenerator();
    const int Ns = pow(2,G_L-1);
    int ** pout = new2d<int>(Ns,4);
    int ** pstate = new2d<int>(Ns,2) ;

    FILE* fptr;
    if(argc > 2)
        fptr = fopen(FILENAME[atoi(argv[2])],"r+b");
    else
        fptr = fopen(__SEQ__,"r+b");

    assert(fptr!=NULL);
    rewind(fptr);

    double* PSNR = MALLOC(double,f) ;//(double*) malloc(sizeof(double)*f);

    trellis(G,G_N,G_L,Ns,pout,pstate);

    for(int i = 0 ; i < len ; ++i)
        snr[i] = __SNR_S+i;

    int weight_type = 0;
    // 0: EEP
    // 1: UEP
    if(argc > 1)
        weight_type = atoi(argv[1]);
    if(argc > 3)
        f = atoi(argv[3]);

    for(int i = 0 ; i < len; ++i){
        fseek(fptr,h*w*3/2*__SKIP,SEEK_SET);

//        for(int j = 0 ; j < PXL ; ++j)
//            weights[j] = 1;//DERIVE_INTRA_SNR[i][j];//(weight_type?DERIVE_SNR[i][j]:1);

        // decode
//        direct_system(__TAG__,fptr,h,w,f,G,pout,pstate,snr[i],weights,PSNR,NULL);
//        intra_system(__TAG__,fptr,h,w,f,G,pout,pstate,snr[i],weights,PSNR,NULL);
//        inter_system(__TAG__,fptr,h,w,f,G,pout,pstate,snr[i],weights,PSNR,NULL);
        if(argc > 2)
            DECODE(FILENAME[atoi(argv[2])],fptr,h,w,f,G,pout,pstate,snr[i],PSNR,weight_type,NULL);
        else
            DECODE(__TAG__,fptr,h,w,f,G,pout,pstate,snr[i],PSNR,weight_type,NULL);

        // print result PSNR
#if __PSNR__
        double avg_psnr = 0.0;
        for(int j = 0 ; j < f ; ++j){
            printf("frame#%3d PSNR = %lf\n",j+1,PSNR[j]);
            avg_psnr += PSNR[j];
        }
        printf("AVERAGE PSNR = %lf\n",avg_psnr/f);
#else
        printf("===========END============\n");
#endif

#if __OUTPUT_SEQ__
        //yuv_write(__TAG__,_Y,U,V,f,h,w);
#endif

        write_psnr_info(PSNR,snr[i]);
#ifdef MEM_MGR_H
        printMem();
#endif
    }

    // free memory
    delete2d<int>(G);
    delete2d<int>(pstate);
    delete2d<int>(pout);
    DELETE(PSNR);//free(PSNR);
    fclose(fptr);
#ifdef MEM_MGR_H
    freeMemMgr();
#endif

}
