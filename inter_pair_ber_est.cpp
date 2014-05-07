#include "frame.h"
#include "uep_predict_utils.h"

const char* FILENAME[4] = {__FOREMAN, __HALL, __STEFAN, __AKIYO};

int main(int argc,char* argv[]){

    startRandom();

    const int h = __HEIGHT, w = __WIDTH, mbSize = 8;
    int f = __FRAME ;
    const int mv_len = h*w/mbSize/mbSize;
    const int len = __SNR_E-__SNR_S+1;
    double snr[len];
    double EbN0 ;
#ifdef MEM_MGR_H
    initMemMgr(h*w*PXL*sizeof(double)*3);
#endif
    FILE* fptr;
    if(argc > 2)
        fptr = fopen(FILENAME[atoi(argv[2])],"r+b");
    else
        fptr = fopen(__SEQ__,"r+b");
    assert(fptr!=NULL);
    rewind(fptr);

    Frame* frame = new Frame(h,w,TYPE_Y,0);
    Frame* frame_prev = new Frame(h,w,TYPE_Y,1);

    double PSNR[f] ;//= MALLOC(double,f);
    double* weights = MALLOC(double,PXL);
    int** mv = new2d<int>(2,mv_len,0);
    int** mv_prev = new2d<int>(2,mv_len,0);

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
        EbN0 = pow(10,snr[i]/10);
        for(int j = 0 ; j < PXL ; ++j){
            weights[j] = 1;//(weight_type?DERIVE_INTRA_SNR[i][j]:1);
        }

        for(int j = 0 ; j < f - 1; j+=2){
            frame_prev->read(fptr);
            frame->read(fptr);

//            motionEstES(frame->Y,frame_prev->Y,h,w,8,5,mv_prev);
//            motionEstES(frame_prev->Y,frame->Y,h,w,8,5,mv);

            PSNR[j] = inter_basic_psnr_est(frame->img_bp,frame_prev->img_bp,mv_prev,h,w,weights,EbN0);
            PSNR[j+1] = inter_basic_psnr_est(frame_prev->img_bp,frame->img_bp,mv,h,w,weights,EbN0);

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

    }

    // free memory

    fclose(fptr);
    delete frame;
    delete frame_prev;
#ifdef MEM_MGR_H
    freeMemMgr();
#else
    DELETE(weights);
    DELETE(frr_prev);
    delete2d<int>(mv);
    delete2d<int>(mv_prev);
#endif
}
