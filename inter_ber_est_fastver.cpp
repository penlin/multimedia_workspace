/*
* 72% timeslape usage than original version , but only record the PSNR information,  BER information would not print
*/
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
    double EbN0[len] ;

    FILE* fptr;
    if(argc > 2)
        fptr = fopen(FILENAME[atoi(argv[2])],"r+b");
    else
        fptr = fopen(__SEQ__,"r+b");
    assert(fptr!=NULL);
    rewind(fptr);

    Frame* frame = new Frame(h,w,TYPE_Y,0);
    Frame* frame_prev = new Frame(h,w,TYPE_Y,1);
    Frame* frame_next = new Frame(h,w,TYPE_Y,2);

    double PSNR[len][f] ;//= MALLOC(double,f);
    double* weights = MALLOC(double,PXL);
    int** mv = new2d<int>(2,mv_len,0);
    int** mv_prev = new2d<int>(2,mv_len,0);

    for(int i = 0 ; i < len ; ++i){
        snr[i] = __SNR_S+i;
        EbN0[i] = pow(10,snr[i]/10);
    }

    int weight_type = 0;
    // 0: EEP
    // 1: UEP
    if(argc > 1)
        weight_type = atoi(argv[1]);
    if(argc > 3)
        f = atoi(argv[3]);

    fseek(fptr,h*w*3/2*__SKIP,SEEK_SET);
    for(int j = 0 ; j < PXL ; ++j)
        weights[j] = 1;//(weight_type?DERIVE_INTRA_SNR[i][j]:1);
    frame->read(fptr);
    frame_next->read(fptr);
    motionEstES(frame->Y,frame_next->Y,h,w,8,5,mv);

    for(int i = 0 ; i < len ; ++i)
        PSNR[i][0] = inter_psnr_est(frame->img_bp,frame_next->img_bp,mv,h,w,weights,EbN0[i]);

    for(int j = 1 ; j < f - 1; ++ j){
        frame_prev->copy(frame);
        frame->copy(frame_next);
        frame_next->read(fptr);

        motionEstES(frame->Y,frame_prev->Y,h,w,8,5,mv_prev);
        motionEstES(frame->Y,frame_next->Y,h,w,8,5,mv);

        for(int i = 0 ; i < len ; ++i)
            PSNR[i][j] = inter_psnr_est(frame->img_bp,frame_prev->img_bp,mv_prev,h,w,weights,EbN0[i],frame_next->img_bp,mv);
//              PSNR[j] = inter_psnr_est(frame_prev->img_bp,frame->img_bp,mv_prev,h,w,weights,EbN0);

    }

    motionEstES(frame_next->Y,frame->Y,h,w,8,5,mv_prev);

    for(int i = 0 ; i < len ; ++i)
        PSNR[i][f-1] = inter_psnr_est(frame_next->img_bp,frame->img_bp,mv_prev,h,w,weights,EbN0[i]);

    // print result PSNR
    double avg_psnr[len] ;
    for(int i = 0 ; i < len ; ++i)
        avg_psnr[i] = 0.0;
    for(int i = 0 ; i < f ; ++i){
        for(int j = 0 ; j < len ; ++j){
            printf("%lf, ",PSNR[j][i]);
            avg_psnr[j] += PSNR[j][i];
        }
        printf("\n");
    }
    printf("AVERAGE PSNR = ");
    for(int i = 0 ; i < len ; ++i)
        printf("%lf, ",avg_psnr[i]/f);
    printf("\n");


    // free memory
    free(weights);
    delete2d<int>(mv);
    delete2d<int>(mv_prev);

    fclose(fptr);
    delete frame;
    delete frame_prev;
    delete frame_next;
}
