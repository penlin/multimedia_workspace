/*
* 72% timeslape usage than original version , but only record the PSNR information,  BER information would not print
*/
#include "frame.h"
#include "uep_predict_utils.h"

const char* FILENAME[4] = {__FOREMAN, __HALL, __STEFAN, __AKIYO};

void calInterEn(int*** img_bp, int*** img_bp_ref, int** mv1, const int &imgh, const int &imgw, double* eEn, int*** img_bp_ref2= NULL, int** mv2 = NULL ){

    int** img_ref = new2d<int>(imgh,imgw);
    int** img_ref2;
    double En = 0.0, En1 = 0.0, beta, beta_prev;
    if(img_bp_ref2!=NULL)
        img_ref2 = new2d<int>(imgh,imgw);

    for(int i = 0 , j = 0, k = 0; i < PXL ; ++i){

        motionComp(img_bp_ref[i],mv1,imgh,imgw,8,img_ref);
        En = En1 = 0.0;
        if(img_bp_ref2==NULL){
            inter_beta_estimation(img_bp[i],beta,img_ref,imgh,imgw);
            for(j= 0 ; j < imgh ; ++j)
                for(k=0 ; k < imgw ; ++k)
                    En+= (2*img_bp[i][j][k]-1)*(2*img_ref[j][k]-1);

            eEn[i] = En*beta;///(imgh*imgw);
        } else {
            motionComp(img_bp_ref2[i],mv2,imgh,imgw,8,img_ref2);
            inter2_beta_estimation(img_bp[i], img_ref2, img_ref, beta, beta_prev, imgh, imgw);
            for(j= 0 ; j < imgh ; ++j)
                for(k=0 ; k < imgw ; ++k){
                    En+= (2*img_bp[i][j][k]-1)*(2*img_ref[j][k]-1);
                    En1+= (2*img_bp[i][j][k]-1)*(2*img_ref2[j][k]-1);
                }

            eEn[i] = (En*beta_prev+En1*beta);
        }

        eEn[i] /= (imgh*imgw);
    }
    delete2d<int>(img_ref);
    if(img_bp_ref2!=NULL){
        delete2d<int>(img_ref2);
    }

}

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
    double* eEn = MALLOC(double,PXL);
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
//    motionEstES(frame->Y,frame_next->Y,h,w,8,5,mv);
    calInterEn(frame->img_bp,frame_next->img_bp,mv,h,w,eEn);
    for(int i = 0 ; i < len ; ++i)
        PSNR[i][0] = mrf_psnr_est(weights,EbN0[i],eEn);

    for(int j = 1 ; j < f - 1; ++ j){
        frame_prev->copy(frame);
        frame->copy(frame_next);
        frame_next->read(fptr);

//        motionEstES(frame->Y,frame_prev->Y,h,w,8,5,mv_prev);
//        motionEstES(frame->Y,frame_next->Y,h,w,8,5,mv);

        calInterEn(frame->img_bp,frame_prev->img_bp,mv_prev,h,w,eEn,frame_next->img_bp,mv);
        for(int i = 0 ; i < len ; ++i)
            PSNR[i][j] = mrf_psnr_est(weights,EbN0[i],eEn);
    }

//    motionEstES(frame_next->Y,frame->Y,h,w,8,5,mv_prev);
    calInterEn(frame_next->img_bp,frame->img_bp,mv_prev,h,w,eEn);
    for(int i = 0 ; i < len ; ++i)
        PSNR[i][f-1] = mrf_psnr_est(weights,EbN0[i],eEn);

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
    free(eEn);
    delete2d<int>(mv);
    delete2d<int>(mv_prev);

    fclose(fptr);
    delete frame;
    delete frame_prev;
    delete frame_next;
}
