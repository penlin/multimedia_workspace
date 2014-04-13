#include <stdio.h>
#include <stdlib.h>
#include "build_value.h"
#include "io_utils.h"
#include "data_alloc.h"
#include "utils.h"
#include "image_process_utils.h"
#include "exit_chart_utils.h"
#include "motion_estimation_utils.h"
#include "mrf_decoder_utils.h"

int main(int argc,char* argv[]){

    startRandom();

    // control SNR step
    const double s_snr = -1;
    const double step = 0.05;
    const int snr_size = 500;

    // control the certain Frame and the bp
    int frame = 2;

    if(argc >= 4)
        frame = strtod(argv[3],NULL);

    // for video encode
    const int puncture = 0;                     // puncture or not
    const double rate = 1/(double)(2-puncture);       // code rate
    const double a = 1;                         // Fading amplitude. a=1 -> AWGN channel
    double EbN0,L_c,sigma;

    // for basical info
    char buffer[50];
    const int h = __HEIGHT, w = __WIDTH, f = frame+2 ;
    const int lm = h*w;
    const int lu = lm+(G_L-1);
    int*** Y = new3d<int>(f,h,w);

    int** MV_prev = new2d<int>(2,lm/64,0);
    int** MV = new2d<int>(2,lm/64,0);

    double** Lu = new2d<double>(PXL,lu,0);
    double** Lu_prev = new2d<double>(PXL,lu,0);
    double** Lu_next = new2d<double>(PXL,lu,0);

    // frame buffer
    int*** imgr_bp = new3d<int>(PXL,h,w);
    int*** imgr_bp_prev = new3d<int>(PXL,h,w);
    int*** imgr_bp_next = new3d<int>(PXL,h,w);
    int*** img = new3d<int>(PXL,h,w);
    int*** img_prev = new3d<int>(PXL,h,w);
    int*** img_next = new3d<int>(PXL,h,w);
    int** imgr = new2d<int>(h,w);
    int** imgr_prev = new2d<int>(h,w);
    int** imgr_next = new2d<int>(h,w);

    double** Lu_c = new2d<double>(PXL,lu);  //channel decoder output

    // buffer for Ia, Ie
    double** Le = new2d<double>(PXL,lm); // source extrinsic informatiom
    double** Le_s_inter = new2d<double>(PXL,lm,0); // source extrinsic information
    double** Le_s_inter_next = new2d<double>(PXL,lm,0); // source extrinsic information
    double** Le_s_intra = new2d<double>(PXL,lm,0); // source extrinsic information
    double** Ia = new2d<double>(PXL,snr_size);
    double** Ie = new2d<double>(PXL,snr_size);
    int idx = 0;
    double tmp_sum = 0, tmp ;

    double* beta_t = (double*) malloc(sizeof(double)*PXL);
    double* beta_s = (double*) malloc(sizeof(double)*PXL);
    double* beta_t_prev = (double*) malloc(sizeof(double)*PXL);

    // read YUV
    sprintf(buffer,"%s_cif.yuv",argv[1]);
    yuv_read(buffer,h,w,f,Y,NULL,NULL);


    img2bp_frame(Y[frame],h,w,imgr_bp);
    img2bp_frame(Y[frame-1],h,w,imgr_bp_prev);
    img2bp_frame(Y[frame+1],h,w,imgr_bp_next);

    printf("processing ...  %2.2f%%",0);
    for(double snr = s_snr; idx < snr_size ; snr+=step, idx++ ){
        printf("\rprocessing ... %2.2f%%",100*(snr-s_snr)/step/snr_size);
        // encode
        EbN0 = pow(10,snr/10);      // convert Eb/N0[dB] to normal number
        L_c = 4*a*EbN0*rate;           // reliability value of the channel
        sigma = 1/sqrt(2*rate*EbN0);   // standard deviation of AWGN noise

        for(int t_lvl = 0 ; t_lvl <PXL ; ++t_lvl){
            tmp_sum = 0;
            for(int i = 0 ; i < lm ; ++i){
                tmp = (2*imgr_bp[t_lvl][i/w][i%w] - 1);
                Lu[t_lvl][i] = 0.5*L_c*( tmp + sigma*gaussian_noise());
                tmp_sum+=log2(1 + exp(-Lu[t_lvl][i]*tmp));

                tmp = (2*imgr_bp_prev[t_lvl][i/w][i%w] - 1);
                Lu_prev[t_lvl][i] = 0.5*L_c*( tmp + sigma*gaussian_noise());
            }

            Ia[t_lvl][idx] = 1 - tmp_sum/lm;
        }

        // decode

        for(int i = 0, j=0, t_lvl=0; i <h ; ++i)
            for(j = 0 ; j < w ; ++j)
                for(t_lvl=0 ; t_lvl < PXL ; ++t_lvl){
                    img[t_lvl][i][j] = ((Lu[t_lvl][j+i*w]>=0)?1:0);
                    img_prev[t_lvl][i][j] = ((Lu_prev[t_lvl][j+i*w]>=0)?1:0);
                }

        bin2dec_img(img,h,w,imgr);
        bin2dec_img(img_prev,h,w,imgr_prev);

        motionEstES(imgr,imgr_prev,h,w,8,5,MV_prev);

        intra_inter_beta_estimation(img,img_prev,MV_prev,beta_s,beta_t,h,w,8);

        motionComp(Lu_prev,MV_prev,h,w,8,Lu_c);
        mrf_siso_inter(Lu,Lu_c,beta_t,h,w,Le_s_inter,0);
        mrf_siso_intra(Lu,beta_s,h,w,Le_s_intra,0);

        for(int i = 0, t_lvl = 0 ; i < lm ; ++i)
            for(t_lvl = 0 ; t_lvl < PXL ; ++t_lvl)
                Le[t_lvl][i] = (Le_s_inter[t_lvl][i] + Le_s_intra[t_lvl][i]);  // Le_s_prev in the next frame

        // prepare for next frame decode

        // encode

        for(int t_lvl = 0 ; t_lvl <PXL ; ++t_lvl){
            for(int i = 0 ; i < lm ; ++i){
                tmp = (2*imgr_bp_next[t_lvl][i/w][i%w] - 1);
                Lu_next[t_lvl][i] = 0.5*L_c*( tmp + sigma*gaussian_noise());  // reuse as Lu_next
            }
        }

        for(int i = 0, j=0, t_lvl=0; i <h ; ++i)
            for(j = 0 ; j < w ; ++j)
                for(t_lvl=0 ; t_lvl < PXL ; ++t_lvl){
                    img_next[t_lvl][i][j] = ((Lu_next[t_lvl][j+i*w]>=0)?1:0);
                }

        bin2dec_img(img_next,h,w,imgr_next);

        motionEstES(imgr,imgr_next,h,w,8,5,MV);

        if(frame==1){
            intra_inter_beta_estimation(imgr_bp,imgr_bp_next,MV,beta_s,beta_t,h,w,8);
            for(int i = 0 ; i < PXL ; ++i)
                beta_t_prev[i] = 0;
        }
        else
            intra_inter2_beta_estimation(imgr_bp_next,imgr_bp,imgr_bp_prev,MV,beta_s,beta_t,beta_t_prev,h,w,8);

        motionComp(Lu_next,MV,h,w,8,Lu_c);
        mrf_siso_inter(Lu,Lu_c,beta_t,h,w,Le_s_inter_next,0);
        mrf_siso_intra(Lu,beta_s,h,w,Le_s_intra,0);

        for(int i = 0, t_lvl = 0 ; i < lm ; ++i)
            for(t_lvl = 0 ; t_lvl < PXL ; ++t_lvl)
                Le[t_lvl][i] = (Le_s_inter_next[t_lvl][i] + Le_s_intra[t_lvl][i] + Le_s_inter[t_lvl][i]*beta_t_prev[t_lvl]);


        // compute Ie
        for(int t_lvl = 0 , i = 0; t_lvl < PXL ; ++t_lvl){
            tmp_sum = 0;
            for(i=0; i < lm ; ++i){
                tmp_sum+=log2(1 + exp(-Le[t_lvl][i]*(2*imgr_bp[0][0][i+t_lvl*lm]-1)));

            }
            Ie[t_lvl][idx] = 1 - tmp_sum/lm;
        }


    }

    printf("\r100.00%% completed!\n");

    FILE *file = fopen("output/exit_curve_ps.txt","a+");

    for(int t_lvl = 0 ; t_lvl < PXL ; ++t_lvl){
        fprintf(file,"==== v2 ========================\n%s: frame#%d,bp#%d\nIa=\n",argv[1],frame+1,t_lvl+1);
        for(int i = 0 ; i < snr_size ; ++i)
            fprintf(file,"%lf,",Ia[t_lvl][i]);
        fprintf(file,"\nIe=\n");
        for(int i = 0 ; i < snr_size ; ++i)
            fprintf(file,"%lf,",Ie[t_lvl][i]);
        fprintf(file,"\n\n");

    }
    fclose(file);

    write_ps_for_matlab(Ia,Ie,snr_size);

    // free memory
    deleteY(Y);

    delete3d<int>(imgr_bp);
    delete3d<int>(imgr_bp_prev);
    delete3d<int>(imgr_bp_next);
    delete3d<int>(img_next);
    delete3d<int>(img_prev);
    delete3d<int>(img);
    delete2d<int>(imgr_next);
    delete2d<int>(imgr_prev);
    delete2d<int>(imgr);

    delete2d<double>(Lu_c);

    delete2d<double>(Lu_next);
    delete2d<double>(Lu);
    delete2d<double>(Lu_prev);

    delete2d<double>(Ia);
    delete2d<double>(Ie);

    delete2d<int>(MV_prev);
    delete2d<int>(MV);

    delete2d<double>(Le);
    delete2d<double>(Le_s_inter);
    delete2d<double>(Le_s_inter_next);
    delete2d<double>(Le_s_intra);

    free(beta_s);
    free(beta_t);
    free(beta_t_prev);
}
