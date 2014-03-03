#include <stdio.h>
#include <stdlib.h>
#include "build_value.h"
#include "io_utils.h"
#include "video_encode.h"
#include "data_alloc.h"
#include "utils.h"
#include "image_process_utils.h"


int main(int argc,char* argv[]){

    startRandom();

    // control SNR step
    const double s_snr = -0.5;
    const double step = 0.02;
    const int snr_size = 500;

    // control the certain Frame and the bp
    int frame = 2;
    int t_lvl = 0;
    double _snr = 0;
    if(argc >= 3 )
        _snr = strtod(argv[2],NULL);

    if(argc >= 4)
        frame = strtod(argv[3],NULL);

    if(argc >= 5)
        t_lvl = strtod(argv[4],NULL);

    // for video encode
    const int puncture = 0;                     // puncture or not
    const double rate = 1/(double)(2-puncture);       // code rate
    const double a = 1;                         // Fading amplitude. a=1 -> AWGN channel
    double EbN0,L_c,sigma;

    // for basical info
    char buffer[50];
    const int h = __HEIGHT, w = __WIDTH, f = frame+1 ;
    const int lm = h*w;
    const int lu = lm+(G_L-1);
    int ** G = getGenerator();
    double * Ly = (double*) malloc(sizeof(double)*2*lu);
    int*** Y = new3d<int>(f,h,w);
    int* map_out = (int*) malloc(sizeof(int)*lm);

    double* Lu = (double*) malloc(sizeof(double)*lu);
    for(int i=0; i<lu ;++i)
        Lu[i] = 0;

    // pstate, pout
    const int Ns = pow(2,G_L-1);
    int ** pout = new2d<int>(Ns,4);
    int ** pstate = new2d<int>(Ns,2) ;
    double* Le1 = (double*)malloc(sizeof(double)*lu);
    double* Le2 = (double*)malloc(sizeof(double)*lu);

    trellis(G,G_N,G_L,Ns,pout,pstate);

    // frame buffer
    int*** imgr_bp = new3d<int>(PXL,h,w);
    double** Lu_c = new2d<double>(PXL,lu);  //channel decoder output

    // buffer for Ia, Ie
    double* Le = (double*) malloc(sizeof(double)*lm);
    double* Ia_pc = (double*) malloc(sizeof(double)*snr_size);
    double* Ie_pc = (double*) malloc(sizeof(double)*snr_size);
    int idx = 0;
    double tmp_sum = 0, tmp ;

    // read YUV
    sprintf(buffer,"%s_cif.yuv",argv[1]);
    yuv_read(buffer,h,w,f,Y,NULL,NULL);

//  video_encode(Y,f,h,w,_snr,G,Ly,map_out) for frame and t_lvl
    EbN0 = pow(10,_snr/10);      // convert Eb/N0[dB] to normal number
    L_c = 4*a*EbN0*rate;           // reliability value of the channel
    sigma = 1/sqrt(2*rate*EbN0);   // standard deviation of AWGN noise

    int* x = (int*) malloc(sizeof(int)*2*lu);
    img2bp_frame(Y[frame],h,w,imgr_bp);
    random_sequence(0,lm-1,map_out);
    rsc_encode(G,G_L,imgr_bp[t_lvl],map_out,w,lm,1,x);

    for(int i = 0 ; i < 2*lu ; ++i)
        Ly[i] = 0.5*L_c*((2*x[i] - 1)+ sigma*gaussian_noise());

    printf("processing ...%2.2f%%\n",0);
    for(double snr = s_snr; idx < snr_size ; snr+=step, idx++ , tmp_sum=0){
        printf("\tprocessing ...%2.2f%%\n",100*(snr-s_snr)/step/snr_size);
        // encode
        EbN0 = pow(10,snr/10);      // convert Eb/N0[dB] to normal number
        L_c = 4*a*EbN0*rate;           // reliability value of the channel
        sigma = 1/sqrt(2*rate*EbN0);   // standard deviation of AWGN noise

//        img2bp_frame(Y[frame],h,w,imgr_bp);
        for(int i = 0 ; i < lm ; ++i){
            tmp = (2*imgr_bp[t_lvl][map_out[i]/w][map_out[i]%w] - 1);
            Lu[i] = 0.5*L_c*( tmp + sigma*gaussian_noise());
            tmp_sum+=log2(1 + exp(-Lu[i]*tmp));
        }

        Ia_pc[idx] = 1 - tmp_sum/lm;

        // decode
        computeLe(Lu,Le1,Le2,lu);
        logmap(Ns, lu, 1, Ly, Le1, Le2, pstate, pout, Lu_c[t_lvl]);

        for(int i = 0 ; i < lm ; ++i)
            Le[map_out[i]] = Lu_c[t_lvl][i] - Lu[i];

        tmp_sum = 0;
        for(int i=0; i < lm ; ++i)
            tmp_sum+=log2(1 + exp(-Le[i]*(2*imgr_bp[t_lvl][i/w][i%w]-1)));

        Ie_pc[idx] = 1 - tmp_sum/lm;

    }
    printf("\r100%% completed!\n");

    FILE *file = fopen("output/exit_curve_pc.txt","a+");

    fprintf(file,"============================\n%s:SNR=%lf, frame#%d,bp#%d\nIa=\n",argv[1],_snr,frame+1,t_lvl+1);
    for(int i = 0 ; i < snr_size ; ++i)
        fprintf(file,"%lf,",Ia_pc[i]);
    fprintf(file,"\nIe=\n");
    for(int i = 0 ; i < snr_size ; ++i)
        fprintf(file,"%lf,",Ie_pc[i]);
    fprintf(file,"\n\n");

    fclose(file);

    write_pc_for_matlab(Ia_pc,Ie_pc,snr_size);

    // free memory
    free(Ly);
    free(map_out);
    delete2d<int>(G);
    deleteY(Y);

    delete3d<int>(imgr_bp);
    delete2d<double>(Lu_c);
    delete2d<int>(pstate);
    delete2d<int>(pout);

    free(Le);
    free(Lu);
    free(Le1);
    free(Le2);

    free(Ia_pc);
    free(Ie_pc);
    free(x);
}
