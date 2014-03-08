#include "build_value.h"
#include "io_utils.h"
#include "data_alloc.h"
#include "utils.h"
#include "rsc_encode.h"
#include "channel_code_utils.h"
#include "algs_direct.h"

#define DECODE __ALGO__
#define BCJR 0

void generateBitstream(int* bitstream, const int &lm, int* map, int* encoded_sequence, int** G){

    random_sequence(0,lm-1,map);
    for(int i = 0 ; i < lm ; ++i)
        bitstream[i] = rand()%2;
#if BCJR
    rsc_encode(G,G_L,bitstream,map,lm,1,encoded_sequence);
#endif // BCJR

}

void AWGN(int* bitstream, const int &lm, int* map, double* receive_sequence, const double &EbN0dB){
    const double a = 1;                         // Fading amplitude. a=1 -> AWGN channel
    const double EbN0 = pow(10,EbN0dB/10);      // convert Eb/N0[dB] to normal number
    const double sigma = 1/sqrt(2*EbN0);          // standard deviation of AWGN noise

    for(int i = 0 ; i < lm ; ++i)
        receive_sequence[i] = a*(2*bitstream[map[i]] - 1)+ sigma*gaussian_noise();  // add noise   + sigma*gaussian_noise()
}


int main(int argc,char* argv[]){

    startRandom();

    const int h = __HEIGHT, w = __WIDTH;
    const int lm = h*w;
    double snr = __SNR;
    int n_err = 0 , Nbit = 0;
    int * map_out = (int*)malloc(sizeof(int)*lm);
    int * bitstream = (int*)malloc(sizeof(int)*lm);

# if BCJR
    const int lu = lm+(G_L-1);
    int ** G = getGenerator();
    int * x  = (int*)malloc(sizeof(int)*2*lu);
    double * Ly = (double*)malloc(sizeof(double)*2*lu);
    const int Ns = pow(2,G_L-1) ;
    int ** pout = new2d<int>(Ns,4);
    int ** pstate = new2d<int>(Ns,2) ;
    trellis(G,G_N,G_L,Ns,pout,pstate);
#else
    double* y = (double*)malloc(sizeof(double)*lm);
#endif

    while(Nbit < 1000000 || n_err < 10){
        // encode
#if BCJR
        generateBitstream(bitstream,lm,map_out,x,G);
        generate_Ly(x,lu,snr,Ly);
#else
        generateBitstream(bitstream,lm,map_out,NULL,NULL);
        AWGN(bitstream,lm,map_out,y,snr);
#endif


        // decode
        Nbit += lm ;
#if BCJR
        DECODE(bitstream,Ly,map_out,lm,lu,pstate,pout,n_err);
#else
        for(int i = 0 ; i < lm ; ++i){
            if((2*bitstream[i]-1)*y[i] < 0 )
                n_err += 1;
        }
#endif
    }
    printf("SNR=%.1f dB, BER=%.5f (%d,%d)\n",snr,(double)n_err/Nbit,n_err,Nbit);

    // free memory
    free(map_out);
    free(bitstream);
#if BCJR
    free(x);
    free(Ly);
    delete2d<int>(G);
    delete2d<int>(pstate);
    delete2d<int>(pout);
#else
    free(y);
#endif // BCJR
}
