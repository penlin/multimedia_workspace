#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "rsc_encode.h"
#include "utils.h"
#include "data_alloc.h"
#include "channel_code_utils.h"
#include "image_process_utils.h"
#include "build_value.h"


/*
*  @Penlin: main function to encode Y into Ly and map_out through channel (Gaussian noise)
*
*  @param: Y            Y value of the video through all frame  [n_frame*imgh*imgw]
*  @param: EbN0dB       a const value for SNR in dB unit        [1]
*  @param: Ly           output log likelihood matrix            [n_frame*PXL*(2*lu)]
*  @param: map_out      bits index map for interleave           [n_frame*PXL*(imgh*imgw)]
*
*/
void video_encode(int*** Y,const int &n_frame,const int &imgh, const int &imgw, const double &EbN0dB, int** G_ptr, double *** Ly , int*** map_out){

    /* Channel code parameter  */

    //parameters
    const int puncture = 0;                     // puncture or not
    double rate = 1/(double)(2-puncture);       // code rate
    const double a = 1;                         // Fading amplitude. a=1 -> AWGN channel
    const int lm = imgh*imgw;                   // Frame size
    const int lu = lm + 2;                      // Length of message bit sequence
    const double EbN0 = pow(10,EbN0dB/10);      // convert Eb/N0[dB] to normal number
    const double L_c = 4*a*EbN0*rate;           // reliability value of the channel
    const double sigma = 1/sqrt(2*rate*EbN0);   // standard deviation of AWGN noise

    // outputfile pre-allocating
    int **img;
    int ***img_bp = new3d<int>(PXL,imgh,imgw);

    int *x = (int*) malloc(sizeof(int)*2*lu);

    for (int f = 0 ; f < n_frame ; ++f){
        printf("Encoding frame#%d\n",f+1);

        img = Y[f];
        img2bp_frame(img,imgh,imgw,img_bp);

        for(int t_lvl = 0 ; t_lvl < PXL ; ++t_lvl){
            // interleave
            random_sequence(0,lm-1,map_out[f][t_lvl]);

            rsc_encode(G_ptr,G_L,img_bp[t_lvl],map_out[f][t_lvl],imgw,lm,1,x);

            for(int i = 0 ; i < 2*lu ; ++i){
                Ly[f][t_lvl][i] = 0.5*L_c*((2*x[i] - 1)+ sigma*gaussian_noise());  // add noise   + sigma*gaussian_noise()
            }
        }
    }

    free(x);
    delete3d<int>(img_bp);

}

