#ifndef __ALGS_DIRECT_H
#define __ALGS_DIRECT_H

#include <stdio.h>
#include <stdlib.h>
#include "data_alloc.h"
#include "build_value.h"


/**
*   @Penlin: algorithm for direct decode
*
*   @param: Ly_in   log likelihood value of bit planes                  [n_frame*PXL*(2*lu)]
*   @param: map_in  mapping for deinterleave                            [n_frame*PXL*lm]
*
*   @param: lu, G
*
**/

void direct_decode(int* bitstream, double* Ly, int* map_out, const int &lm, const int &lu, int** pstate, int** pout, int &n_err){

    // param initial
    const int Ns = pow(2,G_L-1);
    double* Le1 = (double*)malloc(sizeof(double)*lu);
    double* Le2 = (double*)malloc(sizeof(double)*lu);

    for(int i=0; i<lu ;++i)
        Le1[i]  = Le2[i] = 0.5;

    double* Lu_c = (double*)malloc(sizeof(double)*lu);  //channel decoder output

    // decoding
    BCJR_decoding(Ns, lu, 1, Ly, Le1, Le2, pstate, pout, Lu_c);

    // deinterleave
#if __STATUS__
    printf("deinterleave ...%lf\n",getCurrentTime());
#endif
    deinterleave(Lu_c,map_out,lm);

    for(int i = 0 ; i < lm ; ++i){
        if((2*bitstream[i]-1)*Lu_c[i] < 0 )
            n_err += 1;
    }

    free(Lu_c);
    free(Le1);
    free(Le2);
}


#endif // __ALGS_DIRECT_H
