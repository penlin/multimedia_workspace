#ifndef __CHANNEL_CODE_H
#define __CHANNEL_CODE_H

#if __BCJR__==__MAP__
#include "map_util.h"
#else
#include "logmap_util.h"
#endif

#include "math.h"
#include "build_value.h"

#define G_N 2
#define G_L 3
#define G_M (G_L-1)

const static int G[G_N][G_L] = {{1,1,1},{1,0,1}};      // code generator
static int Ns;
static int** G_ptr;
static int** pout;
static int** pstate;



/**
* Copyright 1998, Yufei Wu, MPRG lab, Virginia Tech. for academic use
* Log_MAP algorithm using straightforward method to compute branch cost
* Input: Ly     = scaled received bits Ly=0.5*L_c*y=(2*a*rate*Eb/N0)*y
*        G      = code generator for the RSC code a in binary matrix
*        Lu     = extrinsic information from the previous decoder.
*        ind_dec= index of decoder=1/2 (assumed to be terminated/open)
* Output: L_A   = ln (P(x=1|y)/P(x=-1|y)), i.e., Log-Likelihood Ratio
*              (soft-value) of estimated message input bit at each level
*tic
*
*   @Penlin: modified to C++ version
*
*   @param Ly
*   @param G
*   @param Lu
*   @param ind_dec
*   @param lu, L
*
*/

void BCJR_decoding(const int &lu, const int &ind_dec, double* Ly, double* Le1, double* Le2, double* LA){

#if __BCJR__==__MAP__
    BCJR_map(Ns, lu, ind_dec, Ly,  Le1, Le2,  pstate, pout, LA);
#else
    BCJR_logmap(Ns, lu, ind_dec, Ly,  Le1, Le2,  pstate, pout, LA);
#endif

}

void trellis(){

    int8* state_b = MALLOC(int8,G_M);
    int8 d_k, a_k, bN, b1;
    int8** out = new2d<int8>(2,2);
    int8** state = new2d<int8>(2,G_M);
    int8** nstate = new2d<int8>(Ns,2);
    int8** nout = new2d<int8>(Ns,4);

    // Set up next_out and next_state matrices for RSC code generator G
    for(int state_i = 0,i = 0, input_bit = 0 ; state_i < Ns ; ++state_i){
        deci2binl(state_i,G_M,state_b);
        for(i = 1 ; i < G_M ; ++i )
            state[0][i] = state[1][i] = state_b[i-1];


        for(input_bit = 0 ; input_bit <= 1; ++input_bit){
            d_k = input_bit;
            a_k = (G_ptr[0][0]*d_k + InnerProduct(G_ptr[0],state_b,1,G_M,0,G_M-1))%2;
            out[d_k][0] = d_k;
            out[d_k][1] = (G_ptr[1][0]*a_k + InnerProduct(G_ptr[1],state_b,1,G_M,0,G_M-1))%2;

            state[d_k][0] = a_k;
        }

        for(i = 0 ; i < 4; ++i)
            nout[state_i][i] = 2*out[i/2][i%2] -1;

        bin2deci(state,2,G_M,nstate[state_i]);
    }

    // Possible previous states having reached the present state
    // with input_bit=0/1
    for(int input_bit = 0, state_i = 0, i = 0 ; input_bit < 2; ++input_bit){
        bN = input_bit*G_N;
        b1 = input_bit;
        for(state_i = 0 ; state_i < Ns ; ++state_i){
            pstate[nstate[state_i][b1]][b1] = (int) state_i;
            for(i = bN ; i < (G_N+bN) ; ++i)
                pout[nstate[state_i][b1]][i] = (int) nout[state_i][i];
        }
    }


    delete2d<int8>(out);
    delete2d<int8>(state);
    delete2d<int8>(nstate);
    delete2d<int8>(nout);
    DELETE(state_b);
}

int** getGeneratorPrepare(){
    G_ptr = new2d<int>(G_N,G_L);
    for(int i = 0,j = 0 ; i < G_N ; ++i)
        for(j = 0 ; j < G_L ; ++j)
            G_ptr[i][j] = G[i][j];

    Ns = (int) pow(2,G_M);
    pout = new2d<int>(Ns,4);
    pstate = new2d<int>(Ns,2);
    trellis();
    return G_ptr;
}

void clearChannel(){
    delete2d<int>(G_ptr);
    delete2d<int>(pstate);
    delete2d<int>(pout);
}

/**
*   @Penlin: 2D
*/

void deinterleave(double** Lu, int** map, const int &lm){
#if __INTERLEAVE__
    double* tmp = MALLOC(double,lm);//(double*) malloc(sizeof(double)*lm);

    for(int t_lvl = 0 ; t_lvl < PXL; ++t_lvl){
        for(int i = 0 ; i < lm ; ++i)
            tmp[i] = Lu[t_lvl][i];

        for(int i = 0 ; i < lm ; ++i)
            Lu[t_lvl][map[t_lvl][i]] = tmp[i];
    }

    DELETE(tmp);//free(tmp);
#endif
}

/** @depriciated 2013/07/23
*   @Penlin: 1D test version
*/

void deinterleave(int* Lu, int* map, const int &lm){
#if __INTERLEAVE__
    int* tmp = MALLOC(int,lm);//(int*) malloc(sizeof(double)*lm);

    for(int i = 0 ; i < lm ; ++i)
        tmp[i] = Lu[i];

    for(int i = 0 ; i < lm ; ++i)
        Lu[map[i]] = tmp[i];

    DELETE(tmp);//free(tmp);
#endif
}


/**
*   @Penlin: 2D interleave
*
*/

void interleave(double** Lu, int** map, const int &lm ){
#if __INTERLEAVE__
    double* tmp = MALLOC(double,lm);//(double*) malloc(sizeof(double)*lm);

    for(int t_lvl = 0 ; t_lvl < PXL; ++t_lvl){
        for(int i = 0 ; i < lm ; ++i)
            tmp[i] = Lu[t_lvl][map[t_lvl][i]];

        for(int i = 0 ; i < lm ; ++i)
            Lu[t_lvl][i] = tmp[i];
    }

    DELETE(tmp);//free(tmp);
#endif
}

#endif // __CHANNEL_CODE_H
