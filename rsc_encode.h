/*
function x = rsc_encode(G,m,termination)
% Copyright 1998, Yufei Wu, MPRG lab, Virginia Tech. for academic use
% encodes a binary data block m (0/1) with a RSC (recursive systematic
% convolutional) code defined by generator matrix G, returns the output
% in x (0/1), terminates the trellis with all-0 state if termination>0
if nargin<3, termination = 0; end
[N,L] = size(G); % Number of output bits, Constraint length
M = L-1; % Dimension of the state
lu = length(m)+(termination>0)*M; % Length of the input
lm = lu-M; % Length of the message
state = zeros(1,M); % initialize the state vector
% To generate the codeword
x = [];
for i = 1:lu
   if termination<=0 || (termination>0 && i<=lm)
     d_k = m(i);
   elseif termination>0 && i>lm
     d_k = rem(G(1,2:L)*state.',2);
   end
   a_k = rem(G(1,:)*[d_k state].',2);
   xp = rem(G(2,:)*[a_k state].',2); % 2nd output (parity) bits
   state = [a_k state(1:M-1)]; % Next sttate
   %disp(state)
   x = [x [d_k; xp]]; % since systematic, first output is input bit
end
*/
#ifndef __RSC_ENCODE_H
#define __RSC_ENCODE_H
#include <stdio.h>
#include <stdlib.h>
#include "utils.h"



// updated version, direct use imgr_bp and the random map instead block
void rsc_encode(int ** G, const int &L, PIXEL** imgr, const int &mask, const int* map, const size_t &lm, const int &termination, uint8* output){

    int M = L-1, lu, i, j;
    uint8 *state = MALLOC(uint8,M);
    uint8 d_k = 0, a_k = 0, xp=0;

    for(i=0;i<M;++i)
        state[i] = 0;

    lu = (lm + ((termination>0)?M:0));

    for(i=0;i<lu;++i){
        if(termination <= 0 || ( termination > 0 && i < lm ))
            d_k = ((imgr[0][map[i]]&mask) > 0);
        else if(termination > 0 && i >= lm)
            d_k = InnerProductUint8(G[0], state, 1 ,L-1,0, M-1)%2;

        a_k = (G[0][0]*d_k+InnerProductUint8(G[0],state,1,L-1,0,M-1))%2;
        xp =  (G[1][0]*a_k+InnerProductUint8(G[1],state,1,L-1,0,M-1))%2;

        for(j=M-1;j>0;--j)
           state[j] = state[j-1];

        state[0] = a_k;

        output[2*i] = d_k;
        output[1+2*i] = xp;

    }

    DELETE(state);
}


#endif // __RSC_ENCODE_H
