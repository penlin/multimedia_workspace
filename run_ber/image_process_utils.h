#ifndef __IMAGE_PROCESS_UTILS_H
#define __IMAGE_PROCESS_UTILS_H

#define PXL 8

int ORDER[PXL] = {1 , 2, 4, 8 ,16 ,32, 64, 128};


int* deci2binl(int x, const int &l){
    int* y = (int*)malloc(sizeof(int)*l);
    int bit = 0;

    while(x>=1){
        y[l-(bit++)-1] = x%2;
        x/=2;
    }

    for(int i = 0 ; i < (l-bit) ; ++i )
        y[i] = 0;

    return y;
}

/*
% z=bin2deci(x)
% e.g. x=[ 1     1     0     1     0     0     1     0     0     0
%          1     1     0     1     1     0     0     0     1     0]'
% z =
*
*  @Penlin
*  @param x  [A*B]
*/

int* bin2deci(int** x, const int &A, const int &B){
    int* z = (int*)malloc(sizeof(int)*A);
    for(int i = 0 ; i < A; ++i){
        z[i] = 0;
        for(int j = 0 ; j < B ; ++j)
            z[i] += ORDER[j]*x[i][B-j-1];
    }

    return z;
}

int bin2deci(int* x, const int &B){
    int z = 0 ;
    for(int j = 0 ; j < B ; ++j)
        z += ORDER[j]*x[B-j-1];

    return z;
}


double computePSNR(int** imgr, int** imgO, const int &imgh, const int &imgw){
    double psnr = 0;
    int tmp = 0;

    for(int i = 0 ; i < imgh ; ++i){
        for(int j = 0 ; j < imgw ; ++j){
            tmp = imgr[i][j] - imgO[i][j];
            psnr += (tmp*tmp);
        }
    }

    psnr = psnr/imgh/imgw;
    if(psnr>0)
        psnr = 10*log10(65025/psnr); // psnr = 10*log10( (2^PXL-1)^2 / psnr )
    else
        psnr = -1;   // Inf

    return psnr;
}


#endif // __IMAGE_PROCESS_UTILS_H
