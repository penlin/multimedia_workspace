#ifndef __IMAGE_PROCESS_UTILS_H
#define __IMAGE_PROCESS_UTILS_H

#define PXL 8

int ORDER[PXL] = {1 , 2, 4, 8 ,16 ,32, 64, 128};


 /* @Penlin : tranform a pixels frame into 8-bit planes
 *  @param: img     input frame         [imgh*imgw]
 *  @param: img_bp  output bit planes   [PXL*imgh*imgw]
 */

void img2bp_frame(PIXEL ** img, const size_t &imgh, const size_t &imgw, int *** img_bp){
    int mask = 1, pixel;
    for(int y = 0,x = 0,j = 0; y < imgh ; ++y){
        for(x = 0 ; x < imgw ; ++x){
            pixel = (int) img[y][x];
            for(j=PXL-1 ; j >= 0 ; --j){
                img_bp[j][y][x] = (int)((pixel & mask)?1:0);
                pixel = ( pixel >> 1 );
            }
        }
    }
}

/*
*   @Penlin: Converts a given decimal number x into a binary number of l bits
*
*   @param  x
*   @param  l
*   @output y
*/
//
void deci2binl(const int &x, const int &l, uint8* y){

    for(int i = l-1, bit = 1 ; i >= 0 ; --i, bit <<=1 )
        y[i] = ((x&bit) > 0);
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

void bin2deci(uint8** x, const int &A, const int &B, uint8* z){
    for(int i = 0 , j = 0; i < A; ++i){
        z[i] = 0;
        for(j = 0 ; j < B ; ++j)
            z[i] += ORDER[j]*x[i][B-j-1];
    }
}


int bin2deci(int* x, const int &B){
    int z = 0 ;
    for(int j = 0 ; j < B ; ++j)
        z += ORDER[j]*x[B-j-1];

    return z;
}


void Lu2dec_img(double** Lu, const size_t &lm, PIXEL** imgr){
    for(size_t i = 0 ; i < lm ; ++i)
        imgr[0][i] = (PIXEL)((Lu[0][i]>=0)*ORDER[PXL-1]);

    for(size_t i=0,k=1 ;k <PXL ; ++k){
        for( i=0 ; i < lm; ++i){
//            printf("#(i,k)=(%d,%d)\n",i,k);
            imgr[0][i] = (PIXEL)(imgr[0][i]+(Lu[k][i]>=0)*ORDER[PXL-k-1]);
        }
    }
}

void Lu2dec_img(double** Lu, const size_t &lm, PIXEL** imgr, int** map_out){
    size_t  i = 0;
    for(i = 0 ; i < lm ; ++i)
        imgr[0][map_out[0][i]] = (PIXEL)((Lu[0][i]>=0)*ORDER[PXL-1]);

    for(int t_lvl =PXL-1; t_lvl > 0 ; --t_lvl)
        for(i = 0 ; i < lm ; ++i)
            imgr[0][map_out[t_lvl][i]] = (PIXEL) (imgr[0][map_out[t_lvl][i]]+(Lu[t_lvl][i]>=0)*ORDER[PXL-t_lvl-1]);

}

double computePSNR(PIXEL** imgr, PIXEL** imgO, const size_t &imgh, const size_t &imgw){
    double psnr = 0;
    int tmp = 0;

    for(size_t i = 0,j = 0 ; i < imgh ; ++i){
        for(j = 0 ; j < imgw ; ++j){
            tmp = (int)imgr[i][j] - (int)imgO[i][j];
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

void computeBER(int** imgr, int** imgO, const int &lm, double* ber){

    int mask = 1, value;
    for(int i = 0 ; i < PXL ; ++i)
        ber[i] = 0.0;

    for(int i = 0 ; i < lm ; ++i){
        value = imgr[0][i] ^ imgO[0][i];
        for(int t_lvl = 0; t_lvl < PXL ; ++t_lvl){
            ber[PXL-1-t_lvl] += ((value & (1<<t_lvl))>0);
        }
    }

    for(int i = 0 ; i < PXL ; ++i){
        ber[i] /= lm;
        printf("%lf, ",ber[i]);
    }
}

#endif // __IMAGE_PROCESS_UTILS_H
