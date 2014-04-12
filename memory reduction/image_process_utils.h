#ifndef __IMAGE_PROCESS_UTILS_H
#define __IMAGE_PROCESS_UTILS_H

#define PXL 8

int ORDER[PXL] = {1 , 2, 4, 8 ,16 ,32, 64, 128};


 /* @Penlin : tranform a pixels frame into 8-bit planes
 *  @param: img     input frame         [imgh*imgw]
 *  @param: img_bp  output bit planes   [PXL*imgh*imgw]
 */

void img2bp_frame(int ** img, const int &imgh, const int &imgw, int *** img_bp){
    int mask = 1, pixel;

    for(int y = 0; y < imgh ; ++y){
        for(int x = 0 ; x < imgw ; ++x){
            pixel = (int) img[y][x];
            for(int j=PXL-1 ; j >= 0 ; --j){
                img_bp[j][y][x] = (int)((pixel & mask)?1:0);
                pixel = ( pixel >> 1 );
            }
        }
    }
}

void img2bp_frame_dpcm(int ** img, const int &imgh, const int &imgw, int *** img_bp){
    int mask = 1, pixel;

    for(int y = 0; y < imgh ; ++y){
        for(int x = 0 ; x < imgw ; ++x){
            pixel = (int)img[y][x];
            for(int j=PXL-1 ; j >= 0 ; --j){
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

int** bin2dec_img(int*** img_bp, const int &imgh, const int& imgw){

    int ** imgr  = new2d<int>(imgh,imgw);

    for(int i=0 ; i<imgh ; ++i){
        for(int j=0 ; j<imgw; ++j){
            imgr[i][j] = 0;
            for(int k=PXL-1;k>=0;--k){
                imgr[i][j] += ORDER[PXL-k-1]*img_bp[k][i][j];
 //               printf("%d ",img_bp[i][j][k]);
            }
        }
    }
//    printf("\n\n\n");
    return imgr;
}

void bin2dec_img(int*** img_bp, const int &imgh, const int& imgw, int** imgr){

    for(int i=0 ; i<imgh ; ++i){
        for(int j=0 ; j<imgw; ++j){
            imgr[i][j] = 0;
            for(int k=PXL-1;k>=0;--k)
                imgr[i][j] += ORDER[PXL-k-1]*img_bp[k][i][j];
        }
    }
}

void bin2dec_img_dpcm(int*** img_bp, const int &imgh, const int& imgw, int** imgr){

    for(int i=0 ; i<imgh ; ++i){
        for(int j=0 ; j<imgw; ++j){
            imgr[i][j] = 0;
            for(int k=PXL-1;k>0;--k)
                imgr[i][j] += ORDER[PXL-k-1]*img_bp[k][i][j];
            imgr[i][j]-=128*img_bp[0][i][j];
            imgr[i][j]*=2;
        }
    }
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
