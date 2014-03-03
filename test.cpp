#include <stdio.h>
#include <time.h>
#include "utils.h"
#include "data_alloc.h"
#include "io_utils.h"
#include "video_encode.h"
#include "algs_direct.h"
#include "algs_intra.h"
#include "algs_inter.h"
#include "algs_inter_intra.h"
#include "image_process_utils.h"
#include "channel_code_utils.h"
#include "logmap_util.h"
#include "mrf_decoder_utils.h"
#include "motion_estimation_utils.h"
#include "build_value.h"


// OK
void test_inter_intra(){
    int*** Y;
    double *** Ly;
    int *** map_out;
    int h = 288, w = 352, f = 5 ;
    int lm = h*w;
    int lu = lm+2;
    int snr = 3;
    int ** G = getGenerator();

    Y = yuv_read("stefan_cif.yuv",h,w,f);

    Ly = new3d<double>(f,PXL,2*lu);
    map_out = new3d<int>(f,PXL,lm);

    video_encode(Y,f,h,w,snr,G,Ly,map_out);

    double* PSNR = (double*) malloc(sizeof(double)*f);
    inter_intra("stefan",Y,Ly,map_out,h,w,f,lu,G,PSNR);

    delete3d<double>(Ly);
    delete3d<int>(map_out);
    delete2d<int>(G);
    deleteY(Y);
}


// OK
void test_inter(){
    int*** Y;
    double *** Ly;
    int *** map_out;
    int h = 288, w = 352, f = 5 ;
    int lm = h*w;
    int lu = lm+2;
    int snr = 0;
    int ** G = getGenerator();

    Y = yuv_read("stefan_cif.yuv",h,w,f);

    Ly = new3d<double>(f,PXL,2*lu);
    map_out = new3d<int>(f,PXL,lm);

    video_encode(Y,f,h,w,snr,G,Ly,map_out);

    double* PSNR = (double*) malloc(sizeof(double)*f);
    inter("stefan",Y,Ly,map_out,h,w,f,lu,G,PSNR);

    delete3d<double>(Ly);
    delete3d<int>(map_out);
    delete2d<int>(G);
    deleteY(Y);
}


// OK
void test_intra(){
    int*** Y;
    double *** Ly;
    int *** map_out;
    int h = 288, w = 352, f = 5 ;
    int lm = h*w;
    int lu = lm+2;
    int snr = 0;
    int ** G = getGenerator();

    Y = yuv_read("stefan_cif.yuv",h,w,f);

    Ly = new3d<double>(f,PXL,2*lu);
    map_out = new3d<int>(f,PXL,lm);

    video_encode(Y,f,h,w,snr,G,Ly,map_out);

    double* PSNR = (double*) malloc(sizeof(double)*f);
//    direct_decode("stefan",Y,Ly,map_out,h,w,f,lu,G,PSNR);
    intra("stefan",Y,Ly,map_out,h,w,f,lu,G,PSNR);

    delete3d<double>(Ly);
    delete3d<int>(map_out);
    delete2d<int>(G);
    deleteY(Y);
}


// OK
void test_direct_decode(){
    int*** Y;
    double *** Ly;
    int *** map_out;
    int h = 288, w = 352, f = 5 ;
    int lm = h*w;
    int lu = lm+2;
    int snr = 0;
    int ** G = getGenerator();

    Y = yuv_read("stefan_cif.yuv",h,w,f);

    Ly = new3d<double>(f,PXL,2*lu);
    map_out = new3d<int>(f,PXL,lm);

    video_encode(Y,f,h,w,snr,G,Ly,map_out);

    double* PSNR = (double*) malloc(sizeof(double)*f);
    direct_decode("stefan",Y,Ly,map_out,h,w,f,lu,G,PSNR);

    delete3d<double>(Ly);
    delete3d<int>(map_out);
    delete2d<int>(G);
    deleteY(Y);
}


// OK
void test_img2bp(){
    int imgh = 10;
    int imgw = 10;
    int** img = new2d<int>(imgh,imgw);
    int*** img_bp = new3d<int>(imgh,imgw,PXL);

    for(int i = 0 ; i < imgh ; ++i){
        for(int j = 0 ; j < imgw ; ++j){
            img[i][j] = rand()%256;
            printf("%d ",img[i][j]);
        }
        printf("\n");
    }

    printf("\n\n\n");

    img2bp_frame(img,imgh,imgw,img_bp);

    for(int i = 0 ; i < PXL ; ++i){
        for(int j = 0 ; j < imgh ; ++j){
            for(int k = 0; k < imgw ; ++k)
                printf("%d ",img_bp[j][k][i]);
            printf("\n");
        }
        printf("\n\n");
    }

    delete2d<int>(img);
    delete3d<int>(img_bp);
}

// OK
void test_mrf_param_est(){
    int imgh = 20;
    int imgw = 20;
    int** img = new2d<int>(imgh,imgw);

    for(int i = 0 ; i <imgh ; ++i)
        for(int j = 0 ; j < imgw ; ++j)
                img[i][j] = 1;

    for(int j = 0 ; j < imgw ; ++j)
            img[5][j] = 0;

    double beta;

    intra_beta_estimation(img,beta,imgh,imgw);

    printf("beta=%lf\n",beta);
    delete2d<int>(img);

}

// OK
void test_interleave(){

    int h = 5;
    int w = 8;
    int lm = h*w;
    int* map_out = (int*) malloc(sizeof(int)*lm);
    int** arr = new2d<int>(h,w);
    int* tmp = (int*) malloc(sizeof(int)*lm);

    for(int idx = 0 ; idx < 5 ; ++idx){
        printf("original matrix:\n");
        for(int i = 0 ; i < h ; ++i){
            for(int j = 0 ; j < w ; ++j){
                arr[i][j] = rand()%lm;
                printf("%d\t",arr[i][j]);
            }
            printf("\n");
        }

        printf("\n\nrandom map:\n");
        random_sequence(0,lm-1,map_out);

        for(int i = 0 ; i < lm ; ++i){
            printf("%d ",map_out[i]);
            tmp[i] = arr[map_out[i]/w][map_out[i]%w];
        }

        printf("\n\nmatrix after deinterleave:\n");

        deinterleave(tmp,map_out,lm);

        for(int i = 0 ; i < h ; ++i){
            for(int j = 0 ; j < w ; ++j){
                printf("%d\t",tmp[i*w+j]);
            }
            printf("\n");
        }
        printf("\n\n");

    }

    free(map_out);
    delete2d<int>(arr);
    free(tmp);

}


// OK
void test_psnr(){
    int*** Y;
    int h = 288, w = 352, f = 5 ;
    int** imgr = new2d<int>(h,w);
    int** imgO = new2d<int>(h,w);

    Y = yuv_read("stefan_cif.yuv",h,w,f);

    for(int i = 0 ; i < f ; ++i){
        for(int y = 0 ; y < h ; ++y)
            for(int x = 0 ; x < w ; ++x ){
                imgr[y][x] = Y[y][x][i];
                imgO[y][x] = (imgr[y][x]+5)%256;
            }
        printf("psnr[%d] is %lf \n",i,computePSNR(imgr,imgO,h,w));
    }

    delete2d<int>(imgr);
    delete2d<int>(imgO);
    deleteY(Y);
}

// OK
void test_deci_binl(){
    int M = 8;
    int Ns = pow(2,M);
    int * state;

    for(int i = 0 ; i < Ns ; ++i){
        state = deci2binl(i,M);
        printf("i=%d,and state=",i);
        for(int j = 0 ; j < M ; ++j)
            printf("%d ",state[j]);

        printf("\nand the backward is %d\n",bin2deci(state,M));

        free(state);
    }
}


// OK
void test_bin2dec_img(){
    int h = 2;
    int w = 2;
    int*** img_bp = new3d<int>(h,w,PXL);

    for(int i = 0 ; i < h ; ++i)
        for(int j = 0 ; j < w ; ++j)
            for(int k = 0 ; k < PXL ; ++k)
                img_bp[i][j][k] = rand()%2;

    for(int i = 0 ; i < PXL ; ++i)
        printf("%d ",img_bp[0][0][i]);
    printf("\n\n");

    int** imgr = bin2dec_img(img_bp,h,w);

    for(int i = 0 ; i < h ; ++i){
        for(int j = 0 ; j < w ; ++j){
            printf("%d ",imgr[i][j]);
        }
        printf("\n");
    }

    delete2d<int>(imgr);
    delete3d<int>(img_bp);

}



// OK
void test_encode(){
    int*** Y;
    double *** Ly;
    int *** map_out;
    int h = 288, w = 352, f = 5 ;
    int lm = h*w;
    int lu = lm+2;

    Y = yuv_read("stefan_cif.yuv",h,w,f);

    Ly = new3d<double>(f,PXL,2*lu);
    map_out = new3d<int>(f,PXL,lm);

    video_encode(Y,f,h,w,0,getGenerator(),Ly,map_out);

    int cn=0 ,cp=0;
    for(int i = 0 ; i < 2*lu ; ++i){
        if(Ly[0][0][i]>0)
            cp++;
        else
            cn++;
    }
    printf("positive:%d and nagtive:%d\n",cp,cn);

    delete3d<double>(Ly);
    delete3d<int>(map_out);
    deleteY(Y);
}

// OK
void test_yuv_read(){
    int*** Y;
    int h = 288, w = 352, f = 1;

    Y = yuv_read("stefan_cif.yuv",h,w,f);
    printf("file read done %d\n",Y[0][0][0]);
    for(int t = 0 ; t < f ; ++t){
        for(int i = 50 ; i < 60; ++i){
            for(int j = 50 ; j < 60 ; ++j){
                printf("%d\t", Y[i][j][t]);
            }
            printf("\n");
        }
        printf("\n\n\n");
    }

    deleteY(Y);

}

// OK
void test_random_shuffle(){
    int l = 0, h = 9 , len = h-l+1;
    int * order = (int*) malloc(sizeof(int)*len);
    for(int i = 0 ; i < 10 ; ++i){
        random_sequence(l,h,order);
        for(int j = 0 ; j < len ; ++j)
            printf("%d ",order[j]);
        printf("\n");
    }
}

// OK
void test_gaussain_noise(){
    for(int i = 0 ; i < 10; ++i){
        for(int j = 0 ; j < 10; ++j){
            printf("%lf ",gaussian_noise());
        }
        printf("\n");
    }

}

// OK
void test_trellis(){
    int** G = getGenerator();
    int Ns = 4;
    int ** pout = new2d<int>(Ns,4);
    int ** ps = new2d<int>(Ns,2) ;

    trellis(G,2,3,Ns,pout,ps);

    printf("pout:\n");
    for(int i = 0 ; i < Ns ; ++i){
        for(int j = 0 ; j < 4 ; ++j)
            printf("%d, ",pout[i][j]);
        printf("\n");
    }

    printf("pstate:\n");
    for(int i = 0 ; i < Ns ; ++i){
        for(int j = 0 ; j < 2 ; ++j)
            printf("%d, ",ps[i][j]);
        printf("\n");
    }

    delete2d<int>(G);
    delete2d<int>(pout);
    delete2d<int>(ps);

}


// OK
int M[] = {0,1,0,1,1,1,1,1,0};

void test_rsc_encode(){
    int lm = 9;
    int** G = getGenerator();
    int* m = (int*) malloc(sizeof(int)*lm);
    double* Lu = (double*) malloc(sizeof(double)*(lm+2));
    double* Ly = (double*) malloc(sizeof(double)*2*(lm+2));

    for(int i  = 0 ; i < lm+2 ; ++i)
        Lu[i] = 0;

    for(int i = 0 ; i < lm ; ++i)
        m[i] = M[i];

    int* x = (int*) malloc(sizeof(int)*2*(lm+2));
//    int* x = rsc_encode(G,3,m,lm,1);
    rsc_encode(G,3,m,lm,1,x);

    for(int i = 0 ; i < 2; ++i){
        for(int j = 0 ; j < (lm+2) ; ++j){
            Ly[i+j*2] = x[i+j*2];
            printf("%d ",x[i+j*2]);
        }
        printf("\n");
    }

    delete2d<int>(G);
    free(x);
    free(m);
    free(x);
    free(Lu);
}




// bug? delete2d<int> would crash?!
void test_motionComp(){
    int h = 3;
    int w = 6;
    int mb = 3;
    const int imgI[3][6] = {{0 , 1 , 0 , 0 , 1 , 1},{0, 0 , 1 ,1 ,1 ,1},{0, 0 , 1 ,1 ,0 ,1}};
    const int MV[2][3] = {{0,0,0},{0,0,0}};

    int ** im = new2d<int>(h,w,0);
    int ** comp = new2d<int>(h,w,0);
    int ** mv = new2d<int>(2,h*w/mb/mb);
    for(int i = 0 ; i < h ; ++i)
        for(int j = 0 ; j < w ; ++j){
            im[i][j] = imgI[i][j];
            mv[0][i*w+j] = MV[0][i*w+j];
            mv[1][i*w+j] = MV[1][i*w+j];
        }

    printf("imgI:\n");
    for(int i = 0 ; i < h ; ++i){
        for(int j = 0 ; j < w ; ++j)
            printf("%d ",im[i][j]);
        printf("\n");
    }


    motionComp(im,mv,h,w,mb,comp);

    printf("\ncomp:\n");
    for(int i = 0 ; i < h ; ++i){
        for(int j = 0 ; j < w ; ++j)
            printf("%d ",comp[i][j]);
        printf("\n");
    }

    printf("delete comp\n");
    delete2d<int>(comp);
    printf("delete im\n");
    delete2d<int>(im);
    printf("delete mv\n");
    delete2d<int>(mv);
    printf("done\n");
}

// bug? delete2d<int> would crash?!
void test_LemotionComp(){
    int h = 3;
    int w = 6;
    int mb = 3;
    const double LE[8][18] =
    {{0 , 1 , 0 , 0 , 1 , 1 ,0, 0 , 1 ,1 ,1 ,1 ,0, 0 , 1 ,1 ,1 ,1},{0, 0 , 1 ,1 ,1 ,1, 0 , 1 , 1 ,0, 0 , 1, 0 , 1 , 1 ,0, 0 , 1},
    {0 , 1 , 0 , 0 , 1 , 1 ,0, 0 , 1 ,1 ,1 ,1 ,0, 0 , 1 ,1 ,1 ,1},{0, 0 , 1 ,1 ,1 ,1, 0 , 1 , 1 ,0, 0 , 1, 0 , 1 , 1 ,0, 0 , 1},
    {0 , 1 , 0 , 0 , 1 , 1 ,0, 0 , 1 ,1 ,1 ,1 ,0, 0 , 1 ,1 ,1 ,1},{0, 0 , 1 ,1 ,1 ,1, 0 , 1 , 1 ,0, 0 , 1, 0 , 1 , 1 ,0, 0 , 1},
    {0 , 1 , 0 , 0 , 1 , 1 ,0, 0 , 1 ,1 ,1 ,1 ,0, 0 , 1 ,1 ,1 ,1},{0, 0 , 1 ,1 ,1 ,1, 0 , 1 , 1 ,0, 0 , 1, 0 , 1 , 1 ,0, 0 , 1}};
    const int MV[2][3] = {{0,0,0},{0,0,0}};

    double ** le = new2d<double>(PXL,w*h,0);
    double ** comp = new2d<double>(PXL,w*h,0);
    int ** mv = new2d<int>(2,h*w/mb/mb);
    for(int i = 0 ; i < PXL ; ++i)
        for(int j = 0 ; j < w*h ; ++j){
            le[i][j] = LE[i][j];
        }
    for(int i = 0 ; i<h*w/mb/mb ; ++i){
        mv[0][i] = MV[0][i];
        mv[1][i] = MV[1][i];
    }

    printf("le:\n");
    for(int i = 0 ; i < PXL ; ++i){
        for(int j = 0 ; j < w*h ; ++j)
            printf("%lf ",le[i][j]);
        printf("\n");
    }


    motionComp(le,mv,h,w,mb,comp);

    printf("\ncomp:\n");
    for(int i = 0 ; i < PXL ; ++i){
        for(int j = 0 ; j < w*h ; ++j)
            printf("%lf ",comp[i][j]);
        printf("\n");
    }

    printf("delete comp\n");
    delete2d<double>(comp);
    printf("delete le\n");
    delete2d<double>(le);
    printf("delete mv\n");
    delete2d<int>(mv);
    printf("done\n");
}

// OK
void test_siso_inter(){
    int h = 2, w = 3;

    double** a = new2d<double>(PXL,h*w,2);
    double** b = new2d<double>(PXL,h*w,3);
    double** c = new2d<double>(PXL,h*w);

    double* beta = (double*)malloc(sizeof(double)*PXL);

    for(int i = 0 ; i<PXL ; ++i)
        beta[i] = i*1.5;

    mrf_siso_inter(a,b,beta,h,w,c,1);

    print2DMatrix<double>(c,0,PXL,0,h*w);

    delete2d<double>(a);
    delete2d<double>(b);
    delete2d<double>(c);

    free(beta);
}

// OK
void test_time(){
    for(int i = 0 ; i < 100 ; ++i)
        printf("time:%lf\n",getCurrentTime());
}
// OK
int test_M[3][3] = {{4,3,3},{3,3,2},{3,2,3}};
void test_inv3x3(){

    int** x = new2d<int>(3,3);
    printf("Matrix A=\n");
    for(int i = 0 ,j =0 ; i <3 ; ++i){
        for(j=0 ; j < 3 ; ++j){
            x[i][j] = test_M[i][j];
            printf("%d ",x[i][j]);
        }
        printf("\n");
    }

    printf("\n\nA^-1 =\n");

    double** y = inv3x3<int>(x);

    print2DMatrix<double>(y,0,3,0,3);

    delete2d<int>(x);
    delete2d<double>(y);
}

// OK
void test_yuv_write(){

    int h = 288, w = 352, f = 30 ;
    int*** Y = new3d<int>(f,h,w,0);
    int*** U = new3d<int>(f,h/2,w/2.0);
    int*** V = new3d<int>(f,h/2,w/2,0);

    yuv_read("stefan_cif.yuv",h,w,f,Y,NULL,NULL);
    yuv_write("test",Y,NULL,NULL,f,h,w);

    deleteY(Y);
    deleteY(U);
    deleteY(V);
}


int main(int argc,char* argv[]){
    srand(1);
//    test_random_shuffle();
//    test_encode();
//    test_mrf_param_est();
//    test_yuv_read();
//    test_gaussain_noise();
//    test_bin2dec_img();
//    test_deci_binl();
//    test_psnr();
//    test_interleave();
//    test_trellis();
//    test_rsc_encode();
//    test_direct_decode();
//    test_img2bp();
//    test_intra();
//    test_LemotionComp();
//    test_motionComp();
//    test_siso_inter();
//    test_inter();
//    test_time();
//    test_inv3x3();
//    test_inter_intra();
//    test_yuv_write();


    return 0;
}

