#include <stdio.h>
#include <stdlib.h>
#include "data_types.h"
#include "data_alloc.h"

void test2();
int test1();

int main(){
    test2();
    return 0;
}

void test2(){
    size_t w = 352 , h = 288, PXL = 8, f = 20, mbSize = 8, Ns = 4;
    size_t lm = w*h;
    size_t lu = lm+2;
    initMemMgr(w*h*PXL*sizeof(double)*5);
    for(int i = 0 ; i < 1000; ++i){
        double*** a = new3d<double>(PXL,lm,4);
//        DELETE(a[0][0]);
//        DELETE(a[0]);
        DELETE(a);
    }

    printMem();
    freeMemMgr();
}

int test1(){

    size_t w = 352 , h = 288, PXL = 8, f = 20, mbSize = 8, Ns = 4;
    size_t lm = w*h;
    size_t lu = lm+2;
    initMemMgr(w*h*PXL*sizeof(double)*5);

    double** p0 = new2d<double>(h,w,0.5);
//    int** G =new2d<int>(2,3,2);
//    int ** pout = new2d<int>( 4,4);
//    int ** pstate = new2d<int>(4,2) ;
//    int* y = MALLOC(int,2);
//    double* PSNR = MALLOC(double,f);
//    int** out = new2d<int>(2,2);
//    int** Y = new2d<int>(h,w,0);
//    int*** img_bp = new3d<int>(PXL,h,w,2);
//    int* x = MALLOC(int,2*lu);
//    int** U =  new2d<int>(h/2,w/2,0);
//    pixel*  buffer = MALLOC(pixel,lm*3/2);
//    double** Ly = new2d<double>(PXL,2*lu,2);
//    double** Lu_c = new2d<double>(PXL,lu);
//    double** Le_c = new2d<double>(PXL,lm);
//    int** map = new2d<int>(PXL,lm);
//    double* beta = MALLOC(double,PXL);
//    int** MV = new2d<int>(2,lm/(mbSize*mbSize),0);
//    double** Alpha = new2d<double>(lu,Ns,0);
//    double** Beta = new2d<double>(lu+1,Ns,0);
//    double*** gam = new3d<double>(Ns,Ns,lu+1,1);
//    int** hash_key = new2d<int>(h-2,w-2);
//    double** x3 = new2d<double>(3,3);


    printMem();

    DELETE(p0);
//    DELETE(y);
//    DELETE(PSNR);
//    DELETE(x);
//    DELETE(buffer);
//    DELETE(beta);
//    delete2d<int>(G);
//    delete2d<int>(pout);
//    delete2d<int>(pstate);
//    delete2d<int>(out);
//    delete2d<int>(Y);
//    delete2d<int>(U);
//    delete2d<int>(map);
//    delete2d<int>(MV);
//    delete2d<int>(hash_key);
//    delete2d<double>(Ly);
//    delete2d<double>(Lu_c);
//    delete2d<double>(Le_c);
//    delete2d<double>(Alpha);
//    delete2d<double>(Beta);
//    delete2d<double>(x3);
//    delete3d<int>(img_bp);
//    delete3d<double>(gam);

    printMem();
    freeMemMgr();
    return 0;
}
