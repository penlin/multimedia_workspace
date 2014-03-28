#include <stdio.h>
#include <stdlib.h>
#define PXL 8
#include "uep_predict_utils.h"
#include <math.h>
#include "data_alloc.h"


int main(int argc,char* argv[]){

//    double* alloc = (double*)malloc(sizeof(double)*486604800*10);
//    double* alloc = new double[486604800];
    double* weights = (double*)malloc(sizeof(double)*PXL);
    double value, tmpSum = 0;
//    double *** Ly = new3d<double>(300,8,2*288*352);
    weight_predict_minMSE(weights,pow(10,4.0/10));
//    weight_predict_minPOWER(weights,205);
//    for(int i = 0 ; i < PXL ; ++i){
//        cut(weights[i],MAX_GAMMA,MIN_GAMMA);
//        interp2(weights[i],value);
//        printf("f(r%d=%lf)=%lf , contribution:%lf\n",i,weights[i],value,value/ORDER2[i]);
//        tmpSum+=(value/ORDER2[i]);
//    }
//
//    printf("sum = %lf\n",tmpSum);
    free(weights);
//    delete(alloc);
//    free(alloc);
//    delete3d<double>(Ly);
    return 0;
}
