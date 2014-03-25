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
//    double *** Ly = new3d<double>(300,8,2*288*352);
    weight_predict_minMSE(weights,pow(10,-1.0/10));

    free(weights);
//    delete(alloc);
//    free(alloc);
//    delete3d<double>(Ly);
    return 0;
}
