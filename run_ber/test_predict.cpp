#include <stdio.h>
#include <stdlib.h>
#define PXL 8
#include "uep_predict_utils.h"


int main(int argc,char* argv[]){

    double* weights = (double*)malloc(sizeof(double)*PXL);
    weight_predict_minMSE(weights,1);

    free(weights);
    return 0;
}
