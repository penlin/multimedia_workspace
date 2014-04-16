#include "frame.h"
#include "uep_predict_utils.h"

int main(int argc,char* argv[]){

    double* weights = (double*)malloc(sizeof(double)*PXL);
    double snr = 0.0;
    FILE* fptr = fopen(__SEQ__,"r+b");
    Frame* frame = new Frame(__HEIGHT,__WIDTH,0,0);
    frame->read(fptr);
    if(argc>1)
        snr = atof(argv[1]);
    weight_predict_Intra_minMSE(frame->img_bp,__HEIGHT,__WIDTH,weights,pow(10,snr/10));

    free(weights);
    delete frame;
    fclose(fptr);
    return 0;
}
