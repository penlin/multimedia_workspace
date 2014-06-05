#include "build_value.h"
#include "data_types.h"
#include "io_utils.h"
#include "uep_predict_utils.h"

const double Distortion[6] = {2072.430287, 655.36, 65.536, 6.5536, 0.65536, 0.065536};


int main(int argc,char* argv[]){



    const int len = __SNR_E-__SNR_S+1;
    double snr[len];
    double mse[len];

    double* weights = MALLOC(double,PXL);

    for(int i = 0 ; i < len ; ++i)
        mse[i] = Distortion[i];//__SNR_S+i;

    int weight_type = 0;
    // 0: EEP
    // 1: UEP
    if(argc > 1)
        weight_type = atoi(argv[1]);

    for(int i = 0 ; i < len; ++i){

        if(weight_type){
            weight_predict_minPOWER(weights,mse[i]);
            snr[i] = 0.0;
            for(int j = 0 ; j < PXL ; ++j)
                snr[i] += weights[j];
            snr[i] /= PXL;
            for(int j = 0 ; j < PXL ; ++j)
                weights[j]/=snr[i];
            snr[i] = 10*log10(snr[i]);

        }else{
            for(int j = 0 ; j < PXL ; ++j)
                weights[j] = 1;
            interpf(mse[i]/21845,snr[i]);
            snr[i] = 10*log10(snr[i]);
        }

        printf("Gamma_0 = %lf dB :[",snr[i]);
        for(int j = 0 ; j < PXL ; ++j)
            printf("%lf, ",weights[j]);
        printf("]\n");

    }


    // free memory
    DELETE(weights);

}
