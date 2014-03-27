#ifndef __MOTION_ESTIMATION_SOFT_UTILS_H
#define __MOTION_ESTIMATION_SOFT_UTILS_H
#include <stdio.h>
#include <stdlib.h>

#define BOUND 0.01

//static const int ORDERS[8] = {128 , 64, 32, 16, 8, 4, 2, 1};

int llr_bp_to_img(const double& llr, const int &t_lvl){

    return ORDERS[t_lvl]*((llr>=BOUND)?1:
                          (-llr>=BOUND)?0:
                          (llr>=0)?0.67:
                          0.33);
}


#endif // __MOTION_ESTIMATION_SOFT_UTILS_H
