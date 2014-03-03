#ifndef __BUILD_VALUE_H
#define __BUILD_VALUE_H

/*********************
*** CONSTANT VALUE ***
**********************/
#define PXL 8

// decode algorithm function name
#define __DIRECT__      direct_decode
#define __INTRA__       intra
#define __INTER__       inter_test
#define __INTER_INTRA__ inter_intra

// output sequence type
#define __YUV420__    1
#define __Y__         2

// BCJR decoding algorithm
#define __MAP__     1
#define __LOG_MAP__ 2


/************************
*** BASIC INFORMATION ***
*************************/

// basic decoding information
#define Niter  2                        // #iteration
#define __HEIGHT    288
#define __WIDTH     352
#define __FRAME     40
#define __SNR       0
#define __SNR_S     0
#define __SNR_E     8

/*********************
*** DECODE CONTROL ***
**********************/

// control the decode algorithm
#define __ALGO__    __DIRECT__

#define __RANDOM__      1               // total control for the interleave and AWGN
#define __INTERLEAVE__  1&&__RANDOM__   // generate random order map
#define __NOISE__       1&&__RANDOM__   // generate gaussian noise

#define __DEBUG__       1               // debug msg
#define __STATUS__      0&&__DEBUG__    // msg for status right now
#define __BETA__        0&&__DEBUG__    // msg for beta msg
#define __PSNR__        1&&__DEBUG__    // msg for cPSNR in procedure

#define __OUTPUT_SEQ__      1           // control if output the decoded sequence
#define __OUTPUT_TYPE__     __Y__*__OUTPUT_SEQ__

#define __EXIT_INFO__       0

#define __BCJR__    __MAP__

#endif // __BUILD_VALUE_H
