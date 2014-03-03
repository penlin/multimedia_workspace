#ifndef __MRF_DECODER_UTILS_H
#define __MRF_DECODER_UTILS_H

#include "beta_estimation_utils.h"
#include "motion_estimation_utils.h"
#define PXL 8

/**
*   @Penlin : MRF model beta estimation for intra decoder
*
*   @param imgr_bp      [PXL*imgh*imgw]
*   @param beta         [n_frame*PXL]
*
*/

// main function
void intra_beta_estimation(int*** imgr_bp , double* beta, const int &imgh, const int &imgw){

    for(int t_lvl = 0 ; t_lvl < PXL ; ++t_lvl){
        intra_beta_estimation(imgr_bp[t_lvl],beta[t_lvl],imgh,imgw);

    }

}


/**
*   @Penlin:  beta = inter_beta_estimation(img_bp, img_bp_ref)
*
*   @param imgr_bp, beta, imgr_bp_ref , motionVect
*/
// main function
void inter_beta_estimation(int*** imgr_bp , int*** img_bp_ref, int** motionVect, double* beta, const int &imgh, const int &imgw, const int &mbSize){

    int ** img_ref = new2d<int>(imgh,imgw);

    for(int t_lvl = 0 ; t_lvl < PXL ; ++t_lvl){

        motionComp(img_bp_ref[t_lvl],motionVect,imgh,imgw,mbSize,img_ref);

        inter_beta_estimation(imgr_bp[t_lvl],beta[t_lvl],img_ref,imgh,imgw);

    }

    delete2d<int>(img_ref);
}

/**
*   @Penlin:  beta = inter_beta_estimation(img_bp, img_bp_ref)
*
*   @param imgr_bp, beta_t, beta_t2, imgr_bp_ref, imgr_bp_ref2
*/
// main function
void inter2_beta_estimation(int*** imgr_bp , int*** img_bp_ref, int*** img_bp_ref2, int** motionVect, double* beta_t, double* beta_t2, const int &imgh, const int &imgw, const int &mbSize){

    int ** img_ref = new2d<int>(imgh,imgw);

    for(int t_lvl=0; t_lvl < PXL; ++t_lvl){

        motionComp(imgr_bp[t_lvl],motionVect,imgh,imgw,mbSize,img_ref);

        inter2_beta_estimation(img_bp_ref[t_lvl],img_ref,img_bp_ref2[t_lvl],beta_t[t_lvl],beta_t2[t_lvl],imgh,imgw);

    }

    delete2d<int>(img_ref);
}

/**
*   @Penlin: [beta_s,beta_t] = intra_inter_beta_estimation(img_bp, img_bp_ref)
*
*   @param imgr_bp, beta_t, beta_s, imgr_bp_ref
*/

void intra_inter_beta_estimation(int*** imgr_bp, int*** img_bp_ref, int** motionVect, double* beta_s, double* beta_t, const int &imgh, const int &imgw, const int &mbSize){

    int ** img_ref = new2d<int>(imgh,imgw);

    for(int t_lvl=0; t_lvl < PXL; ++t_lvl){

        motionComp(img_bp_ref[t_lvl],motionVect,imgh,imgw,mbSize,img_ref);

        intra_inter_beta_estimation(imgr_bp[t_lvl],img_ref,beta_s[t_lvl],beta_t[t_lvl],imgh,imgw);

    }

    delete2d<int>(img_ref);
}

/**
*   @Penlin: [beta_s beta_t beta_t2] = intra_inter2_beta_estimation(img_bp, img_bp_ref, img_bp_ref2)
*
*   @param imgr_bp, imgr_bp_ref, img_bp_ref2
*   @param [output]beta_s, beta_t, beta_t2
*/

void intra_inter2_beta_estimation(int*** imgr_bp, int*** img_bp_ref,int*** img_bp_ref2, int** motionVect, double* beta_s, double* beta_t, double* beta_t2, const int &imgh, const int &imgw, const int &mbSize){

    int ** img_ref = new2d<int>(imgh,imgw);

    for(int t_lvl=0; t_lvl < PXL; ++t_lvl){

        motionComp(imgr_bp[t_lvl],motionVect,imgh,imgw,mbSize,img_ref);

        intra_inter2_beta_estimation(img_bp_ref[t_lvl],img_ref,img_bp_ref2[t_lvl],beta_s[t_lvl],beta_t[t_lvl],beta_t2[t_lvl],imgh,imgw);

    }

    delete2d<int>(img_ref);
}

/**  ======================================================================  **/
/**  ==== estimation methods is above, the mrf_siso methods is below  =====  **/
/**  ======================================================================  **/

/**
*   @Penlin: Le_mrf = mrf_siso_intra(L,beta)
*
*   @param Le_c, beta
*   @param imgh, imgw
*   @param comp 1 or 0
*
*   @output Lu_s        [PXL*lm]
*/

void mrf_siso_intra(double** Le_c, double* beta, const int &imgh, const int &imgw, double** Lu_s, const int &comp){

//    const int lm = imgh*imgw;
    double** p0 = new2d<double>(imgh,imgw);
    int t_lvl, i, j;

    for(t_lvl = 0 ; t_lvl < PXL ; ++t_lvl){

        for(i = 0 ; i < imgh ; ++i)
            for(j = 0 ; j < imgw ; ++j){
                Lu_s[t_lvl][i*imgw+j] = Le_c[t_lvl][i*imgw+j]*comp;
                p0[i][j] = 1/(exp(Le_c[t_lvl][i*imgw+j])+1);
            }

        for(i = 1 ; i < imgh-1 ; ++i )
            for(j = 1 ; j < imgw-1 ; ++j)
                Lu_s[t_lvl][i*imgw+j] += (2*beta[t_lvl]*(2 - p0[i-1][j] - p0[i+1][j] - p0[i][j-1] - p0[i][j+1]));

    }

    delete2d<double>(p0);
}

/**
*   @Penlin: Le_mrf = mrf_siso_inter(L,L_ref,beta)
*
*   @param Le_c, Le_ref, beta
*   @param imgh, imgw
*   @param comp 1 or 0
*
*   @output Lu_s        [PXL*lm]
*/

void mrf_siso_inter(double** Le_c, double** Le_c_ref, double* beta, const int &imgh, const int &imgw, double** Lu_s, const int &comp){

    int i, t_lvl ,lm = imgh*imgw;

    for(t_lvl = 0 ; t_lvl < PXL ; ++t_lvl)
        for(i=0 ; i < lm ; ++i)
            Lu_s[t_lvl][i] = (Le_c[t_lvl][i]*comp + beta[t_lvl]*(1-2/(1+exp(Le_c_ref[t_lvl][i]))));

}


#endif // __MRF_DECODER_UTILS_H
