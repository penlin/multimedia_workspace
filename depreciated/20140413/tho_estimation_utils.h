#ifndef __THO_ESTIMATION_H
#define __THO_ESTIMATION_H

#define PXL 8

/**  ======================================================================  **/
/**  ==================== for inter beta estimation =======================  **/
/**  ======================================================================  **/

// for inter_beta_estimation
void inter_tho_estimation(int** imgr_bp , double &beta, int** img_bp_ref, const int &imgh, const int &imgw){

    int sum = 0, i, j;
//    int** his_idx = new2d<int>(imgh,imgw);
//    int his[4] = {1,1,1,1};

    for(i=0 ; i<imgh; ++i )
        for(j=0; j<imgw ; ++j)
            sum+=(imgr_bp[i][j]==img_bp_ref[i][j]);

    beta = (double)sum/imgw/imgw;

//    printf("beta=%f\n",beta);
/*
    for(i=0;i<4;++i)
        printf("his[%d]=%d\n",i,his[i]);
    printf("sum=%d\n",sum);
*/
//    beta = ( sum*(log(his[3])-log(his[1])) + (imgh*imgw-sum)*(log(his[0]) - log(his[2])) )/(double)(imgh*imgw);

//    delete2d<int>(his_idx);
}


/**
*   @Penlin:  beta = inter_beta_estimation(img_bp, img_bp_ref)
*
*   @param imgr_bp, beta, imgr_bp_ref , motionVect
*/
// main function
void inter_tho_estimation(int*** imgr_bp , int*** img_bp_ref, int** motionVect, double* beta, const int &imgh, const int &imgw, const int &mbSize){

//    int ** img_ref = new2d<int>(imgh,imgw);

    for(int t_lvl = 0 ; t_lvl < PXL ; ++t_lvl){

//        motionComp(img_bp_ref[t_lvl],motionVect,imgh,imgw,mbSize,img_ref);

        inter_tho_estimation(imgr_bp[t_lvl],beta[t_lvl],img_bp_ref[t_lvl],imgh,imgw);

    }

//    delete2d<int>(img_ref);
}


#endif // __THO_ESTIMATION_H
