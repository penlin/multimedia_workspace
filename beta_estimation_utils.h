#ifndef __BETA_ESTIMATION_H
#define __BETA_ESTIMATION_H

#define PXL 8


/**  ======================================================================  **/
/**  ==================== for intra beta estimation =======================  **/
/**  ======================================================================  **/

// for intra_beta_estimation
void intra_beta_estimation(int** imgr_bp, double &beta, const int &imgh, const int &imgw){

    const int N = 4;
    int i, j;
    int** hash_key = new2d<int>(imgh-2,imgw-2);
    int** his_idx = new2d<int>(imgh-2,imgw-2);
    int* his = (int*) malloc(sizeof(int)*32);

    for(i=0 ; i < 32; ++i)
        his[i] = 1;

    for(i=0 ; i<imgh-2 ; ++i){
        for(j=0; j<imgw-2 ; ++j){
            hash_key[i][j] = (8*imgr_bp[i][j+1]+4*imgr_bp[i+2][j+1]+2*imgr_bp[i+1][j]+imgr_bp[i][j+2]);
            his_idx[i][j] = (16*imgr_bp[i+1][j+1] + hash_key[i][j]);
            his[his_idx[i][j]]++;
        }
    }

    double aa = 0, ab = 0 , a = 0, b = 0;
    for(i=0; i < imgh-2 ; ++i){
        for(j=0; j<imgw-2 ; ++j){
            a = 2*(imgr_bp[i][j+1]+imgr_bp[i+2][j+1]+imgr_bp[i+1][j]+imgr_bp[i][j+2])-N;
            b = (log(his[hash_key[i][j]+16]) - log(his[hash_key[i][j]]));

            aa+= (a*a);
            ab+= (a*b);
        }
    }

    beta = ab/aa;

    delete2d<int>(hash_key);
    delete2d<int>(his_idx);
    free(his);

}

/**  ======================================================================  **/
/**  ==================== for inter beta estimation =======================  **/
/**  ======================================================================  **/

// for inter_beta_estimation
void inter_beta_estimation(int** imgr_bp , double &beta, int** img_bp_ref, const int &imgh, const int &imgw){

    int sum = 0, i, j;
    int** his_idx = new2d<int>(imgh,imgw);
    int his[4] = {1,1,1,1};

    for(i=0 ; i<imgh; ++i )
        for(j=0; j<imgw ; ++j){
            sum+= img_bp_ref[i][j];
            his_idx[i][j] = (2*imgr_bp[i][j] + img_bp_ref[i][j]);
            if(his_idx[i][j]>=4)
                printf("his_idx out of bound:(%d,%d), imgr_bp = %d, img_bp_ref = %d\n",i,j,imgr_bp[i][j],img_bp_ref[i][j]);
            else
                his[his_idx[i][j]]++;
        }
/*
    for(i=0;i<4;++i)
        printf("his[%d]=%d\n",i,his[i]);
    printf("sum=%d\n",sum);
*/
    beta = ( sum*(log(his[3])-log(his[1])) + (imgh*imgw-sum)*(log(his[0]) - log(his[2])) )/(double)(imgh*imgw);

    delete2d<int>(his_idx);
}

// for inter2_beta_estimation
void inter2_beta_estimation(int** imgr_bp , int** img_bp_ref, int** img_bp_ref2, double &beta_t, double &beta_t2, const int &imgh, const int &imgw){

    int i, j ;
    int** hash_key = new2d<int>(imgh,imgw);
    int** his_idx = new2d<int>(imgh,imgw);
    int his[8] = {1,1,1,1,1,1,1,1};

    for(i=0 ; i < imgh ; ++i){
        for(j=0 ; j<imgw  ; ++j){
            hash_key[i][j] = ( 2*img_bp_ref[i][j] + img_bp_ref2[i][j] );
            his_idx[i][j] = ( 4*imgr_bp[i][j] + hash_key[i][j]);
            his[his_idx[i][j]]++;
        }
    }

    //[beta_s beta_t] = sum((b.')*a)/sum(a.^2);
    double aa = 0, aa2 = 0 , a_a = 0, ab = 0 , ab2 = 0 , b = 0;
    for(i = 0 ; i < imgh ; ++i ){
        for(j=0 ; j <imgw ; ++j){
            b = (log(his[hash_key[i][j]+4])-log(his[hash_key[i][j]]));

            ab+=((2*img_bp_ref[i][j]-1)*b);
            ab2+=((2*img_bp_ref2[i][j]-1)*b);

            aa+=((2*img_bp_ref[i][j]-1)*(2*img_bp_ref[i][j]-1));
            aa2+=((2*img_bp_ref2[i][j]-1)*(2*img_bp_ref2[i][j]-1));

            a_a+=((2*img_bp_ref2[i][j]-1)*(2*img_bp_ref[i][j]-1));
        }
    }

    b = aa*aa2 - a_a*a_a;

    beta_t = (aa2*ab-a_a*ab2)/b ;
    beta_t2 = (aa*ab2-a_a*ab)/b;

    delete2d<int>(hash_key);
    delete2d<int>(his_idx);

}


/**  ======================================================================  **/
/**  ================= for inter intra beta estimation ====================  **/
/**  ======================================================================  **/

// for intra_inter_beta_estimation
void intra_inter_beta_estimation(int** imgr_bp , int** img_bp_ref, double &beta_s, double &beta_t, const int &imgh, const int &imgw){

    const int N = 5;
    int i, j ,range = 64;
    int** hash_key = new2d<int>(imgh-2,imgw-2);
    int** his_idx = new2d<int>(imgh-2,imgw-2);
    int* his = (int*) malloc(sizeof(int)*range);

    for(i=0 ; i< range; ++i)
        his[i] = 1;

    for(i=0 ; i < imgh-2 ; ++i){
        for(j=0 ; j<imgw-2 ; ++j){
            hash_key[i][j] = ( 16*imgr_bp[i][j+1] + 8*imgr_bp[i+2][j+1] + 4*imgr_bp[i+1][j] + 2*imgr_bp[i][j+2] + img_bp_ref[i+1][j+1]);
            his_idx[i][j] = ( 32*imgr_bp[i+1][j+1] + hash_key[i][j]);
            his[his_idx[i][j]]++;
        }
    }
/*
    for(int i = 0 ,cnt = 0; i < range ; ++i)
        printf("his[%d]=%d ... %d\n",i,his[i],(cnt+=his[i]));
*/

    // beta = sum((b.')*a)/sum(a.^2);
    double aa = 0, aa2 = 0 , ab = 0 , ab2 = 0, a = 0 ,b = 0 , a_a = 0;
    for(i = 0 ; i < imgh-2 ; ++i ){
        for(j=0 ; j <imgw-2 ; ++j){
            b = (log(his[hash_key[i][j]+32])-log(his[hash_key[i][j]]));
            a = (imgr_bp[i][j+1] + imgr_bp[i+2][j+1] + imgr_bp[i+1][j] + imgr_bp[i][j+2]);

            ab+=((2*img_bp_ref[i+1][j+1]-1)*b); // a_t
            ab2+=((2*a-(N-1))*b);   // a_s

            aa+=((2*img_bp_ref[i+1][j+1]-1)*(2*img_bp_ref[i+1][j+1]-1));  //a_t
            aa2+=((2*a-(N-1))*(2*a-(N-1))); //a_s

            a_a+=((2*img_bp_ref[i+1][j+1]-1)*(2*a-(N-1)));
        }
    }

    b = aa*aa2 - a_a*a_a;

    beta_t = (aa2*ab-a_a*ab2)/b ;
    beta_s = (aa*ab2-a_a*ab)/b;

    delete2d<int>(hash_key);
    delete2d<int>(his_idx);

    free(his);

}

// for intra_inter2_beta_estimation
void intra_inter2_beta_estimation(int** imgr_bp , int** img_bp_ref, int** img_bp_ref2, double &beta_s, double &beta_t, double &beta_t2, const int &imgh, const int &imgw){

    int i, j ,range = 128;
    int** hash_key = new2d<int>(imgh-2,imgw-2);
    int** his_idx = new2d<int>(imgh-2,imgw-2);
    int* his = (int*) malloc(sizeof(int)*range);

    for(i=0 ; i< range; ++i)
        his[i] = 1;

    for(i=0 ; i < imgh-2 ; ++i){
        for(j=0 ; j<imgw-2 ; ++j){
            hash_key[i][j] = ( 32*imgr_bp[i][j+1] + 16*imgr_bp[i+2][j+1] + 8*imgr_bp[i+1][j] + 4*imgr_bp[i][j+2] + 2*img_bp_ref[i+1][j+1] + img_bp_ref2[i+1][j+1] );
            his_idx[i][j] = ( 64*imgr_bp[i+1][j+1] + hash_key[i][j]);
            his[his_idx[i][j]]++;
        }
    }

    // beta = sum((b.')*a)/sum(a.^2);
    double  ab = 0 , ab2 = 0 , ab3 = 0, a = 0 ,b = 0;
    double** X = new2d<double>(3,3,0);

    for(i = 0 ; i < imgh-2 ; ++i ){
        for(j=0 ; j <imgw-2 ; ++j){
            b = (log(his[hash_key[i][j]+64])-log(his[hash_key[i][j]]));
            a = (imgr_bp[i][j+1] + imgr_bp[i+2][j+1] + imgr_bp[i+1][j] + imgr_bp[i][j+2]);
/*
            ab+=((2*img_bp_ref[i+1][j+1]-1)*b); // a_t
            ab2+=((2*img_bp_ref2[i+1][j+1]-1)*b);   // a_t2
            ab3+=((2*a-4)*b);  // 3a_s

            aa+=((2*img_bp_ref[i+1][j+1]-1)*(2*img_bp_ref[i+1][j+1]-1));  //a_t
            aa2+=((2*img_bp_ref2[i+1][j+1]-1)*(2*img_bp_ref2[i+1][j+1]-1)); //a_t2
            aa3+=((2*a-4)*(2*a-4)); // a_s
*/
            ab+=((2*img_bp_ref[i+1][j+1]-1)*b); // a_t
            ab2+=((2*img_bp_ref2[i+1][j+1]-1)*b);   // a_t2
            ab3+=((2*a-4)*b);  // 3a_s

            X[0][0]+=((2*img_bp_ref[i+1][j+1]-1)*(2*img_bp_ref[i+1][j+1]-1));  //a_t
            X[1][1]+=((2*img_bp_ref2[i+1][j+1]-1)*(2*img_bp_ref2[i+1][j+1]-1)); //a_t2
            X[2][2]+=((2*a-4)*(2*a-4)); // a_s

            X[0][1]+=((2*img_bp_ref[i+1][j+1]-1)*(2*img_bp_ref2[i+1][j+1]-1));  //a_t*a_t2
            X[0][2]+=((2*img_bp_ref[i+1][j+1]-1)*(2*a-4)); //a_t*a_s
            X[1][2]+=((2*img_bp_ref2[i+1][j+1]-1)*(2*a-4)); // a_s*a_t2

        }
    }

    X[1][0] = X[0][1]/=1000; X[2][0] = X[0][2]/=1000 ; X[2][1] = X[1][2]/=1000;
    X[0][0] /= 1000; X[1][1] /= 1000 ; X[2][2] /= 1000;

    double** invX = inv3x3(X);
    beta_t = (invX[0][0]*ab + invX[0][1]*ab2 + invX[0][2]*ab3)/1000;
    beta_t2 = (invX[1][0]*ab + invX[1][1]*ab2 + invX[1][2]*ab3)/1000;
    beta_s = (invX[2][0]*ab + invX[2][1]*ab2 + invX[2][2]*ab3)/1000;
//    beta_t = ab/aa ;
//    beta_t2 = ab2/aa2;
//    beta_s = ab3/aa3;

    delete2d<int>(hash_key);
    delete2d<int>(his_idx);
    delete2d<double>(X);
    delete2d<double>(invX);
    free(his);

}

#endif // __BETA_ESTIMATION_H
