#include "build_value.h"
#include "io_utils.h"
#include "video_encode.h"
#include "data_alloc.h"
#include "utils.h"
#include "uep_predict_utils.h"
#include "mrf_decoder_utils.h"


void ber_analysis(const char* filename, int*** Y, const int &imgh, const int &imgw, const int &n_frame, double** beta){

    int i,j,k, y, f, t_lvl,n,pixel, mask = 1;
    const int offset_u = imgh*imgw, offset_v = imgh*imgw*5/4;
    int** imgO;
    int*** img_bp = new3d<int>(PXL,imgh,imgw);
    double correlation = 0.0, est = 0.0;
//    double* beta = (double*)malloc(sizeof(double)*PXL);

    FILE *pFile;
    unsigned char * buffer;

    pFile = fopen(filename,"r+b");
    assert(pFile!=NULL);
    rewind(pFile);

    buffer = (unsigned char*)malloc(sizeof(unsigned char)*imgh*imgw*3/2);

    for(i = 0, j = 0, k = 0 ; i < n_frame ; ++i) {
        //printf("read frame #%d\n",i+1);
        fread(buffer,1,imgw*imgh*3/2,pFile);
        for(j = 0 ; j < imgh ; ++j){
            for(k = 0 ; k < imgw ; ++k){
                Y[i][j][k] = (int) buffer[k+j*imgw];
                pixel = Y[i][j][k];
                for(t_lvl=PXL-1 ; t_lvl >= 0 ; --t_lvl){
                    img_bp[t_lvl][j][k] = (int)((pixel & mask)?1:0);
                    pixel = ( pixel >> 1 );
                }
            }
        }

        intra_beta_estimation(img_bp,beta[i],imgh,imgw);
        for(t_lvl=0;t_lvl<PXL;++t_lvl){
            printf("%lf,",beta[i][t_lvl]);
        }
        printf("\n");

    }


    fclose(pFile);
    free(buffer);
//    free(beta);
    delete3d<int>(img_bp);

}

void write_analysis_info(double**beta, const int &imgh, const int &imgw, const int &n_frame){

    FILE* f_ptr = fopen("beta_analysis.txt","a+");

    if(f_ptr==NULL){
        printf("====  Warning! Error open PSNR file ====\n");
        return;
    }

    fprintf(f_ptr,"=============================\n%s\n",getLocalTimenDate());

    for(int i = 0 ; i < n_frame ; ++i){
        for(int j = 0 ; j < PXL ; ++j){
            fprintf(f_ptr,"%lf, ",beta[i][j]);
        }
        fprintf(f_ptr, ";\n");
    }
    printf("\n\n");

    fclose(f_ptr);
}


int main(int argc,char* argv[]){

    startRandom();

    const int h = __HEIGHT, w = __WIDTH, f = __FRAME ;
    const int lm = h*w;
    double snr = __SNR;
    int*** Y = new3d<int>(f,h,w);
    double** beta = new2d<double>(f,PXL,0);
//    double* weights = (double*)malloc(sizeof(double)*PXL);
//    for(int i = 0 ; i < PXL ; ++i)
//        weights[i] = 1;

    ber_analysis("../sequence/stefan_cif.yuv",Y,h,w,f,beta);
    write_analysis_info(beta,h,w,f);



    // free memory
    deleteY(Y);
    delete2d<double>(beta);
//    free(weights);
}
