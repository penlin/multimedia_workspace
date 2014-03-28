#include "build_value.h"
#include "io_utils.h"
#include "video_encode.h"
#include "data_alloc.h"
#include "utils.h"
#include "uep_predict_utils.h"


const double DERIVE_SNR_N1[PXL] = {3.015249, 2.409902, 1.813491, 0.728476, 0.032169, 0.000001, 0.000001, 0.000001 };
const double DERIVE_SNR_0[PXL] = {2.6493 , 2.233, 1.7531, 1.2123,  0.1237 , 0.01, 0.01, 0.0086 }; // 29 dB
const double DERIVE_SNR_1[PXL] = {2.395 , 1.9641, 1.5781, 1.2324,  0.76924 , 0.045221, 0.0079442, 0.0079442 }; // 29 dB
const double DERIVE_SNR_2[PXL] = {2.0161 , 1.8493, 1.5252, 1.221,  0.9304 , 0.43326, 0.018426, 0.0063103 }; // 29 dB
const double TEST[PXL] = {1.977668, 1.538618, 0.943222, 0.049886, 0.01, 0.01, 0.01, 0.01};
const double TEST_TEST[PXL]= {3.476999, 2.704871, 1.657564, 0.087482, 0.018022, 0.018022, 0.018022, 0.018022};

const double fr[PXL] = {0.0026   , 0.0097 ,  0.0484  ,  0.2058   , 0.2298  ,  0.2298  ,  0.2298  ,  0.2298 };
void bit_analysis(const char* filename, int*** Y, const int &imgh, const int &imgw, const int &n_frame, int** bit_cnt){

    int i,j,k, y, f, t_lvl,n,pixel, mask = 1;
    const int offset_u = imgh*imgw, offset_v = imgh*imgw*5/4;
    int** imgO;
    int*** img_bp = new3d<int>(PXL,imgh,imgw);
    double correlation = 0.0, est = 0.0;


    FILE *pFile;
    unsigned char * buffer;

    pFile = fopen(filename,"r+b");
    assert(pFile!=NULL);
    rewind(pFile);

    buffer = (unsigned char*)malloc(sizeof(unsigned char)*imgh*imgw*3/2);

    for(i = 0, j = 0, k = 0 ; i < n_frame ; ++i) {
        printf("read frame #%d\n",i+1);
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

        for(t_lvl=0 ; t_lvl<PXL ; ++t_lvl){
            for(j=0 ; j < imgh ; ++j){
                for(k=0 ; k< imgw ; ++k){
                    bit_cnt[i][t_lvl]+=img_bp[t_lvl][j][k];
                }
            }
        }


        correlation= 0.0;
        for(j=0 ; j < imgh ; ++j)
            for(k=0 ; k< imgw ; ++k)
                for(t_lvl=0 ; t_lvl<PXL ; ++t_lvl)
                    for(n=0; n < PXL ; ++n){
                        if(n==t_lvl){
                            est = fr[t_lvl];
                        }else{
                            est = fr[t_lvl]*fr[n]*(1-2*img_bp[t_lvl][j][k])*(1-2*img_bp[n][j][k]);
                        }
                        correlation+=(1<<(PXL-t_lvl-1))*(1<<(PXL-n-1))*est;
                        //correlation+=img_bp[t_lvl][j][k]*img_bp[n][j][k];
                    }

        correlation = correlation/imgw/imgh;
        printf("#%d: E[bl,bn] = %lf\n",i,correlation);

    }


    fclose(pFile);
    free(buffer);
    delete3d<int>(img_bp);

}

void write_analysis_info(int** bit_cnt, const int &imgh, const int &imgw, const int &n_frame){

    FILE* f_ptr = fopen("bit_analysis.txt","a+");

    if(f_ptr==NULL){
        printf("====  Warning! Error open PSNR file ====\n");
        return;
    }

    fprintf(f_ptr,"=============================\n%s\n",getLocalTimenDate());

    for(int i = 0 ; i < n_frame ; ++i){
        for(int j = 0 ; j < PXL ; ++j){
            fprintf(f_ptr,"%d ",bit_cnt[i][j]);
        }
        fprintf(f_ptr, ";");
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
    int** bit_cnt = new2d<int>(f,PXL,0);
    double* weights = (double*)malloc(sizeof(double)*PXL);
    for(int i = 0 ; i < PXL ; ++i)
        weights[i] = 1;

    bit_analysis("../sequence/akiyo_cif.yuv",Y,h,w,f,bit_cnt);
    write_analysis_info(bit_cnt,h,w,f);



    // free memory
    deleteY(Y);
    delete2d<int>(bit_cnt);
    free(weights);
}
