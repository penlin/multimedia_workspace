#ifndef __IO_UTILS_H
#define __IO_UTILS_H
#include<stdio.h>
#include<stdlib.h>
#include<assert.h>
#include "data_alloc.h"
#include "build_value.h"
#include "utils.h"


/**#   _______          _________ _______          _________
#  (  ___  )|\     /|\__   __/(  ____ )|\     /|\__   __/
#  | (   ) || )   ( |   ) (   | (    )|| )   ( |   ) (
#  | |   | || |   | |   | |   | (____)|| |   | |   | |
#  | |   | || |   | |   | |   |  _____)| |   | |   | |
#  | |   | || |   | |   | |   | (      | |   | |   | |
#  | (___) || (___) |   | |   | )      | (___) |   | |
#  (_______)(_______)   )_(   |/       (_______)   )_(
#  **/

/**
*  @Penlin: write .yuv video with Y value matrix (or U , V)
*  @param: filename
*  @param: imgh
*  @param: imgw
*  @param: n_frame
*  @param: Y       Y value matrix [n_frame*imgh*imgw]
**/

void yuv_write(const char* filename, int*** Y, int*** U, int*** V, const int &n_frame, const int &imgh, const int &imgw){

#if !__OUTPUT_SEQ__
    return;
#endif
    const int lm = imgh*imgw;
    FILE *pFile;
    unsigned char* buffer;

    char file[50];
#if __OUTPUT_TYPE__ == __YUV420__
    sprintf(file,"%s/decoded_yuv420_%s.yuv","sequence",filename);
#else
    sprintf(file,"%s/decoded_y_%s.yuv","sequence",filename);
#endif
    pFile = fopen(file,"w+b");
    assert(pFile!=NULL);

#if __OUTPUT_TYPE__ == __YUV420__
    buffer = (unsigned char*)malloc(sizeof(unsigned char)*(lm*3/2));
    const int offset_u = lm, offset_v = lm*5/4;
#else
    buffer = (unsigned char*)malloc(sizeof(unsigned char)*lm);
#endif
    for(int i = 0 , j = 0 , k = 0; i < n_frame ; ++i){
        printf("write frame #%d\n",i+1);
        for(j = 0 ; j < imgh ; ++j){
            for(k = 0 ; k < imgw ; ++k)
                buffer[k+j*imgw] = (unsigned char) Y[i][j][k];
        }

#if __OUTPUT_TYPE__ == __YUV420__
        if(U==NULL && i > 0){
            fwrite(buffer,1,(lm*3/2),pFile);
            continue;
        }

        for(j = 0 ; j < imgh/2 ; ++j){
            for(k = 0 ; k < imgw/2 ; ++k){
                if(U!=NULL){
                    buffer[k+j*imgw/2+offset_u] = (unsigned char) U[i][j][k];
                    buffer[k+j*imgw/2+offset_v] = (unsigned char) V[i][j][k];
                }else{
                    buffer[k+j*imgw/2+offset_u] = (unsigned char) 0;
                    buffer[k+j*imgw/2+offset_v] = (unsigned char) 0;
                }
            }
        }

        fwrite(buffer,1,(lm*3/2),pFile);
#else
        fwrite(buffer,1,lm,pFile);
#endif
    }

    fclose(pFile);
    free(buffer);
}



void write_exit_info(double*** Ie_c, double*** Ie_c_prev, double*** Ie_s, double*** Ie_s_prev ){

    FILE* f_exit = fopen("output/exit.txt","a+");

    if(f_exit==NULL){
        printf("====  Warning! Error open EXIT Chart file ====\n");
        return;
    }

    fprintf(f_exit,"=============================\n%s%d %d\n",getLocalTimenDate(),__FRAME,Niter);

    for(int frame = 0 , iter = 0 , t_lvl = 0 ; frame < __FRAME ; ++frame){
        for(t_lvl = 0 ; t_lvl < PXL ; ++t_lvl){
            fprintf(f_exit,"frame#%d bp#%d (Ie_c, Ie_s, Ie_c_prev, Ie_s_prev)\n",frame+1,t_lvl+1);
            for(iter = 0 ; iter < Niter ; ++iter)
                fprintf(f_exit,"%lf,",Ie_c[frame][iter][t_lvl]);
            fprintf(f_exit,"\n");
            for(iter = 0 ; iter < Niter ; ++iter)
                fprintf(f_exit,"%lf,",Ie_s[frame][iter][t_lvl]);
            fprintf(f_exit,"\n");
            for(iter = 0 ; iter < Niter ; ++iter)
                fprintf(f_exit,"%lf,",Ie_c_prev[frame][iter][t_lvl]);
            fprintf(f_exit,"\n");
            for(iter = 0 ; iter < Niter ; ++iter)
                fprintf(f_exit,"%lf,",Ie_s_prev[frame][iter][t_lvl]);
            fprintf(f_exit,"\n\n");
        }
    }

    fclose(f_exit);

}


void write_psnr_info(double** PSNR){

    FILE* f_psnr = fopen("output/psnr.txt","a+");

    if(f_psnr==NULL){
        printf("====  Warning! Error open PSNR file ====\n");
        return;
    }

//    fseek(f_psnr,0,SEEK_END);
    fprintf(f_psnr,"=============================\n%s\n",getLocalTimenDate());

    for(int snr = __SNR_S , frame = 0; snr <= __SNR_E ; ++snr){
        fprintf(f_psnr,"SNR=%d\n",snr);
        for(frame = 0 ; frame < __FRAME ; ++frame){
            fprintf(f_psnr,"frame #%3d, PSNR=%lf\n",frame+1,PSNR[snr-__SNR_S][frame]);
        }
    }

    fprintf(f_psnr,"\nTotal Time Comsume:%lf\n\n",getCurrentTime());

    fclose(f_psnr);
}



void write_psnr_info(double* const PSNR){

    FILE* f_psnr = fopen("output/psnr.txt","a+");

    if(f_psnr==NULL){
        printf("====  Warning! Error open PSNR file ====\n");
        return;
    }

    fprintf(f_psnr,"=============================\n%s\n",getLocalTimenDate());
    fprintf(f_psnr,"SNR=%d\n",__SNR);

    for(int i = 0 ; i < __FRAME ; ++i){
        fprintf(f_psnr,"frame#%3d PSNR = %lf\n",i+1,PSNR[i]);
    }

    fprintf(f_psnr,"\nTotal Time Comsume:%lf\n\n",getCurrentTime());

    fclose(f_psnr);
}


void write_power_info(double* const weights, const double &psnr){

    FILE* f_power = fopen("output/power.txt","a+");

    if(f_power==NULL){
        printf("====  Warning! Error open POWER file ====\n");
        return;
    }

    fprintf(f_power,"=============================\n%s\n",getLocalTimenDate());
    fprintf(f_power,"AVG SNR = %d: ",__SNR);

    for(int i = 0 ; i < PXL-1 ; ++i)
        fprintf(f_power,"%.2f,",weights[i]);

    fprintf(f_power,"%.2f:%lf\n",weights[PXL-1],psnr);

//    fprintf(f_power,"\nTotal Time Comsume:%lf\n\n",getCurrentTime());

    fclose(f_power);
}


void write_pc_for_matlab(double* Ia, double* Ie, const int &ln){
    FILE* f = fopen("output/pc.m","w+");

    fprintf(f,"function p = pc()\nIa=[");
    for(int i = 0 ; i < ln ; ++i)
        fprintf(f,"%lf,",Ia[i]);
    fprintf(f,"];\nIe=[");
    for(int i = 0 ; i < ln ; ++i)
        fprintf(f,"%lf,",Ie[i]);
    fprintf(f,"];\np=polyfit(Ia,Ie,3);\nxi = 0:0.01:1;\nyi = polyval(p,xi);\nplot(Ia, Ie,'o','LineWidth',2,'Color',[0 0 0]);\nhold on;\nplot(xi,yi,'LineWidth',1);\naxis([0,1,0,1]);\ngrid on;\nhold off;");
    fclose(f);
}

void write_ps_for_matlab(double** Ia, double** Ie, const int &ln){
    char buffer[50];
    for(int t_lvl = 0 ; t_lvl < PXL ; ++t_lvl){
        sprintf(buffer,"output/ps_bp%d.m",t_lvl+1);
        FILE* f = fopen(buffer,"w+");

        fprintf(f,"function p = ps_bp%d()\nIa=[",t_lvl+1);
        for(int i = 0 ; i < ln ; ++i)
            fprintf(f,"%lf,",Ia[t_lvl][i]);
        fprintf(f,"];\nIe=[");
        for(int i = 0 ; i < ln ; ++i)
            fprintf(f,"%lf,",Ie[t_lvl][i]);
        fprintf(f,"];\np=polyfit(Ia,Ie,3);\nxi = 0:0.01:1;\nyi = polyval(p,xi);\nplot(Ie, Ia,'o','LineWidth',2,'Color',[0 0 0]);\nhold on;\nplot(yi,xi,'LineWidth',1);\naxis([0,1,0,1]);\ngrid on;\nhold off;");
        fclose(f);
    }
}

/**#  _________ _        _______          _________
#  \__   __/( (    /|(  ____ )|\     /|\__   __/
#     ) (   |  \  ( || (    )|| )   ( |   ) (
#     | |   |   \ | || (____)|| |   | |   | |
#     | |   | (\ \) ||  _____)| |   | |   | |
#     | |   | | \   || (      | |   | |   | |
#  ___) (___| )  \  || )      | (___) |   | |
#  \_______/|/    )_)|/       (_______)   )_(
#  **/

/**
*  @Penlin: read .yuv video into Y value matrix
*  @param: filename
*  @param: imgh
*  @param: imgw
*  @param: n_frame
*  @param: Y       Y value matrix [n_frame*imgh*imgw]
**/

int*** yuv_read(const char* filename, const int &imgh, const int &imgw, const int &n_frame){
    FILE *pFile;
    unsigned char * buffer;
//    int*** Y = new3d<int>(imgh,imgw,n_frame);
    int*** Y = new3d<int>(n_frame,imgh,imgw);

    char file[50];
    sprintf(file,"%s/%s","sequence",filename);

    pFile = fopen(file,"r+b");
    assert(pFile!=NULL);
    rewind(pFile);

    buffer = (unsigned char*)malloc(sizeof(unsigned char)*imgh*imgw);
    for(int i = 0 ; i < n_frame ; ++i) {
        printf("read frame #%d\n",i+1);
        fread(buffer,1,imgw*imgh,pFile);
        for(int j = 0 ; j < imgh ; ++j){
            for(int k = 0 ; k < imgw ; ++k)
                Y[i][j][k] = (int) buffer[k+j*imgw];
        }
        fseek(pFile,imgw*imgh/2,SEEK_CUR);
    }

    fclose(pFile);
    free(buffer);
    return Y;
}

void yuv_read(const char* filename, const int &imgh, const int &imgw, const int &n_frame, int*** Y, int*** U, int*** V){
    FILE *pFile;
    unsigned char * buffer;

    char file[50];
    sprintf(file,"%s/%s","sequence",filename);

    pFile = fopen(file,"r+b");
    assert(pFile!=NULL);
    rewind(pFile);

    buffer = (unsigned char*)malloc(sizeof(unsigned char)*imgh*imgw*3/2);
    const int offset_u = imgh*imgw, offset_v = imgh*imgw*5/4;

    for(int i = 0, j = 0, k = 0 ; i < n_frame ; ++i) {
        printf("read frame #%d\n",i+1);
        fread(buffer,1,imgw*imgh*3/2,pFile);
        for(j = 0 ; j < imgh ; ++j){
            for(k = 0 ; k < imgw ; ++k)
                Y[i][j][k] = (int) buffer[k+j*imgw];
        }

        if(U==NULL)
            continue;
        for(j = 0 ; j < imgh/2 ; ++j){
            for(k = 0 ; k < imgw/2 ; ++k){
                U[i][j][k] = (int) buffer[offset_u + k+j*imgw/2];
                V[i][j][k] = (int) buffer[offset_v + k+j*imgw/2];
            }
        }
    }

    fclose(pFile);
    free(buffer);
}

void yuv_random_read(const char* filename, const int &imgh, const int &imgw, const int &n_frame, int*** Y, int*** U, int*** V){
    FILE *pFile;
    unsigned char * buffer;
    // initial _firstFrame at random(0,#Frame-n_frame)
    int _firstFrame = 0 ;
    char file[50];
    sprintf(file,"%s/%s","sequence",filename);

    pFile = fopen(file,"r+b");
    assert(pFile!=NULL);
    rewind(pFile);
    fseek(pFile,0L,SEEK_END);
  //  _firstFrame = 50 ;//random_select_from(0,ftell(pFile)/(imgh*imgw*3/2)-n_frame);
    printf("first frame is %d\n",_firstFrame);
    fseek(pFile,imgh*imgw*3/2*_firstFrame,SEEK_SET);
    buffer = (unsigned char*)malloc(sizeof(unsigned char)*imgh*imgw*3/2);
    const int offset_u = imgh*imgw, offset_v = imgh*imgw*5/4;

    for(int i = 0, j = 0, k = 0 ; i < n_frame ; ++i) {
        printf("read frame #%d\n",i+1);
        fread(buffer,1,imgw*imgh*3/2,pFile);
        for(j = 0 ; j < imgh ; ++j){
            for(k = 0 ; k < imgw ; ++k)
                Y[i][j][k] = (int) buffer[k+j*imgw];
        }

        if(U==NULL)
            continue;
        for(j = 0 ; j < imgh/2 ; ++j){
            for(k = 0 ; k < imgw/2 ; ++k){
                U[i][j][k] = (int) buffer[offset_u + k+j*imgw/2];
                V[i][j][k] = (int) buffer[offset_v + k+j*imgw/2];
            }
        }
    }

    fclose(pFile);
    free(buffer);
}



void deleteY(int *** Y){
    delete3d<int>(Y);
}

#endif // __IO_UTILS_H
