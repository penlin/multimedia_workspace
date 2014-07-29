#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include "data_types.h"
#include "io_utils.h"
#include "frame.h"

const int width = 352;
const int height = 288;
const int nframe = 300;

int main(int argc,char* argv[]){

    int lm = width*height;
    //char buf_ori[50],buf_eep[50],buf_uep[50],buf_out[50];
    FILE* fori;
    FILE* fdir;
    FILE* feep;
    FILE* fuep;
    FILE* fout;
    unsigned char* buf_ori = (unsigned char*)malloc(sizeof(unsigned char)*lm);
    unsigned char* buf_dir = (unsigned char*)malloc(sizeof(unsigned char)*lm);
    unsigned char* buf_eep = (unsigned char*)malloc(sizeof(unsigned char)*lm);
    unsigned char* buf_uep = (unsigned char*)malloc(sizeof(unsigned char)*lm);
    unsigned char* buffer = (unsigned char*)malloc(sizeof(unsigned char)*lm*4);
    fori = fopen(argv[1],"r+b");
    fdir = fopen(argv[2],"r+b");
    feep = fopen(argv[3],"r+b");
    fuep = fopen(argv[4],"r+b");
    fout = fopen("sequence\\output.yuv","w+b");
    if(fori==NULL || fdir == NULL || feep ==NULL || fuep ==NULL || fout == NULL){
        printf("files cannot open!\n");
        return 0;
    }

    for(int i = 0 ; i < nframe ; ++i){
        if(i%90 == 0){
            fseek(fdir,0,SEEK_SET);
        }
        fread(buf_ori,1,lm,fori);
        fread(buf_dir,1,lm,fdir);
        fread(buf_eep,1,lm,feep);
        fread(buf_uep,1,lm,fuep);
        memset(buffer,0,4*lm);
        for(int j = 0 ; j < height ; ++j){
            for(int k = 0 ; k < width ; ++k){
                /*
                buffer[j*4*width+k] = buf_ori[j*width+k];               // left-up
                buffer[j*4*width+k+width] = buf_dir[j*width+k];         // right-up
                buffer[j*4*width+k+2*width] = buf_eep[j*width+k];      // left-bottom
                buffer[j*4*width+k+3*width] = buf_uep[j*width+k];       // right-bottom
                */
                buffer[j*2*width+k] = buf_ori[j*width+k];               // left-up
                buffer[j*2*width+k+width] = buf_dir[j*width+k];         // right-up
                buffer[(j+height)*2*width+k] = buf_eep[j*width+k];      // left-bottom
                buffer[(j+height)*2*width+k+width] = buf_uep[j*width+k];       // right-bottom
            }
        }
        fwrite(buffer,1,4*lm,fout);
//        fseek(fori,lm/2,SEEK_CUR);
//        fseek(fdir,lm/2,SEEK_CUR);
//        fseek(feep,lm/2,SEEK_CUR);
//        fseek(fuep,lm/2,SEEK_CUR);

    }

    fclose(fori);
    fclose(fdir);
    fclose(feep);
    fclose(fuep);
    fclose(fout);
    free(buf_ori);
    free(buf_dir);
    free(buf_eep);
    free(buf_uep);
    free(buffer);
    return 1;
}


