#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include "data_types.h"
#include "io_utils.h"

#define PXL  8
const int width = 352;
const int height = 288;
int skip = 0;

int main(int argc,char* argv[]){

    int lm = width*height;
    int mask = 1;
    //char buf_ori[50],buf_eep[50],buf_uep[50],buf_out[50];
    char buf[50];
    FILE* fori;
    FILE* fout;
    unsigned char* buf_ori = (unsigned char*)malloc(sizeof(unsigned char)*lm);
    unsigned char* buffer = (unsigned char*)malloc(sizeof(unsigned char)*lm);
    fori = fopen(argv[1],"r+b");
    fout =  fopen("sequence\\output.yuv","w+b");
    if(fori==NULL ){
        printf("files cannot open!\n");
        return 0;
    }

    if(argc > 2){
        skip = atoi(argv[2]);
    }

    fseek(fori,skip,SEEK_SET);
    fread(buf_ori,1,lm,fori);
    fclose(fori);
    for(int i = 0 ; i < PXL ; ++i){
        for(int j = 0 ; j < lm ; ++j){
            buffer[j] = (buf_ori[j] & mask)*255;
        }
        fwrite(buffer,1,lm,fout);
        mask = mask << 1;
    }
    fclose(fout);

    free(buf_ori);
    free(buffer);
    return 1;
}


