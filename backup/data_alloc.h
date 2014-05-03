#ifndef __DATA_ALLOC_H
#define __DATA_ALLOC_H
#include <string.h>

#if __MEM_MGR__
#include "memMgr.h"
#endif

#ifndef MEM_MGR_H
#define MALLOC(T,x)  (T*)malloc(sizeof(T)*(x))
#define DELETE(x)    free((x))
#else
#define MALLOC(T,x)  (T*)_memMgr->allocArr(sizeof(T),(x))
#define DELETE(x)    _memMgr->free((x))

static MemMgr* _memMgr;

void initMemMgr(int _blockSize = 65536){_memMgr = new MemMgr(_blockSize);}

void freeMemMgr(){delete _memMgr;}

void printMem(){_memMgr->print();}
#endif


template<class T>
T** new2d(const int &r, const int &c){
    T **ret = MALLOC(T*,r);//(T**)malloc(sizeof(T*)*r);
    ret[0] = MALLOC(T,r*c);//(T*)malloc(sizeof(T)*r*c);
    for(int i = 1 ; i < r ; ++i)
        ret[i] = ret[i-1] + c;
    return ret;
}

template<class T>
T** new2d(const int &r, const int &c , const T &value){
    int i = 0;
    T **ret = MALLOC(T*,r);//(T**)malloc(sizeof(T*)*r);
    ret[0] = MALLOC(T,r*c);//(T*)malloc(sizeof(T)*r*c);
    for(i = 1 ; i < r ; ++i)
        ret[i] = ret[i-1] + c;
//    for(i = 0 ; i < r*c ; ++i)
//        ret[0][i] = value;
    memset(ret[0],value,r*c);
    return ret;
}

template<class T>
void delete2d(T** a) {
//    free(a[0]);
//    free(a);
    DELETE(a[0]);
    DELETE(a);
}


template<class T>
T*** new3d(const int &r, const int &c, const int &d) {
    T*** ret = MALLOC(T**,r);//(T***)malloc(sizeof(T**)*r);
    ret[0] = MALLOC(T*,r*c);//(T**)malloc(sizeof(T*)*r*c);
    ret[0][0] =  MALLOC(T,r*c*d);//(T*)malloc(sizeof(T)*r*c*d);

    T* tmp = ret[0][0];

    for(int i = 0 ; i < r ; ++i ){
        if(i!=0)
            ret[i] = ret[i-1] + c;
        for(int j = 0 ; j < c ; ++j){
            if(!(i==0 && j==0)) {
                ret[i][j] = tmp + d;
                tmp = ret[i][j];
            }
        }
    }
    return ret;
}

template<class T>
T*** new3d(const int &r, const int &c, const int &d, const T &value) {
    int i, j;
    T*** ret = MALLOC(T**,r);//(T***)malloc(sizeof(T**)*r);
    ret[0] = MALLOC(T*,r*c);//(T**)malloc(sizeof(T*)*r*c);
    ret[0][0] =  MALLOC(T,r*c*d);//(T*)malloc(sizeof(T)*r*c*d);

    T* tmp = ret[0][0];

    for(i = 0 ; i < r ; ++i ){
        if(i!=0)
            ret[i] = ret[i-1] + c;
        for(j = 0 ; j < c ; ++j){
            if(!(i==0 && j==0)) {
                ret[i][j] = tmp + d;
                tmp = ret[i][j];
            }
        }
    }
//    for(i=0; i < r*c*d ; ++i)
//        ret[0][0][i] = value;
    memset(ret[0][0],value,r*c*d);
    return ret;
}

template<class T>
void delete3d(T*** a) {
//    free(a[0][0]);
//    free(a[0]);
//    free(a);
    DELETE(a[0][0]);
    DELETE(a[0]);
    DELETE(a);
}




#endif
