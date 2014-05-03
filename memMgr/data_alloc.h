#ifndef __DATA_ALLOC_H
#define __DATA_ALLOC_H
#include <string.h>
#inlcude "built_value.h"

#if __MEM_MGR__
#include "memMgr.h"
#endif

#ifndef MEM_MGR_H
#define MALLOC(T,x)  (T*)malloc(sizeof(T)*(x))
#define DELETE(x)    free((x))
#else
#define SPLIT         0
#define MALLOC(T,x)  (T*)_memMgr->allocArr(sizeof(T),(x))
#define DELETE(x)    _memMgr->free((x))

static MemMgr* _memMgr;

void initMemMgr(int _blockSize = 65536){_memMgr = new MemMgr(_blockSize);}

void freeMemMgr(){delete _memMgr;}

void printMem(){_memMgr->print();}
#endif

template<class T>
T** new2d(const int &r, const int &c){
#if SPLIT
    T **ret = MALLOC(T*,r);//(T**)malloc(sizeof(T*)*r);
    ret[0] = MALLOC(T,r*c);//(T*)malloc(sizeof(T)*r*c);
    for(int i = 1 ; i < r ; ++i)
        ret[i] = ret[i-1] + c;
#else
    T** ret = (T**)_memMgr->alloc2DMat(sizeof(T),r,c);
    for(int i = 1 ; i < r ; ++i)
        ret[i] = ret[i-1] + c;
#endif
    return ret;
}

template<class T>
T** new2d(const int &r, const int &c , const T &value){
    T** ret = new2d<T>(r,c);
    for(int i = 0 ; i < r*c ; ++i)
        ret[0][i] = value;
    return ret;
}

template<class T>
void delete2d(T** a) {
#if SPLIT
    DELETE(a[0]);
#endif
    DELETE(a);
}


template<class T>
T*** new3d(const int &r, const int &c, const int &d) {
#if SPLIT
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
#else
    T*** ret = (T***)_memMgr->alloc3DMat(sizeof(T),r,c,d);
    if(*(size_t *)ret!=(size_t)&ret[r]){
        ret[0] = (T**)&ret[r];
        ret[0][0] = (T*)&ret[r*(c+1)];
        for(int i = 1 ; i < c ;++i)
            ret[0][i] = ret[0][i-1] + d;
        for(int i = 1,j=0 ,k=c*d; i < r ; ++i){
            ret[i] = ret[i-1]+c;
            for(j=0;j<c;++j){
                ret[i][j] = ret[i-1][j] + k;
            }
        }
    }
#endif
    return ret;
}

template<class T>
T*** new3d(const int &r, const int &c, const int &d, const T &value) {
    T*** ret = new3d<T>(r,c,d);
    for(int i=0; i < r*c*d ; ++i)
        ret[0][0][i] = value;
    return ret;
}

template<class T>
void delete3d(T*** a) {
#if SPLIT
    DELETE(a[0][0]);
    DELETE(a[0]);
#endif
    DELETE(a);
}

#endif
