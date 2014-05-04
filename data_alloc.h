#ifndef __DATA_ALLOC_H
#define __DATA_ALLOC_H
#include <string.h>
#include "build_value.h"

#if __MEM_MGR__
#include "memMgr.h"
#endif

#define SPLIT         0

static int repeatUsage2 = 0;
static int repeatUsage3 = 0;
static int num2D = 0;
static int num3D = 0;

#ifndef MEM_MGR_H
#define MALLOC(T,x)  (T*)malloc(sizeof(T)*(x))
#define DELETE(x)    free((x))
#else

#define MALLOC(T,x)  (T*)_memMgr->allocArr(sizeof(T),(x))
#define DELETE(x)    _memMgr->free((x))

static MemMgr* _memMgr;

void initMemMgr(int _blockSize = 65536){_memMgr = new MemMgr(_blockSize);}

void freeMemMgr(){delete _memMgr;}

void printRepeat(){
    printf("Totaly allocate 2D Matrix : %d \n2D Matrix Re-use : %d\n",num2D,repeatUsage2);
    printf("Totaly allocate 3D Matrix : %d \n3D Matrix Re-use : %d\n",num3D,repeatUsage3);
}

void printMem(){
    _memMgr->print();
//    printRepeat();
}
#endif

template<class T>
T** new2d(const int &r, const int &c){
//    num2D++;
#if __MEM_MGR__ && !SPLIT
    T** ret = (T**)_memMgr->alloc2DMat(sizeof(T),r,c);
    if(*(size_t *)ret!=(size_t)&ret[r] ){
//        repeatUsage2++;
        ret[0] = (T*)&ret[r];
        for(int i = 1 ; i < r ; ++i)
            ret[i] = ret[i-1] + c;
    }
#else
    T **ret = MALLOC(T*,r);
    ret[0] = MALLOC(T,r*c);
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
#if !( __MEM_MGR__ && !SPLIT)
    DELETE(a[0]);
#endif
    DELETE(a);
}


template<class T>
T*** new3d(const int &r, const int &c, const int &d) {
//    num3D++;
#if __MEM_MGR__ && !SPLIT
    T*** ret = (T***)_memMgr->alloc3DMat(sizeof(T),r,c,d);
    if(*(size_t *)ret!=(size_t)&ret[r]){
//        repeatUsage3++;
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
#else
    T*** ret = MALLOC(T**,r);
    ret[0] = MALLOC(T*,r*c);
    ret[0][0] =  MALLOC(T,r*c*d);

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
#if !( __MEM_MGR__ && !SPLIT)
    DELETE(a[0][0]);
    DELETE(a[0]);
#endif
    DELETE(a);
}


#endif
