#ifndef __MOTION_ESTIMATION_UTILS_H
#define __MOTION_ESTIMATION_UTILS_H
#include <stdio.h>
#include <stdlib.h>
#include <limits>
#include <float.h>

#define PXL 8

typedef struct Point {
    int X;
    int Y;
}Point;

const int ORDERS[8] = {128 , 64, 32, 16, 8, 4, 2, 1};

static int height;
static int width;
static int mbSize;
static int p;

template<class T>
static void cal_soft_map(T* const a, T* const b, double &err){
    int i, j;

    err = 0;
    for(i=0 ; i <PXL ; ++i){
        err+= ORDERS[i]*((a[i]>=b[i])?(a[i]-b[i]):(b[i]-a[i]));
    }
}

template<class T>
static void cal_soft_likelihood(T* const a, T* const b, double &err){
    int i, j;
    double ea, eb;
    err = 0;
    for(i=0 ; i <PXL ; ++i){
        //ea = exp(a[i]);
        //eb = exp(b[i]);
        err+= ORDERS[i]*(log((a[i]+b[i])/(a[i]*b[i]+1)));
    }
}
template<class T>
static void cal_soft_mse(T* const a, T* const b, double &err){
    int i, j;

    err = 0;
    for(i=0 ; i <PXL ; ++i){
        for(j=0 ; j < PXL ; ++j)
            if(i!=j)
                err += (a[i]*a[j] - 2*a[i]*b[j] + b[i]*b[j])*ORDERS[i]*ORDERS[j];
            else
                err += (a[i] - 2*a[i]*b[j] + b[j])*ORDERS[i]*ORDERS[j];
    }
}

template<class T>
static void cal_hard_mad(T* const a, T* const b, double &err){
    int i, p1 = 0, p2=0;

    err = 0;
    for(i=0 ; i <PXL ; ++i){
        p1 += ((int)(a[i]+0.5)/1)*ORDERS[i];
        p2 += ((int)(b[i]+0.5)/1)*ORDERS[i];
    }
    if(p1>p2)
        err += (p1 - p2);
    else
        err += (p2 - p1);

}

template<class T>
static void cal_hard_mse(T* const a, T* const b, double &err){
    int i, p1=0, p2=0;

    err = 0;
    for(i=0 ; i <PXL ; ++i){
        p1 += ((int)(a[i]+0.5)/1)*ORDERS[i];
        p2 += ((int)(b[i]+0.5)/1)*ORDERS[i];
    }

    err += (p1-p2)*(p1-p2);

}

/*
*   @Penlin: calculate cost function for window
*
*   @param start,delta : define the window domain in imgP and imgI
*   @param err: [output] the cost for this window
*/
template<class T>
void costFunc(T*** const imgP, T*** const imgI, Point const &start, Point const &delta, double &err){
    int i ,j, sX, sY, eY, eX, dx, dy;
    double tmp = 0;
    err = 0;
    sY = start.Y;
    sX = start.X;

    dx = delta.X;
    dy = delta.Y;

    eY = sY + mbSize ;
    eX = sX + mbSize ;

    for(i= sY ; i < eY ; ++i){
        for(j = sX ; j < eX ; ++j){
            cal_soft_likelihood<T>(imgP[i][j],imgI[i+dy][j+dx],tmp);
            err += tmp;
        }
    }

}

/**
*   @Penlin: Main function for Motion Estimation
*
*   @param ImgP  soft imgr_bp matrix [imgh*imgw*PXL]
*   @param ImgI  soft imgr_bp matrix [imgh*imgw*PXL]
*   @param h,w,mb_size,me_range
*
*   @param mv   output value         [2*mb_count]
*/

// soft version
template<class T>
void motionEstES(T*** const imgP, T*** const imgI,const int &h, const int &w , const int &mb_size, const int &me_range, int** mv){
    int i , j, m, n, curVer, curHor, cnt=0 , nMV;
    Point point, delta;
    double cost, minCost = DBL_MAX;

    height = h;
    width = w;
    mbSize = mb_size;
    p = me_range;
    nMV = height*width/(mbSize*mbSize);

    for(i=0 ; i < height - mbSize + 1 ; i+=mbSize ){
        for(j=0; j < width - mbSize + 1 ; j+=mbSize){

            for(m = -p ; m <= p ; ++m){
                for(n = -p ; n <= p ; ++n){

                    curVer = i + m;
                    curHor = j + n;

                    if((curVer + mbSize ) > height || curVer < 0 || (curHor + mbSize ) > width || curHor < 0 )
                        continue;
                    else{

                        point.Y = i;
                        point.X = j;

                        delta.Y = m;
                        delta.X = n;
                        costFunc<T>(imgP,imgI,point,delta,cost);

                        if( cost < minCost ){
                            minCost = cost;
                            mv[0][cnt] = m;
                            mv[1][cnt] = n;
                        }
                    }

                }
            }

            cnt++;
            minCost = DBL_MAX;
        }
    }

}

// hard version motion estimation
void costFunc(Pixel** const imgP, Pixel** const imgI, Point const &start, Point const &delta, double &err){
    int i ,j, sX, sY, eY, eX, dx, dy;
    double tmp = 0;
    err = 0;
    sY = start.Y;
    sX = start.X;

    dx = delta.X;
    dy = delta.Y;

    eY = sY + mbSize ;
    eX = sX + mbSize ;

    for(i= sY ; i < eY ; ++i){
        for(j = sX ; j < eX ; ++j){
            tmp = ((int16)imgP[i][j] - (int16)imgI[i+dy][j+dx]);
            tmp = (tmp>0?tmp:-tmp);
            err += tmp;
        }
    }

}

// hard version
void motionEstES(Pixel** const imgP, Pixel** const imgI,const int &h, const int &w , const int &mb_size, const int &me_range, int** mv){
    int i , j, m, n, curVer, curHor, cnt=0 ;
    Point point, delta;
    double cost, minCost = -1;

    height = h;
    width = w;
    mbSize = mb_size;
    p = me_range;
//    nMV = height*width/(mbSize*mbSize);

    for(i=0 ; i < height - mbSize + 1 ; i+=mbSize ){
        for(j=0; j < width - mbSize + 1 ; j+=mbSize){

            for(m = -p ; m <= p ; ++m){
                for(n = -p ; n <= p ; ++n){

                    curVer = i + m;
                    curHor = j + n;

                    if((curVer + mbSize ) > height || curVer < 0 || (curHor + mbSize ) > width || curHor < 0 )
                        continue;
                    else{

                        point.Y = i;
                        point.X = j;

                        delta.Y = m;
                        delta.X = n;
                        costFunc(imgP,imgI,point,delta,cost);

                        if(minCost < 0 || cost < minCost ){
                            minCost = cost;
                            mv[0][cnt] = m;
                            mv[1][cnt] = n;
                        }
                    }

                }
            }

            cnt++;
            minCost = -1;
        }
    }

}



/**
*   @Penlin: imgComp = motionComp(imgI, motionVect, mbSize)
*   @param imgI         reference frame
*   @param motionVect   output of ME
*   @param imgh, imgw, mbSize
*
*   @param imgComp      output compansented frame
*/
void motionComp(int8** imgI, int** motionVect, const int &imgh, const int &imgw, const int &mbSize, int8** imgComp){

    int i, j, mbCount = 0,ii,jj,refBlkVer ,refBlkHor ;

//    printf("imgh=%d,imgw=%d,mbSize=%d\n",imgh,imgw,mbSize);

    for(i = 0 ; i <= imgh-mbSize ; i+=mbSize){
        for(j = 0 ; j <= imgw-mbSize ; j+=mbSize){

//            printf("i=%d,j=%d,mv[0]=%d,mv[1]=%d (count=%d)\n",i,j,motionVect[0][mbCount],motionVect[1][mbCount],mbCount);
            refBlkVer = i + motionVect[0][mbCount];
            refBlkHor = j + motionVect[1][mbCount];

            for(ii=0 ; ii < mbSize ; ++ii)
                for(jj=0 ; jj < mbSize; ++jj)
                    imgComp[i+ii][j+jj] = imgI[refBlkVer+ii][refBlkHor+jj];

            mbCount++;
        }
    }

}

/**
*   @Penlin: LeComp = motionComp(Le_c, motionVect, mbSize)
*   @param Le_c         reference Le
*   @param motionVect   output of ME
*   @param imgh, imgw, mbSize
*
*   @param LeComp      output compansented Le
*/
void motionComp(double** Le, int** motionVect, const int &imgh, const int &imgw, const int &mbSize, double** LeComp){

    int i, j, mbCount = 0,ii,jj,refBlkVer ,refBlkHor ,t_lvl;

//    printf("imgh=%d,imgw=%d,mbSize=%d\n",imgh,imgw,mbSize);

    for(t_lvl=0 ; t_lvl<PXL ; ++t_lvl){

        mbCount = 0;
        for(i = 0 ; i <= imgh-mbSize ; i+=mbSize){
            for(j = 0 ; j <= imgw-mbSize ; j+=mbSize){

//                printf("i=%d,j=%d,mv[0]=%d,mv[1]=%d (count=%d)\n",i,j,motionVect[0][mbCount],motionVect[1][mbCount],mbCount);
                refBlkVer = i + motionVect[0][mbCount];
                refBlkHor = j + motionVect[1][mbCount];

                for(ii=0 ; ii < mbSize ; ++ii)
                    for(jj=0 ; jj < mbSize; ++jj)
                        LeComp[t_lvl][(i+ii)*imgw+j+jj] = Le[t_lvl][(refBlkVer+ii)*imgw+refBlkHor+jj];

                mbCount++;
            }
        }
    }

}


#endif // __MOTION_ESTIMATION_UTILS_H
