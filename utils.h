#ifndef __UTILS_H
#define __UTILS_H
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <time.h>

#include "data_alloc.h"
#include "build_value.h"

static time_t tBuffer;
static char* LOCALTIME;

void startRandom(){
    tBuffer = time(NULL);
    LOCALTIME = ctime(&tBuffer);
#if __RANDOM__
    srand(tBuffer);
#else
    srand(1);
#endif

}


/*
*  @Penlin: adopt modified random shuffle algorithm
*  @param: low    low bound [include]
*  @param: high   high bound [include]
*
*   To shuffle an array a of n elements (indices 0 .. n-1):
*   for i from n - 1 downto 1 do
*           j = random integer with 0 <= j <= i
*           exchange a[j] and a[i]
*/
void random_sequence(const int &low, const int &high, int * order){
    int length = high - low + 1, j, tmp;

    for(int i = 0 ; i < length ; ++i)
        order[i] = low + i;

#if __INTERLEAVE__
//    length = 0;
    for(int i = length-1 ; i >= 1 ; --i){
        j = rand()%(i+1);
        tmp = order[j];
        order[j] = order[i];
        order[i] = tmp;
    }
#endif

}



/*  @Penlin : Standard normal distribution noise generator
*/

double gaussian_noise(){
#if __NOISE__
    static double V1, V2, S;
    static int phase = 0 ;
    double X;

    if(phase == 0 ){
        do{
            double U1 = (double)rand() / RAND_MAX;
            double U2 = (double)rand() / RAND_MAX;

            V1 = 2*U1 -1 ;
            V2 = 2*U2 -1 ;
            S = V1*V1 + V2*V2;
        } while( S >= 1 || S == 0);
        X = V1 * sqrt(-2 * log(S) / S);
    }else
        X = V2 * sqrt(-2 * log(S) / S);

    phase = 1 - phase;
    return X;
#else
    return 0;
#endif
}


int random_select_from(const int &from, const int &to ){
    return rand()%(to-from) + from;
}

template<class T>
void copyMatrix(T** a, T**b , const int &lm){
    for(int i = 0 ; i < lm ; ++i)
        a[0][i] = b[0][i];
}


/*
    n1 = sizeof(a) / sizeof(*a);
    n2 = sizeof(b) / sizeof(*b);

    asssert(n1 > end_a && start_a >= 0 && (end_a-start_a) <= n1
            && n2 > end_b && start_b >= 0 && (end_b-start_b) <= n2
            && (end_b-start_a) == (end_a - start_a));
*/
//int InnerProduct(const int *a,const int *b,const int start_a ,const int end_a ,const int start_b ,const int end_b) {
//    int res = 0;
//    int  i , j;
//
//    for(i=start_a, j=start_b;i<=end_a && j<=end_b ;++i,++j)
//        res += a[i]*b[j];
//
//    return res;
//}

int8 InnerProduct(const int *a,const int8 *b,const int &start_a ,const int &end_a ,const int &start_b ,const int &end_b) {
    int8 res = 0;
    int  i , j;

    for(i=start_a, j=start_b;i<=end_a && j<=end_b ;++i,++j)
        res += a[i]*b[j];

    return res;
}

/**
*   @Penlin: print the 2D matrix from i: [i_s,i_e) and j: [j_s,j_e)
*/

template<class T>
void print2DMatrix(T** mat, const int &i_s, const int &i_e, const int &j_s, const int &j_e){
    for(int i = i_s, j = j_s; i < i_e ; ++i){
        for(j=j_s ; j < j_e ; ++j){
            switch(sizeof(T)){
            case 1: // char
                printf("%c ",mat[i][j]);
                break;
            case 4:
                printf("%d ",mat[i][j]);
                break;
            case 8:
                printf("%lf ",mat[i][j]);
                break;
            }
        }
        printf("\n");
    }
}

/**
*   @Penlin: get the current time for test the effiency
*/

double getCurrentTime(){
//    time_t tBuffer;
    return (float)clock()/CLOCKS_PER_SEC;
}

char* getLocalTimenDate(){
    return LOCALTIME;
}


/**
*   @Penlin: get the inverse of a 3x3 matrix
*   the algorithm is reference the method of Matlab (inv)
*/

template<class T>
double** inv3x3(T** inx){

    double** x = new2d<double>(3,3);
    for(int i = 0, j = 0 ; i < 3 ; ++i)
        for(j=0 ; j < 3 ; ++j)
            x[i][j] = inx[i][j];

//    printf("calculate inv3x3\n");
//    print2DMatrix<T>(inx,0,3,0,3);
    int p1 = 0, p2 = 3, p3 = 6,tmp;
    double absx11 = abs(x[0][0]), absx21 = abs(x[1][0]), absx31 = abs(x[2][0]), tmp_T, tmp_T1;
    double* tmp_ptr_T;

    if(absx21 > absx11 && absx21 > absx31){
        // swap row 1 and row 2
        tmp = p1;
        p1 = p2;
        p2 = tmp;
//        p1^=p2^=p1^=p2;
        tmp_ptr_T = x[0];
        x[0] = x[1];
        x[1] = tmp_ptr_T;

    }else if(absx31 > absx11){
        // swap row 1 and row 3
        tmp = p1;
        p1 = p3;
        p3 = tmp;
//        p1^=p3^=p1^=p3;
        tmp_ptr_T = x[0];
        x[0] = x[2];
        x[2] = tmp_ptr_T;
    }

    // First opportunity to compute save 1 / x(1,1).
    x[1][0] /= x[0][0];
    x[2][0] /= x[0][0];
    x[1][1] -= x[1][0]*x[0][1];
    x[2][1] -= x[2][0]*x[0][1];
    x[1][2] -= x[1][0]*x[0][2];
    x[2][2] -= x[2][0]*x[0][2];

    if(abs(x[2][1]) > abs(x[1][1])){
        //swap row 2 and row 3
        tmp = p2;
        p2 = p3;
        p3 = tmp;
//        p2^=p3^=p2^=p3;
        tmp_ptr_T = x[1];
        x[1] = x[2];
        x[2] = tmp_ptr_T;
    }

    // First opportunity to compute and save 1 / x(2,2).
    x[2][1] /= x[1][1];
    x[2][2] -= x[2][1]*x[1][2];
    // Several opportunities here to replace divisions with mults
    // by saved reciprocal values of x(1,1), x(2,2), and x(3,3).
    double** y = new2d<double>(3,3);
    tmp_T = (x[2][1]*x[1][0]-x[2][0])/x[2][2];
    tmp_T1 = (-x[1][0]-x[1][2]*tmp_T)/x[1][1];
    y[p1%3][p1/3] = (1 -x[0][1]*tmp_T1 - x[0][2]*tmp_T)/x[0][0];
    y[(p1+1)%3][(p1+1)/3] = tmp_T1;
    y[(p1+2)%3][(p1+2)/3] = tmp_T;

    tmp_T = -x[2][1]/x[2][2];
    tmp_T1 = (1-x[1][2]*tmp_T)/x[1][1];
    y[p2%3][p2/3] = (-x[0][1]*tmp_T1 - x[0][2]*tmp_T)/x[0][0];
    y[(p2+1)%3][(p2+1)/3] = tmp_T1;
    y[(p2+2)%3][(p2+2)/3] = tmp_T;

    tmp_T = 1/x[2][2];
    tmp_T1 = (-x[1][2]*tmp_T)/x[1][1];
    y[p3%3][p3/3] = (-x[0][1]*tmp_T1 - x[0][2]*tmp_T)/x[0][0];
    y[(p3+1)%3][(p3+1)/3] = tmp_T1;
    y[(p3+2)%3][(p3+2)/3] = tmp_T;

    delete2d<double>(x);
//    printf("done\n");
//    print2DMatrix<double>(y,0,3,0,3);
    return y;
}


double exp4(double x) {
    return (720+x*(720+x*(360+x*(120+x*(30+x*(6+x))))))*0.0013888888f;
}

#endif // __UTILS_H
