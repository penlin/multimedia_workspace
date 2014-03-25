#ifndef __DATA_ALLOC_H
#define __DATA_ALLOC_H

template<class T>
T** new2d(const int &r, const int &c){
    T **ret = new T*[r];
    ret[0] = new T[r*c];
    for(int i = 1 ; i < r ; ++i)
        ret[i] = ret[i-1] + c;
    return ret;
}

template<class T>
T** new2d(const int &r, const int &c , const T &value){
    int i = 0;
    T **ret = new T*[r];
    ret[0] = new T[r*c];
    for(i = 1 ; i < r ; ++i)
        ret[i] = ret[i-1] + c;
    for(i = 0 ; i < r*c ; ++i)
        ret[0][i] = value;
    return ret;
}

template<class T>
void delete2d(T** a) {
    delete a[0];
    delete a;
}


template<class T>
T*** new3d(const int &r, const int &c, const int &d) {
    printf("1");
//    T ***ret = new T**[r];
    T*** ret = (T***)malloc(sizeof(T**)*r);
    printf("2");
//    ret[0] = new T*[r*c];
    ret[0] = (T**)malloc(sizeof(T*)*r*c);
    printf("3");
//    ret[0][0] = new T[r*c*d];
    ret[0][0] = (T*)malloc(sizeof(T)*r*c*d);
    printf("4");

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
    T ***ret = new T**[r];
    ret[0] = new T*[r*c];
    ret[0][0] = new T[r*c*d];

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
    for(i=0; i < r*c*d ; ++i)
        ret[0][0][i] = value;

    return ret;
}

template<class T>
void delete3d(T*** a) {
//    delete a[0][0];
//    delete a[0];
//    delete a;
    free(a[0][0]);
    free(a[0]);
    free(a);
}

#endif
