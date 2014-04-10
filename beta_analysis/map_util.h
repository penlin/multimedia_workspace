#ifndef __BCJR_MAP_H
#define __BCJR_MAP_H

#include<stdio.h>
#include<math.h>
#include "image_process_utils.h"
#include "utils.h"

/*
* @Penlin:
* Le1= P0; Le2=P1;
*/

void computeLe(const double* Lu, double* Le1, double* Le2, const int &lu){

    for(int i = 0 ; i < lu ; ++i){
        Le1[i] = 1/(1+exp(Lu[i]));      // P0
        Le2[i] = 1-Le1[i];              // P1
    }
}


/**
*   @Penlin: BCJR MAP algorithm
*/

void ComputeGamma(double*** gamma,const int &lu ,const int &Ns , int **ps, int **pout, double *Ly, double *Le1, double * Le2)
{
    int i,j;
    double Lyk[2],tmp=0;

    for(i=1;i<(lu+1);++i){
//        Lyk[0] = 1/(1+exp(-Ly[2*i-2])); // P=1
//        Lyk[1] = 1/(1+exp(-Ly[2*i-1])); // P=1
        Lyk[0] = Ly[2*i-2];
        Lyk[1] = Ly[2*i-1];
        tmp = 0;
        for(j=0;j<Ns;++j){
//             gamma[ps[j][0]][j][i] = (1-Lyk[0])*((1-pout[j][1])/2+Lyk[1]*pout[j][1])*Le1[i-1] ;
//             gamma[ps[j][1]][j][i] =  Lyk[0]*((1-pout[j][3])/2+Lyk[1]*pout[j][3])*Le2[i-1];
            gamma[ps[j][0]][j][i] = exp(-Lyk[0]+pout[j][1]*Lyk[1])*Le1[i-1];
            gamma[ps[j][1]][j][i] = exp(Lyk[0]+pout[j][3]*Lyk[1])*Le2[i-1];
//            tmp = gamma[ps[j][0]][j][i] + gamma[ps[j][1]][j][i];
//            gamma[ps[j][0]][j][i]/=tmp;
//            gamma[ps[j][1]][j][i]/=tmp;
        }
    }
}

void ComputeAlpha(double ***gamma, double **alpha,const int &lu,const int &Ns)
{
    int i, j, k;
    double tmp_alpha = 0;

    for(i=1;i<lu;++i){
        tmp_alpha = 0;
        for(j=0;j<Ns;++j){
            alpha[i][j] = 0;
            for(k=0;k<Ns;++k){
                alpha[i][j]+=gamma[k][j][i]*alpha[i-1][k];
            }
            tmp_alpha+=alpha[i][j];
        }

        for(j=0;j<Ns;++j)
            alpha[i][j]/= tmp_alpha;
    }
}

void ComputeBeta(double ***gamma, double **beta,const int &lu, const int &Ns)
{
    int i, j, k;
    double tmp_beta;

    for(i=lu-1;i>=0;--i){
        tmp_beta = 0;
        for(j=0;j<Ns;++j){
            beta[i][j] = 0;
            for(k=0;k<Ns;++k){
                beta[i][j] += gamma[j][k][i+1]*beta[i+1][k];
            }

            tmp_beta+=beta[i][j];
        }
        for(j=0;j<Ns;++j)
            beta[i][j]/= tmp_beta;
    }
}

void ComputeSoftOutput(double *La, double **alpha, double **beta, double ***gamma, int **ps,const int &Ns,const int &lu)
{
    double tmp1, tmp2;
    int i,j;

    for(i=0;i<lu;++i){
        tmp1 = 0;
        tmp2 = 0;
        for(j=0;j<Ns;++j){
//            for(k=0;k<Ns;++k){
                tmp1 += (gamma[ps[j][0]][j][i+1]*alpha[i][ps[j][0]]*beta[i+1][j]);
                tmp2 += (gamma[ps[j][1]][j][i+1]*alpha[i][ps[j][1]]*beta[i+1][j]);
//            }
        }
        La[i] = log(tmp2/tmp1);// - log(tmp1+EPS);
    }
}

/**
*   @Penlin: FBA (forward backward algorithm) to solve BCJR
*
*   @param Ns           [constant]
*   @param lu           [constant]
*   @param ind_dec      [boolean constant]
*   @param Ly           [2*lu]
*   @param Le1          [lu]
*   @param Le2          [lu]
*   @param ps           [Ns*2]
*   @param pout         [Ns*4]
*
*   @output LA
*/

void BCJR_map(const int &Ns, const int &lu, const int &ind_dec, double* Ly, double* Le1, double* Le2, int** ps, int** pout, double* LA){

        int i, j;
        double** Alpha = new2d<double>(lu,Ns,0);
        double** Beta = new2d<double>(lu+1,Ns,0);
        double*** gam = new3d<double>(Ns,Ns,lu+1,0);

        const double _logNs = 1/Ns;

        // initialize the alpha, beta, gamma
        Alpha[0][0] = 0.5;
        for(j=1 ; j<Ns ; ++j)
            Alpha[0][j]= 0;

        Beta[lu][0] = 0.5;
        for(j=1 ; j<Ns ; ++j){
            if(ind_dec){
                Beta[lu][j] = 0;
            }else{
                Beta[lu][j] = _logNs;
            }
        }

        for(i=0 ; i<Ns ; ++i)
            for(j=0 ; j<Ns ; ++j)
                gam[i][j][0]= 0.5;

        // compute gamma at every depth level (stage)
        ComputeGamma(gam,lu,Ns,ps,pout,Ly,Le1,Le2);

        // compute alpha in forward recursive
        ComputeAlpha(gam,Alpha,lu,Ns);

        // compute beta in backward recursive
        ComputeBeta(gam,Beta,lu,Ns);

        // compute the soft output LLR
        ComputeSoftOutput(LA,Alpha,Beta,gam,ps,Ns,lu);

        // free memory
        delete2d<double>(Alpha);
        delete2d<double>(Beta);
        delete3d<double>(gam);
}



#endif // __BCJR_MAP_H
