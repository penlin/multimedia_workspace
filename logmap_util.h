#ifndef __LOG_MAP_H
#define __LOG_MAP_H

#include<stdio.h>
#include<math.h>
#include "image_process_utils.h"
#include "utils.h"

#define EPS 1e-50
#define INFTY -1e2;

/*
* @Penlin:
* Le1=-log(1+exp(Lu)); Le2=Lu+Le1; % ln(exp((u+1)/2*Lu)/(1+exp(Lu)))
*/
/** depreciated @Penlin: 2013/07/30 **/
void computeLe(const double* Lu, double* Le1, double* Le2, const int &lu){

    for(int i = 0 ; i < lu ; ++i){
        Le1[i] = -log(1+exp(Lu[i]));
        Le2[i] = (Le1[i] + Lu[i]);
    }
}


/**
*   @Penlin: BCJR log-MAP algorithm
*/

void ComputeGamma(double*** gamma,const int &lu ,const int &Ns , int **ps, int **pout, double *Ly, double *Le1, double * Le2)
{
    int i,j;
    double Lyk[2];

    for(i=1;i<(lu+1);++i){
        Lyk[0] = Ly[2*i-2];
        Lyk[1] = Ly[2*i-1];
        for(j=0;j<Ns;++j){
             gamma[ps[j][0]][j][i] = -Lyk[0] + Lyk[1]*pout[j][1] + Le1[i-1] ;
             gamma[ps[j][1]][j][i] =  Lyk[0] + Lyk[1]*pout[j][3] + Le2[i-1] ;
        }
    }
}

void ComputeAlpha(double *tmpMax, double ***gamma, double **alpha,const int &lu,const int &Ns)
{
    int i, j, k;
    double tmp_alpha;

    for(i=1;i<lu;++i){
        tmpMax[i] = 1e-10;
        for(j=0;j<Ns;++j){
            tmp_alpha = 0;
            for(k=0;k<Ns;++k){
                tmp_alpha += exp(gamma[k][j][i] + alpha[i-1][k]);
            }
            if(tmp_alpha <= EPS){
                alpha[i][j] = INFTY;
            }else
                alpha[i][j] = log(tmp_alpha);
            if(alpha[i][j]>tmpMax[i])
                tmpMax[i] = alpha[i][j];
        }

        for(j=0;j<Ns;++j)
            alpha[i][j] -= tmpMax[i];
    }
}

void ComputeBeta(double *tmpMax, double ***gamma, double **beta,const int &lu, const int &Ns)
{
    int i, j, k;
    double tmp_beta;

    for(i=lu-1;i>=0;--i){
        for(j=0;j<Ns;++j){
            tmp_beta = 0;
            for(k=0;k<Ns;++k){
                tmp_beta += exp(gamma[j][k][i+1] + beta[i+1][k]);
            }
            if(tmp_beta <= EPS){
                beta[i][j] = INFTY;
            }else{
                beta[i][j] = log(tmp_beta);
            }
        }
        for(j=0;j<Ns;++j)
            beta[i][j] -= tmpMax[i];
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
                tmp1 += exp(gamma[ps[j][0]][j][i+1] + alpha[i][ps[j][0]] + beta[i+1][j]);
                tmp2 += exp(gamma[ps[j][1]][j][i+1] + alpha[i][ps[j][1]] + beta[i+1][j]);
        }
        La[i] = log((tmp2+EPS)/(tmp1+EPS));// - log(tmp1+EPS);
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

void BCJR_logmap(const int &Ns, const int &lu, const int &ind_dec, double* Ly, double* Le1, double* Le2, int** ps, int** pout, double* LA){

//        printf("allocate buffer ... %lf \n",getCurrentTime());

        int i, j;
        double** Alpha = new2d<double>(lu,Ns,0);
        double** Beta = new2d<double>(lu+1,Ns,0);
        double*** gam = new3d<double>(Ns,Ns,lu+1,-1e2);
        double* tmp;

        const double _logNs = -log(Ns);

        // initialize the alpha, beta, gamma
//        printf("initialize ... %lf \n",getCurrentTime());

        for(j=1 ; j<Ns ; ++j)
            Alpha[0][j]= INFTY;

        for(j=1 ; j<Ns ; ++j){
            if(ind_dec){
                Beta[lu][j] = INFTY;
            }else{
                Beta[lu][j] = _logNs;
            }
        }

        for(i=0 ; i<Ns ; ++i)
            for(j=0 ; j<Ns ; ++j)
                gam[i][j][0]= 0;

//        printf("computeGamma ... %lf \n",getCurrentTime());
        // compute gamma at every depth level (stage)
        ComputeGamma(gam,lu,Ns,ps,pout,Ly,Le1,Le2);

        // compute alpha in forward recursive
//        printf("computeAlpha ... %lf \n",getCurrentTime());
        tmp = (double*)malloc(sizeof(double)*lu);
        ComputeAlpha(tmp,gam,Alpha,lu,Ns);

        // compute beta in backward recursive
//        printf("computeBeta ... %lf \n",getCurrentTime());
        ComputeBeta(tmp,gam,Beta,lu,Ns);

        // compute the soft output LLR
//        printf("computeSoftOutput ... %lf \n",getCurrentTime());
        ComputeSoftOutput(LA,Alpha,Beta,gam,ps,Ns,lu);

        // free memory
//        printf("free buffer ... %lf \n",getCurrentTime());
        delete2d<double>(Alpha);
        delete2d<double>(Beta);
        delete3d<double>(gam);
        free(tmp);
//        printf("done!%lf \n",getCurrentTime());
}



#endif // __LOG_MAP_H
