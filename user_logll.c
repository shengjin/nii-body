////////////
#include <stdio.h>
#include <math.h>
#include <stdlib.h>
#include <string.h>
#include "nbody.h"
#include "alloc.h"


//int nbody1(double para0, double para1, double para2, double para3, double para4, double para5, double para6, int nline_data, double *tgrid_days, double **simul_ra_dec);
int nbody2(double para0, double para1, double para2, double para3, double para4, double para5, double para6, double para7, double para8, double para9, double para10, double para11, double para12, double para13, int nline_data, double* tgrid_days, double** simul_ra_dec);

extern int ndim_data; 

//////////////////////////////
// NOTE: change logll_beta in mpi_batch/init if func_prototype is changed.
double logll_beta(double *ptr_one_chain, int nline_data, double *data_NlineNdim, double beta_one);

// N_parm not explicitly declare here because we will use each of parms invery detail, thus for sure we know the number.
double logll_beta(double *ptr_one_chain, int nline_data, double *data_NlineNdim, double beta_one)
{
   /* model description:
    * 
    * single planet , N_parm = 14
    * 
    */
    //
    const double PI = 3.141592653589793;

    double para0;
    double para1;
    double para2;
    double para3;
    double para4;
    double para5;
    double para6;
    double para7;
    double para8;
    double para9;
    double para10;
    double para11;
    double para12;
    double para13;
    double para14;
    
    para0 = *ptr_one_chain;
    para1 = *(ptr_one_chain+1);
    para2 = *(ptr_one_chain+2);
    para3 = *(ptr_one_chain+3);
    para4 = *(ptr_one_chain+4);
    para5 = *(ptr_one_chain+5);
    para6 = *(ptr_one_chain+6);
    para7 = *(ptr_one_chain+7);
    para8 = *(ptr_one_chain+8);
    para9 = *(ptr_one_chain+9);
    para10 = *(ptr_one_chain+10);
    para11 = *(ptr_one_chain+11);
    para12 = *(ptr_one_chain+12);
    para13 = *(ptr_one_chain+13);
    para14 = *(ptr_one_chain+14);

    double *tgrid_days;
    tgrid_days = alloc_1d_double(nline_data);
    //
    double **simul_ra_dec;
    simul_ra_dec = alloc_2d_double(nline_data, 2);
    //
    double **obs_ra_dec; // add noise term latter
    obs_ra_dec = alloc_2d_double(nline_data, 2);
    //
    for (int i=0; i<nline_data; i++)
    {
        tgrid_days[i] = data_NlineNdim[i*ndim_data];
        obs_ra_dec[i][0] = data_NlineNdim[i*ndim_data+1];
        obs_ra_dec[i][1] = data_NlineNdim[i*ndim_data+2];
    }

  
   // nbody1(para0, para1, para2, para3, para4, para5, para6, nline_data, tgrid_days, simul_ra_dec);
    nbody2(para0, para1, para2, para3, para4, para5, para6, para7, para8, para9, para10, para11, para12, para13, nline_data, tgrid_days, simul_ra_dec);
    /*
    if (beta_one == 1.0)
    {
        for (int i=0; i<nline_data; i++)
        {
            printf("%lf %lf %lf %lf %lf\n", tgrid_days[i], obs_ra_dec[i][0], obs_ra_dec[i][1], simul_ra_dec[i][0], simul_ra_dec[i][1]);
        }
    }
    */

    //2, compare caculated curve with i30.dat, calcualte the difference, which is logll 

    double logll;

    double logll_ra, logll_dec;

    double sig_power;
    sig_power = para14*para14;

    double PI2 = 6.28318530716;
    double AC_twice_all;
    AC_twice_all = -0.5*log(PI2)*(double)nline_data*2.0 -0.5*log(sig_power)*(double)nline_data*2.0;

    logll = AC_twice_all;

    for (int i=0; i<nline_data; i++)
    {
        // printf("idx %d ok5 \n", i);
        logll_ra =  - pow((obs_ra_dec[i][0]-simul_ra_dec[i][0]), 2.0)/2.0/sig_power;
        logll_dec =  - pow((obs_ra_dec[i][1]-simul_ra_dec[i][1]), 2.0)/2.0/sig_power;
        // printf("err2 %.3f \n", (pow( data_NlineNdim[i*ndim_data+2], 2.0)));
        logll = logll + logll_ra + logll_dec;
    }

    // final tempering
    logll = logll*beta_one;
    //printf("%lf %lf \n", logll, beta_one);

    free_1d_double(tgrid_days);
    tgrid_days = NULL;
    free_2d_double(simul_ra_dec);
    simul_ra_dec = NULL;
    free_2d_double(obs_ra_dec);
    obs_ra_dec = NULL;
    //    
    return logll;
}















