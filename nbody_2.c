#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>
#include "nbody.h"
#include "alloc.h"

#include <time.h>

double xmass[MAX_BODIES]; // Mass, in unit of Solar mass

int nbody2(double para0, double para1, double para2, double para3, double para4, double para5, double para6, double para7, double para8, double para9, double para10, double para11, double para12, double para13, int nline_data, double* tgrid_days, double** simul_ra_dec);

const double stellar_mass = 1.022; //the stellar mass in unit of Solar mass
const double Mearth_to_Msun = 3e-6;

const double yr_days = 365.242199074;

const double degrees_2pi = 360.0;

const double pc_to_au = 206265.0; /* 1 pc = 206265 AU */

const double tpi = 2.0 * 3.1415926535897932;
const double pi = 3.1415926535897932;
const double rad2uas = (180.0 * 3600.0 * 1e6) / pi; /* 1 rad -> micro-arcsec */

const int nbod = 3; // number of bodies
const int ndim = 6; // dimension of integration

const double h_init_scale = 1.0;  // init step size in days
const double err_tolerance = 1e-6; 
const double t_start = 0.0;


// the distance of the target system from observers on Earth 
const double observer_distance_pc = 5.0; //unit of parsec


// RKF78 in case derivative function does not explicitly depend on time t.
int rkf78_no_t(double* ptr_t, double* ptr_dt, double err_tol, double* x, int n);
// RKF78 in case derivative function depend on time t, e.g., Earth's satelites.
//int rkf78(double* ptr_t, double* ptr_dt, double err_tol, double* x, int n);

// angular conversion 
double zatan(double sina, double cosa);

// orbit elements to Cartesian coordinates
int ele2rat(double a, double e, double ci, double w, double omega, double cm0, int ii, double xout[6]);

// Cartesian to elements
int rat2ele(double xin[6], double* a, double* e, double* ci, double* w, double* omega, double* cm, int ii, int ito);



//int nbody(double para0, double para1, double para2, double para3, double para4, double para5, double para6, int nline_data, double *tgrid_days, double **simul_ra_dec)
int nbody2(double para0, double para1, double para2, double para3, double para4, double para5, double para6, double para7, double para8, double para9, double para10, double para11, double para12, double para13, int nline_data, double* tgrid_days, double** simul_ra_dec)

{
    int i, j;


    /* arrays for elements and masses */

    double p0[MAX_BODIES], e0[MAX_BODIES], ci0[MAX_BODIES], w0[MAX_BODIES];
    double omg0[MAX_BODIES], cm00[MAX_BODIES];
    double raw_mass[MAX_BODIES];     // temperory copy of mass 

    // stellar mass, unit of Solor mass
    xmass[0] = stellar_mass;

    i = 1;
    p0[i] = para0;
    e0[i] = para1;
    ci0[i] = para2;
    w0[i] = para3;
    omg0[i] = para4;
    cm00[i] = para5;
    raw_mass[i] = para6;
    xmass[i] = raw_mass[i] * Mearth_to_Msun;

    i = 2;
    p0[i] = para7;
    e0[i] = para8;
    ci0[i] = para9;
    w0[i] = para10;
    omg0[i] = para11;
    cm00[i] = para12;
    raw_mass[i] = para13;

    xmass[i] = raw_mass[i] * Mearth_to_Msun;
    
    //FILE* faei = fopen("c_aei.txt", "w");
    //FILE* fd = fopen("c_d.txt", "w");
    //FILE* fra_dec = fopen("ra_dec.txt", "w");
    //if (!fra_dec) { perror("ra_dec.txt"); return 1; }
    //if (!fd) { perror("c_d.txt"); return 1; }

    // calc semi_major axes, convert degree to radian 
    double a0[MAX_BODIES], a[MAX_BODIES];
    double e[MAX_BODIES], ci[MAX_BODIES], w[MAX_BODIES], omg[MAX_BODIES], cm0[MAX_BODIES], cm[MAX_BODIES];
    for (i = 1; i < nbod; i++) {
        a0[i] = pow(pow(p0[i] / yr_days, 2) * xmass[0], 1.0 / 3.0);
        a[i] = a0[i];
        e[i] = e0[i];
        ci[i] = ci0[i] * tpi / degrees_2pi; 
        w[i] = w0[i] * tpi / degrees_2pi;
        omg[i] = omg0[i] * tpi / degrees_2pi;
        cm0[i] = cm00[i] * tpi / degrees_2pi;
    }

    /* print check
    printf("Initial a, e, mass (bodies 2..%d):\n", nbod);
    for (i = 1; i < nbod; i++) {
        printf("% .6f  % .6f  % .6f\n", a[i], e[i], xmass[i]);
    }
     */

    //initialzie the state vector
    double x[MAX_BODIES][6] = { {0} };  // x,y,z and velocity x, v_y, v_z
    double xp[6];                       // temporary store the positions of one body

    double xrk[MAX_BODIES * 6];         // 1-dimension array used for rkf78 integrator



    // initialize the cartesian grid
    for (i = 1; i < nbod; i++) {
        ele2rat(a[i], e[i], ci[i], w[i], omg[i], cm0[i], i, xp);
        for (j = 0; j < ndim; j++) x[i][j] = xp[j];
    }

    // calibration of the center of mass
    double dx[6] = { 0 };
    for (j = 0; j < ndim; j++)
    {
        dx[j] = 0.0;
        for (i = 1; i < nbod; i++)
        {
            dx[j] += x[i][j] * xmass[i];
        }
    }
    double tmass = 0.0;
    for (i = 0; i < nbod; i++) tmass += xmass[i];
    for (j = 0; j < ndim; j++) x[0][j] = -dx[j] / tmass;

    // re-calculate the positions relative to the center of mass
    for (i = 1; i < nbod; i++)
    {
        for (j = 0; j < ndim; j++)
        {
            x[i][j] += x[0][j];
        }
    }

    // convert to 1-dimension array for the rkf78 integrator
    for (i = 0; i < nbod; i++)
    {
        for (j = 0; j < ndim; j++)
        {
            xrk[i * ndim + j] = x[i][j];
        }
    }

    // control parameters for the rkf78 integrator
    
    double h = tpi / yr_days * h_init_scale;    //  initial step size
    double err = err_tolerance;                 // err tolerance
    double t = t_start;                     //  initial time
    double tstop;
    tstop = (tgrid_days[nline_data-1])/ yr_days * tpi;    // total integration time

    // control para (1 means yes; or 0 means no): should we generate output in the next timestep  
    int output_next_step = 1;

    int i_out_ingrid = 0;
    double tout = tgrid_days[i_out_ingrid]/ yr_days * tpi;    // first output time

    double dist_au = observer_distance_pc * pc_to_au;

    /* main loop of the integration */
    while (t < tstop) {
        // update the coordinate
        for (i = 0; i < nbod; i++)
        {
            for (j = 0; j < ndim; j++)
            {
                x[i][j] = xrk[i * ndim + j];
            }
        }

        /* 
        for (i = 1; i < nbod; i++)
        {
            for (j = 0; j < ndim; j++)
            {
                xp[j] = x[i][j] - x[0][j];
            }
            rat2ele(xp, &a[i], &e[i], &ci[i], &w[i], &omg[i], &cm[i], i, 0);
        }
	*/

        // 
        if ( (output_next_step == 1) || (i_out_ingrid == nline_data-1) )
	{

            double t_days = (t / tpi) * yr_days; /* t -> days */
            //printf("%.6f %lf", t, tout);

            // compute the RA and DEC in astrometry, unit of micro arcsecond
            double sx = x[0][0];   // stellar X, unit of AU
            double sy = x[0][1];   // stellar Y
            double sz = x[0][2];   // stellar Z

            double star_ra_rad = sx / dist_au;
            double star_dec_rad = sy / dist_au;

            double star_ra_uas = star_ra_rad * rad2uas;
            double star_dec_uas = star_dec_rad * rad2uas;


            simul_ra_dec[i_out_ingrid][0] = star_ra_uas;
            simul_ra_dec[i_out_ingrid][1] = star_dec_uas;
            /* star_ra_uas, star_dec_uas */
            //printf("%20.6f  %20.6f\n", simul_ra_dec[i_out_ingrid][0], simul_ra_dec[i_out_ingrid][1]);
            //fflush(fra_dec);

            if (i_out_ingrid < nline_data -1)
	    {
                i_out_ingrid ++;
                tout = tgrid_days[i_out_ingrid]/ yr_days * tpi;   // next output time
                //printf("%d %lf\n", i_out_ingrid, tgrid_days[i_out_ingrid]);
            }
        }



        if (t+h >= tout)  // the integrator will past the next time grid in the next step
	{
	    h = tout - t;
	    //
            // 
            //rkf78(&t, &h, err, xrk, nbod * ndim);
            rkf78_no_t(&t, &h, err, xrk, nbod * ndim);
	    //
	    output_next_step = 1;
	}
	else
	{
            //rkf78(&t, &h, err, xrk, nbod * ndim);
            rkf78_no_t(&t, &h, err, xrk, nbod * ndim);
	    //
	    output_next_step = 0;
	}

/*
*/

        // stability check
        for (i = 1; i < nbod; i++) {
            if (fabs(a[i] - a0[i]) > 0.1 * a0[i]) { break; }
        }
    }

    /* to compute the running time
    end = clock();
    cpu_time_used = ((double)(end - start)) / CLOCKS_PER_SEC;
    printf("Program run time: %f s\n", cpu_time_used);
    */

    //fclose(faei);
    //fclose(fd);
    //fclose(fra_dec);

    return 0;
}

