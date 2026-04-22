#include <stdio.h>
#include <stdlib.h>
#include <math.h>


// derivative function used in integration
int f_deriv(double *x, double *y);


/////////////////////////////////////////////
// example derivative function used in integration
int yhc(double t, double *x, double *y);


 
// --- scratch buffer that stays alive between calls ---
// 16 slots of size n: x0, x1, y1..y13, one spare
// only reallocated if n grows, so we can reuse it for multiple calls with the same n
 
#define RKF78_NBUF 16
 
static double *rkf78_scratch = NULL;
static int     rkf78_scratch_n = 0;
 
static inline double *rkf78_get_scratch(int n)
{
    if (n <= rkf78_scratch_n) return rkf78_scratch;
 
    free(rkf78_scratch);
    size_t bytes = (size_t)RKF78_NBUF * (size_t)n * sizeof(double);
    size_t aligned_bytes = (bytes + 63u) & ~(size_t)63u;
    rkf78_scratch = (double *)aligned_alloc(64, aligned_bytes);
    if (!rkf78_scratch) {
        fprintf(stderr, "rkf78: aligned_alloc failed for %zu bytes\n", aligned_bytes);
        exit(EXIT_FAILURE);
    }
    rkf78_scratch_n = n;
    return rkf78_scratch;
}
// step controller knobs
// 1/8 exponent because the embedded error is O(h^8)
// safety < 1 so we don't ride right on the edge and get rejected next step
// fac_min/max keep dt from doing anything crazy in one go
 
#define RKF78_SAFETY     0.9
#define RKF78_FAC_MIN    0.2
#define RKF78_FAC_MAX    5.0
#define RKF78_ORDER_INV  0.125
 
//  autonomous version (no explicit time in RHS) — this is what nbody calls not the time dependent one

int rkf78_no_t(double *ptr_t, double *ptr_dt, double err_tol, double *x, int n)
{
    // grab the persistent buffer, carve it into named slices
    // restrict tells the compiler these don't alias so it can vectorize
    double *restrict buf = rkf78_get_scratch(n);
 
    double *restrict x0  = buf + 0  * n;
    double *restrict x1  = buf + 1  * n;
    double *restrict y1  = buf + 2  * n;
    double *restrict y2  = buf + 3  * n;
    double *restrict y3  = buf + 4  * n;
    double *restrict y4  = buf + 5  * n;
    double *restrict y5  = buf + 6  * n;
    double *restrict y6  = buf + 7  * n;
    double *restrict y7  = buf + 8  * n;
    double *restrict y8  = buf + 9  * n;
    double *restrict y9  = buf + 10 * n;
    double *restrict y10 = buf + 11 * n;
    double *restrict y11 = buf + 12 * n;
    double *restrict y12 = buf + 13 * n;
    double *restrict y13 = buf + 14 * n;
 
    // save x so we can rewind if the step gets rejected and not start over
    for (int i = 0; i < n; i++) x0[i] = x[i];
    const double t0 = *ptr_t;
 
    // work with local dt, no point dereferencing a pointer every loop
    double dt = *ptr_dt;
    double err_now = HUGE_VAL; // just needs to be > err_tol to start but i chose big value to avoid any complicatons
    int first_try = 1;
 
    while (err_now > err_tol) {
 
        // rewind x on retry — skip first time, x is already correct
        if (!first_try) {
            for (int i = 0; i < n; i++) x[i] = x0[i];
        }
        first_try = 0;
 
        // stage 1 — no need to copy x into x1, just eval f at x directly
        f_deriv(x, y1);
 
        // stages 2–13: Fehlberg tableau
        // pulling dt/coeff out as a scalar so the inner loop is pure FMA
 
        const double c2 = dt * (2.0 / 27.0);
        for (int i = 0; i < n; i++)
            x1[i] = x[i] + c2 * y1[i];
        f_deriv(x1, y2);
 
        const double c3 = dt / 36.0;
        for (int i = 0; i < n; i++)
            x1[i] = x[i] + c3 * (y1[i] + 3.0 * y2[i]);
        f_deriv(x1, y3);
 
        const double c4 = dt / 24.0;
        for (int i = 0; i < n; i++)
            x1[i] = x[i] + c4 * (y1[i] + 3.0 * y3[i]);
        f_deriv(x1, y4);
 
        const double c5 = dt / 48.0;
        for (int i = 0; i < n; i++)
            x1[i] = x[i] + c5 * (20.0 * y1[i] + 75.0 * (y4[i] - y3[i]));
        f_deriv(x1, y5);
 
        const double c6 = dt / 20.0;
        for (int i = 0; i < n; i++)
            x1[i] = x[i] + c6 * (y1[i] + 5.0 * y4[i] + 4.0 * y5[i]);
        f_deriv(x1, y6);
 
        const double c7 = dt / 108.0;
        for (int i = 0; i < n; i++)
            x1[i] = x[i] + c7 * (-25.0 * y1[i] + 125.0 * y4[i]
                                 - 260.0 * y5[i] + 250.0 * y6[i]);
        f_deriv(x1, y7);
 
        const double c8 = dt / 900.0;
        for (int i = 0; i < n; i++)
            x1[i] = x[i] + c8 * (93.0 * y1[i] + 244.0 * y5[i]
                                 - 200.0 * y6[i] + 13.0 * y7[i]);
        f_deriv(x1, y8);
 
        const double c9 = dt / 90.0;
        for (int i = 0; i < n; i++)
            x1[i] = x[i] + c9 * (180.0 * y1[i] - 795.0 * y4[i]
                                 + 1408.0 * y5[i] - 1070.0 * y6[i]
                                 + 67.0   * y7[i] + 270.0  * y8[i]);
        f_deriv(x1, y9);
 
        const double c10 = dt / 540.0;
        for (int i = 0; i < n; i++)
            x1[i] = x[i] + c10 * (-455.0 * y1[i] + 115.0  * y4[i]
                                  - 3904.0 * y5[i] + 3110.0 * y6[i]
                                  - 171.0  * y7[i] + 1530.0 * y8[i]
                                  - 45.0   * y9[i]);
        f_deriv(x1, y10);
 
        const double c11 = dt / 4100.0;
        for (int i = 0; i < n; i++)
            x1[i] = x[i] + c11 * (2383.0 * y1[i] - 8525.0 * y4[i]
                                  + 17984.0 * y5[i] - 15050.0 * y6[i]
                                  + 2133.0  * y7[i] + 2250.0  * y8[i]
                                  + 1125.0  * y9[i] + 1800.0  * y10[i]);
        f_deriv(x1, y11);
 
        const double c12 = dt / 4100.0;
        for (int i = 0; i < n; i++)
            x1[i] = x[i] + c12 * (60.0 * y1[i] - 600.0 * y6[i]
                                  - 60.0 * y7[i]
                                  + 300.0 * (y9[i] - y8[i] + 2.0 * y10[i]));
        f_deriv(x1, y12);
 
        const double c13 = dt / 4100.0;
        for (int i = 0; i < n; i++)
            x1[i] = x[i] + c13 * (-1777.0 * y1[i] - 8525.0 * y4[i]
                                  + 17984.0 * y5[i] - 14450.0 * y6[i]
                                  + 2193.0  * y7[i] + 2550.0  * y8[i]
                                  + 825.0   * y9[i] + 1200.0  * y10[i]
                                  + 4100.0  * y12[i]);
        f_deriv(x1, y13);
 
        // 8th-order solution — the actual answer we keep
        const double cF = dt / 840.0;
        for (int i = 0; i < n; i++)
            x[i] += cF * (272.0 * y6[i]
                          + 216.0 * (y7[i] + y8[i])
                          +  27.0 * (y9[i] + y10[i])
                          +  41.0 * (y12[i] + y13[i]));
 
        // error: difference between 7th and 8th order
        // just need the max component, L2 norm would be overkill
        const double cE = dt * 41.0 / 810.0;
        double e = 0.0;
        for (int i = 0; i < n; i++) {
            double diff = cE * (y1[i] + y11[i] - y12[i] - y13[i]);
            double a = fabs(diff);
            if (a > e) e = a;
        }
        err_now = e;
 
        // pick optimal next dt instead of blind halving
        // factor = 0.9 * (tol/err)^(1/8), clamped to [0.2, 5.0]
        double factor;
        if (err_now == 0.0) {
            factor = RKF78_FAC_MAX;  // zero error, open the throttle all the way up
        } else {
            factor = RKF78_SAFETY * pow(err_tol / err_now, RKF78_ORDER_INV);
            if (factor < RKF78_FAC_MIN) factor = RKF78_FAC_MIN;
            if (factor > RKF78_FAC_MAX) factor = RKF78_FAC_MAX;
        }
 
        if (err_now > err_tol) {
            // too sloppy, shrink and retry
            dt *= factor;
        } else {
            // good step, advance time, grow dt for next call, done
            *ptr_t = t0 + dt;
            dt *= factor;
            break;
        }
    }
 
    *ptr_dt = dt;
    return 0;
}


/////////////////////////////////////// 
// Pruned version of the RKF7(8) algorithm for cases without time dependence.
int rkf78_no_t_old(double *ptr_t, double *ptr_dt, double err_tol, double *x, int n)
{
    double err_ratio = 1e-3;
    double err_now_init = 1.0;

    double *x0 = (double *)malloc(n * sizeof(double));

    for (int i = 0; i < n; i++) {
        x0[i] = x[i];
    }
    double t0;
    t0 = *ptr_t;

    double err_now;
    err_now = err_now_init;

    while (err_now > err_tol) {
        *ptr_t = t0;
        for (int i = 0; i < n; i++) {
            x[i] = x0[i];
        }

        double *x1 = (double *)malloc(n * sizeof(double));
        double *y1 = (double *)malloc(n * sizeof(double));
        double *y2 = (double *)malloc(n * sizeof(double));
        double *y3 = (double *)malloc(n * sizeof(double));
        double *y4 = (double *)malloc(n * sizeof(double));
        double *y5 = (double *)malloc(n * sizeof(double));
        double *y6 = (double *)malloc(n * sizeof(double));
        double *y7 = (double *)malloc(n * sizeof(double));
        double *y8 = (double *)malloc(n * sizeof(double));
        double *y9 = (double *)malloc(n * sizeof(double));
        double *y10 = (double *)malloc(n * sizeof(double));
        double *y11 = (double *)malloc(n * sizeof(double));
        double *y12 = (double *)malloc(n * sizeof(double));
        double *y13 = (double *)malloc(n * sizeof(double));

        for (int i = 0; i < n; i++) {
            x1[i] = x0[i];
        }
        f_deriv(x1, y1);

        for (int i = 0; i < n; i++) {
            x1[i] = x[i] + *ptr_dt * 2.0 / 27.0 * y1[i];
        }
        f_deriv(x1, y2);

        for (int i = 0; i < n; i++) {
             x1[i] = x[i] + *ptr_dt * (y1[i] + 3.0 * y2[i]) / 36.0;
        }
        f_deriv(x1, y3);

        for (int i = 0; i < n; i++) {
             x1[i] = x[i] + *ptr_dt * (y1[i] + 3.0 * y3[i]) / 24.0;
        }
        f_deriv(x1, y4);

        for (int i = 0; i < n; i++) {
             x1[i] = x[i] + *ptr_dt * (y1[i] * 20.0 + (-y3[i] + y4[i]) * 75.0) / 48.0;
        }
        f_deriv(x1, y5);

        for (int i = 0; i < n; i++) {
             x1[i] = x[i] + *ptr_dt * (y1[i] + y4[i] * 5.0 + y5[i] * 4.0) / 20.0;
        }
        f_deriv(x1, y6);

        for (int i = 0; i < n; i++) {
             x1[i] = x[i] + *ptr_dt * (-y1[i] * 25.0 + y4[i] * 125.0 - y5[i] * 260.0 + y6[i] * 250.0) / 108.0;
        }
        f_deriv(x1, y7);

        for (int i = 0; i < n; i++) {
             x1[i] = x[i] + *ptr_dt * (y1[i] * 93.0 + y5[i] * 244.0 - y6[i] * 200.0 + y7[i] * 13.0) / 900.0;
        }
        f_deriv(x1, y8);

        for (int i = 0; i < n; i++) {
             x1[i] = x[i] + *ptr_dt * (y1[i] * 180.0 - y4[i] * 795.0 + y5[i] * 1408.0 - y6[i] * 1070.0 + y7[i] * 67.0 + y8[i] * 270.0) / 90.0;
        }
        f_deriv(x1, y9);

        for (int i = 0; i < n; i++) {
             x1[i] = x[i] + *ptr_dt * (-y1[i] * 455.0 + y4[i] * 115.0 - y5[i] * 3904.0 + y6[i] * 3110.0 - y7[i] * 171.0 + y8[i] * 1530.0 - y9[i] * 45.0) / 540.0;
        }
        f_deriv(x1, y10);

        for (int i = 0; i < n; i++) {
             x1[i] = x[i] + *ptr_dt * (y1[i] * 2383.0 - y4[i] * 8525.0 + y5[i] * 17984.0 - y6[i] * 15050.0 + y7[i] * 2133.0 + y8[i] * 2250.0 + y9[i] * 1125.0 + y10[i] * 1800.0) / 4100.0;
        }
        f_deriv(x1, y11);

        for (int i = 0; i < n; i++) {
             x1[i] = x[i] + *ptr_dt * (y1[i] * 60.0 - y6[i] * 600.0 - y7[i] * 60.0 + (y9[i] - y8[i] + 2.0 * y10[i]) * 300.0) / 4100.0;
        }
        f_deriv(x1, y12);

        for (int i = 0; i < n; i++) {
             x1[i] = x[i] + *ptr_dt * (-y1[i] * 1777.0 - y4[i] * 8525.0 + y5[i] * 17984.0 - y6[i] * 14450.0 + y7[i] * 2193.0 + y8[i] * 2550.0 + y9[i] * 825.0 + y10[i] * 1200.0 + y12[i] * 4100) / 4100.0;
        }
        f_deriv(x1, y13);

        for (int i = 0; i < n; i++) {
             x[i] = x[i] + *ptr_dt * (y6[i] * 272.0 + (y7[i] + y8[i]) * 216.0 + (y9[i] + y10[i]) * 27.0 + (y12[i] + y13[i]) * 41.0) / 840.0;
        }

        err_now = 0;
        for (int i = 0; i < n; i++) {
             double diff;
             diff = (y1[i] + y11[i] - y12[i] - y13[i]) * (*ptr_dt) * 41.0 / 810.0;
             err_now = fmax(err_now, fabs(diff));
        }

        *ptr_t = *ptr_t + *ptr_dt;

        if (err_now > err_tol) {
             *ptr_dt = *ptr_dt / 2.0;
        }

        if (err_now < err_tol * err_ratio) {
             *ptr_dt = *ptr_dt * 2.0;
        }

        free(x1);
	x1 = NULL;
        free(y1);
	y1 = NULL;
        free(y2);
	y2 = NULL;
        free(y3);
	y3 = NULL;
        free(y4);
	y4 = NULL;
        free(y5);
	y5 = NULL;
        free(y6);
	y6 = NULL;
        free(y7);
	y7 = NULL;
        free(y8);
	y8 = NULL;
        free(y9);
	y9 = NULL;
        free(y10);
	y10 = NULL;
        free(y11);
	y11 = NULL;
        free(y12);
	y12 = NULL;
        free(y13);
	y13 = NULL;
    }

    free(x0);
    x0 = NULL;

    return 0;
}




//
//////////////////// Full RKF7(8) algorithm
int rkf78(double *ptr_t, double *ptr_dt, double err_tol, double *x, int n) {
    // Runge Kutta Fehlberg 7(8) numerical integration method
    
    /*
     t  (*ptr_t):  start time of this integration (t_init)
     dt (*ptr_dt): initial step size attempted to use (may fail due to large err)
     *
     ptr_t:        pointer to t
     ptr_dt:       pointer to dt
     err_tol:      truncation error tolerance [non-dimensional]
     *x:           integration vector at time = ti
     n:            number of differential equations (dimension of x)
     *
     RETURN:  *x with updated values at time tf.
     *
     */

    // local constants
    // 
    // criteria in the tiny err case to increase step size
    double err_ratio = 1e-4; 
    // set manually err_now to big value to start the loop
    double err_now_init = 1.0;


    //////////////////////////////
    // Create a copy of initial x and t, now x0 and t0
    //
    // allocate arrays
    double *x0 = (double *)malloc(n * sizeof(double));
    // copy initial values
    for (int i = 0; i < n; i++) 
    {
            x0[i] = x[i];
    }
    //
    // integration start time at each step
    double t0;
    t0 = *ptr_t;


    ////////////////////////////
    //  
    // err for the working step size
    double err_now;
    // set it to 1.0 at the beginning
    err_now = err_now_init;  


    /////////////////////////////////////////////////
    // Main loop
    //
    while (err_now > err_tol) 
    {

        ///////////////////
	// re-initialize the t, t1, *x
	*ptr_t = t0;
        double t1;
	t1 = *ptr_t;
	//
        for (int i = 0; i < n; i++) 
	{
            x[i] = x0[i];
        }
        
	// allocate arrays of x and integration coefficients
        double *x1 = (double *)malloc(n * sizeof(double));
        double *y1 = (double *)malloc(n * sizeof(double));
        double *y2 = (double *)malloc(n * sizeof(double));
        double *y3 = (double *)malloc(n * sizeof(double));
        double *y4 = (double *)malloc(n * sizeof(double));
        double *y5 = (double *)malloc(n * sizeof(double));
        double *y6 = (double *)malloc(n * sizeof(double));
        double *y7 = (double *)malloc(n * sizeof(double));
        double *y8 = (double *)malloc(n * sizeof(double));
        double *y9 = (double *)malloc(n * sizeof(double));
        double *y10 = (double *)malloc(n * sizeof(double));
        double *y11 = (double *)malloc(n * sizeof(double));
        double *y12 = (double *)malloc(n * sizeof(double));
        double *y13 = (double *)malloc(n * sizeof(double));


	//////////////////
	// set x1 = x0
        for (int i = 0; i < n; i++) 
	{
            x1[i] = x0[i];
        }
        // calc y1 use derivative function
        yhc(t1, x1, y1);

        //////////
        for (int i = 0; i < n; i++) 
	{
            x1[i] = x[i] + *ptr_dt * 2.0 / 27.0 * y1[i];
        }
        t1 = *ptr_t + *ptr_dt * 2.0 / 27.0;
	//
        yhc(t1, x1, y2);

        //////////
        for (int i = 0; i < n; i++) 
	{
            x1[i] = x[i] + *ptr_dt * (y1[i] + 3.0 * y2[i]) / 36.0;
        }
        t1 = *ptr_t + *ptr_dt / 9.0;
	//
        yhc(t1, x1, y3);

        //////////
        for (int i = 0; i < n; i++) 
	{
            x1[i] = x[i] + *ptr_dt * (y1[i] + 3.0 * y3[i]) / 24.0;
        }
        t1 = *ptr_t + *ptr_dt / 6.0;
	//
        yhc(t1, x1, y4);

        //////////
        for (int i = 0; i < n; i++) 
	{
            x1[i] = x[i] + *ptr_dt * (y1[i] * 20.0 + (-y3[i] + y4[i]) * 75.0) / 48.0;
        }
        t1 = *ptr_t + *ptr_dt * 5.0 / 12.0;
	//
        yhc(t1, x1, y5);

        //////////
        for (int i = 0; i < n; i++) 
	{
            x1[i] = x[i] + *ptr_dt * (y1[i] + y4[i] * 5.0 + y5[i] * 4.0) / 20.0;
        }
        t1 = *ptr_t + *ptr_dt * 1.0 / 2.0;
	//
        yhc(t1, x1, y6);

        //////////
        for (int i = 0; i < n; i++) 
	{
            x1[i] = x[i] + *ptr_dt * (-y1[i] * 25.0 + y4[i] * 125.0 - y5[i] * 260.0 + y6[i] * 250.0) / 108.0;
        }
        t1 = *ptr_t + *ptr_dt * 5.0 / 6.0;
	//
        yhc(t1, x1, y7);

        //////////
        for (int i = 0; i < n; i++) 
	{
            x1[i] = x[i] + *ptr_dt * (y1[i] * 93.0 + y5[i] * 244.0 - y6[i] * 200.0 + y7[i] * 13.0) / 900.0;
        }
        t1 = *ptr_t + *ptr_dt * 1.0 / 6.0;
	//
        yhc(t1, x1, y8);

        //////////
        for (int i = 0; i < n; i++) 
	{
            x1[i] = x[i] + *ptr_dt * (y1[i] * 180.0 - y4[i] * 795.0 + y5[i] * 1408.0 - y6[i] * 1070.0 + y7[i] * 67.0 + y8[i] * 270.0) / 90.0;
        }
        t1 = *ptr_t + *ptr_dt * 2.0 / 3.0;
	//
        yhc(t1, x1, y9);

        //////////
        for (int i = 0; i < n; i++) 
	{
            x1[i] = x[i] + *ptr_dt * (-y1[i] * 455.0 + y4[i] * 115.0 - y5[i] * 3904.0 + y6[i] * 3110.0 - y7[i] * 171.0 + y8[i] * 1530.0 - y9[i] * 45.0) / 540.0;
        }
        t1 = *ptr_t + *ptr_dt * 1.0 / 3.0;
	//
        yhc(t1, x1, y10);

        //////////
        for (int i = 0; i < n; i++) 
	{
            x1[i] = x[i] + *ptr_dt * (y1[i] * 2383.0 - y4[i] * 8525.0 + y5[i] * 17984.0 - y6[i] * 15050.0 + y7[i] * 2133.0 + y8[i] * 2250.0 + y9[i] * 1125.0 + y10[i] * 1800.0) / 4100.0;
        }
        t1 = *ptr_t + *ptr_dt;
	//
        yhc(t1, x1, y11);

        //////////
        for (int i = 0; i < n; i++) 
	{
            x1[i] = x[i] + *ptr_dt * (y1[i] * 60.0 - y6[i] * 600.0 - y7[i] * 60.0 + (y9[i] - y8[i] + 2.0 * y10[i]) * 300.0) / 4100.0;
        }
        t1 = *ptr_t;
	//
        yhc(t1, x1, y12);

        //////////
        for (int i = 0; i < n; i++) 
	{
            x1[i] = x[i] + *ptr_dt * (-y1[i] * 1777.0 - y4[i] * 8525.0 + y5[i] * 17984.0 - y6[i] * 14450.0 + y7[i] * 2193.0 + y8[i] * 2550.0 + y9[i] * 825.0 + y10[i] * 1200.0 + y12[i] * 4100) / 4100.0;
        }
        t1 = *ptr_t + *ptr_dt;
	//
        yhc(t1, x1, y13);

        //////////
        for (int i = 0; i < n; i++) 
	{
            x[i] = x[i] + *ptr_dt * (y6[i] * 272.0 + (y7[i] + y8[i]) * 216.0 + (y9[i] + y10[i]) * 27.0 + (y12[i] + y13[i]) * 41.0) / 840.0;
        }

        /////////
	// find the largest error
        err_now = 0;
        //
        for (int i = 0; i < n; i++)
	{
            double diff;
            diff = (y1[i] + y11[i] - y12[i] - y13[i]) * (*ptr_dt) * 41.0 / 810.0;
            err_now = fmax(err_now, fabs(diff));
        }

        *ptr_t = *ptr_t + *ptr_dt;

        if (err_now > err_tol) 
	{
            *ptr_dt = *ptr_dt / 2.0;
        }

        if (err_now < err_tol * err_ratio)
	{
            *ptr_dt = *ptr_dt * 2.0;
        }
	
	//
        /////////////////////
        // free the allocated arrays
        free(x1);
	x1 = NULL;
        free(y1);
	y1 = NULL;
        free(y2);
	y2 = NULL;
        free(y3);
	y3 = NULL;
        free(y4);
	y4 = NULL;
        free(y5);
	y5 = NULL;
        free(y6);
	y6 = NULL;
        free(y7);
	y7 = NULL;
        free(y8);
	y8 = NULL;
        free(y9);
	y9 = NULL;
        free(y10);
	y10 = NULL;
        free(y11);
	y11 = NULL;
        free(y12);
	y12 = NULL;
        free(y13);
	y13 = NULL;
    }

    free(x0);
    x0 = NULL;

    return 0;
}




