#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>
#include "nbody.h"
#include "alloc.h"

#include <time.h>

//int nbody1(double para0, double para1, double para2, double para3, double para4, double para5, double para6, int nline_data, double *tgrid_days, double **simul_ra_dec);
int nbody2(double para0, double para1, double para2, double para3, double para4, double para5, double para6, double para7, double para8, double para9, double para10, double para11, double para12, double para13, int nline_data, double* tgrid_days, double** simul_ra_dec);

const double tpi = 2.0 * 3.1415926535897932;
const int nbod = 3; // 天体数量
const int ndim = 6; // 积分维度

double xmass[MAX_BODIES]; // 质量（太阳质量单位）
// 新增：观测距离参数（单位：秒差距 pc）
const double observer_distance_pc = 5.0;

// RKF78积分器
int rkf78_no_t(double* ptr_t, double* ptr_dt, double err_tol, double* x, int n);
//int rkf78(double* ptr_t, double* ptr_dt, double err_tol, double* x, int n);

// 角度转换
double zatan(double sina, double cosa);

// 轨道要素转直角坐标
int ele2rat(double a, double e, double ci, double w, double omega, double cm0, int ii, double xout[6]);

// 直角坐标转轨道要素
int rat2ele(double xin[6], double* a, double* e, double* ci, double* w, double* omega, double* cm, int ii, int ito);



//int nbody(double para0, double para1, double para2, double para3, double para4, double para5, double para6, int nline_data, double *tgrid_days, double **simul_ra_dec)
int nbody2(double para0, double para1, double para2, double para3, double para4, double para5, double para6, double para7, double para8, double para9, double para10, double para11, double para12, double para13, int nline_data, double* tgrid_days, double** simul_ra_dec)

{
    int i, j;

    // 恒星质量（1.0太阳质量）
    xmass[0] = 1.022;
    const double yr_days = 365.25;

    //clock_t start, end;
    //double cpu_time_used;
    //start = clock();

    /* 读取初始轨道要素和质量 */
    double p0[MAX_BODIES], e0[MAX_BODIES], ci0[MAX_BODIES], w0[MAX_BODIES];
    double omg0[MAX_BODIES], cm00[MAX_BODIES];
    double raw_mass[MAX_BODIES];     /* 临时存放质量 */

    i = 1;
    p0[i] = para0;
    e0[i] = para1;
    ci0[i] = para2;
    w0[i] = para3;
    omg0[i] = para4;
    cm00[i] = para5;
    raw_mass[i] = para6;
    xmass[i] = raw_mass[i] * 3e-6;

    i = 2;
    p0[i] = para7;
    e0[i] = para8;
    ci0[i] = para9;
    w0[i] = para10;
    omg0[i] = para11;
    cm00[i] = para12;
    raw_mass[i] = para13;
    xmass[i] = raw_mass[i] * 3e-6;
    
    //FILE* faei = fopen("c_aei.txt", "w");
    //FILE* fd = fopen("c_d.txt", "w");
    // 打开赤经赤纬输出文件（单位：微角秒）
    //FILE* fra_dec = fopen("ra_dec.txt", "w");
    //if (!fra_dec) { perror("ra_dec.txt"); return 1; }
    //if (!fd) { perror("c_d.txt"); return 1; }

    /* 从周期计算半长轴，转换角度单位 */
    double a0[MAX_BODIES], a[MAX_BODIES];
    double e[MAX_BODIES], ci[MAX_BODIES], w[MAX_BODIES], omg[MAX_BODIES], cm0[MAX_BODIES], cm[MAX_BODIES];
    for (i = 1; i < nbod; i++) {
        a0[i] = pow(pow(p0[i] / 365.2422, 2) * xmass[0], 1.0 / 3.0);
        a[i] = a0[i];
        e[i] = e0[i];
        ci[i] = ci0[i] * tpi / 360.0; 
        w[i] = w0[i] * tpi / 360.0;
        omg[i] = omg0[i] * tpi / 360.0;
        cm0[i] = cm00[i] * tpi / 360.0;
    }

    /* 输出初始参数 
    printf("Initial a, e, mass (bodies 2..%d):\n", nbod);
    for (i = 1; i < nbod; i++) {
        printf("% .6f  % .6f  % .6f\n", a[i], e[i], xmass[i]);
    }
     */

    /* 状态向量初始化 */
    double x[MAX_BODIES][6] = { {0} };  // 位置和速度
    double xp[6];                       // 临时存储单个天体的坐标
    double xrk[MAX_BODIES * 6];         // 积分器用的1D数组


    // 初始化直角坐标
    for (i = 1; i < nbod; i++) {
        ele2rat(a[i], e[i], ci[i], w[i], omg[i], cm0[i], i, xp);
        for (j = 0; j < ndim; j++) x[i][j] = xp[j];
    }

    // 质心校准
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

    // 调整坐标（相对质心）
    for (i = 1; i < nbod; i++)
    {
        for (j = 0; j < ndim; j++)
        {
            x[i][j] += x[0][j];
        }
    }

    // 转换为1D数组供积分器使用
    for (i = 0; i < nbod; i++)
    {
        for (j = 0; j < ndim; j++)
        {
            xrk[i * ndim + j] = x[i][j];
        }
    }

    // 积分控制参数
    double h = tpi / yr_days / 1.0;    // 初始步长
    double err = 1e-6;                 // 误差 tolerance
    double t = 0.0;                     // 初始时间
    double tstop;
    tstop = (tgrid_days[nline_data-1])/ yr_days * tpi;       // 总积分时间

    //printf("%.6f\n", tstop);
    
    int output_next_step = 1;

    int i_out_ingrid = 0;
    double tout = tgrid_days[i_out_ingrid]/ yr_days * tpi;       // 首次输出时间

    double dist_au = observer_distance_pc * 206265.0; /* 1 pc = 206265 AU */
    double pi = tpi / 2.0;
    double rad2uas = (180.0 * 3600.0 * 1e6) / pi; /* 1 rad -> micro-arcsec */

    /* 主积分循环 */
    while (t < tstop) {
        // 更新坐标
        for (i = 0; i < nbod; i++)
        {
            for (j = 0; j < ndim; j++)
            {
                x[i][j] = xrk[i * ndim + j];
            }
        }

        /* 转换为轨道要素
        for (i = 1; i < nbod; i++)
        {
            for (j = 0; j < ndim; j++)
            {
                xp[j] = x[i][j] - x[0][j];
            }
            rat2ele(xp, &a[i], &e[i], &ci[i], &w[i], &omg[i], &cm[i], i, 0);
        }
	*/

        // 输出数据
        if ( (output_next_step == 1) || (i_out_ingrid == nline_data-1) )
	{

            double t_days = (t / tpi) * yr_days; /* t -> days */
            //printf("%.6f %lf", t, tout);

            /* 计算主星的小角近似 RA/Dec（相对于观测者，单位 µas） */
            double sx = x[0][0];   /* 主星 X (AU) */
            double sy = x[0][1];   /* 主星 Y (AU) */
            double sz = x[0][2];   /* 主星 Z (AU) */

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
                tout = tgrid_days[i_out_ingrid]/ yr_days * tpi;   // next输出时间
                //printf("%d %lf\n", i_out_ingrid, tgrid_days[i_out_ingrid]);
            }
        }



        if (t+h >= tout)  // the integrator will past the next time grid in the next step
	{
	    h = tout - t;
	    //
            // 自适应RK步长积分 
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

        /* 稳定性检查 */
        for (i = 1; i < nbod; i++) {
            if (fabs(a[i] - a0[i]) > 0.1 * a0[i]) { break; }
        }
    }

    /* 计算并输出运行时间
    end = clock();
    cpu_time_used = ((double)(end - start)) / CLOCKS_PER_SEC;
    printf("Program run time: %f s\n", cpu_time_used);
    */

    // 关闭文件
    //fclose(faei);
    //fclose(fd);
    //fclose(fra_dec);

    return 0;
}

