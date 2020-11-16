// This is a modification of original pdm2b_4.41.c and 
// this code is used to generate a python wrapper with Cython. 
// - I have removed all the file inputs and outputs of the pdm2 function, 
// and created "f_array" and "theta_array" as outputs; 
// - I have added "f_min", "f_max", "del_f", and "nbins" as inputs to 
// adjust the range and precision of frequency search.

/*  pdm2.c - period analysis package  */
/* $Id: pdm2.c,v 1.04 2009/12/08 01:19:42 rfs Exp rfs $  */

/*============================================================*
/*    S_TRAN - Copyright (c) 2012, Stellingwerf Consulting    *
This program is free software: you can redistribute it and/or modify
    it under the terms of the GNU General Public License as published by
    the Free Software Foundation, either version 3 of the License, or
    (at your option) any later version.
/*============================================================*/

/*    Version 4.14 - 4/12/2013  -  rich data file version      */

/*    4.14 - change MAXIT to 200          */
/*    4.13 - add check on sigma           */
/*    cleaned up a bit from version 4.10  */

/*    includes Blazhko Analysis coding                              */
/*    will do a series of subsegments of the data of length Nrange  */
/*        and slide them Nstep points to the right each pass        */
/*    summary of amplitudes and periods are in Blaz_sum.dat         */

/*  phase_shift added to shift pdmcurve plots  */

/*   3.04 changes   */
/*  this version reads data from an external text file     */
/*  data file may have 2 or 3 columns: time, data, [sigma] */
/*  separator must be either spaces or commas              */ 

/*  see the S-tran Application Guide and the PDM technical guide   */
/*  for details on setting up a specific problem                   */

/*-------------------100 bin version------------------------------------*/
/*  if less than 1000 pts in a segment, make the bins double wide (50/2) */
/*  bin centers are at 0, 0.1, 0.2, ...= ibin/10.                      */

/* this version has the oringinal PDM settings as defaults  */
/* updated confidence levels based on Beta distribution     */
/* NEW options are: 1) linear mean curve, 2) subharmaonic averages  */
/* META analyses include 1) search for changing period, and           */
/*  2) do a Monte-Carlo nalysis to get the confidence levels correct  */
/* SPLINE fit added 2/1/2011                                  */
/* Converted to read an external data file 2/2/2011  v.04 */
/* d^2/dx^2 corrections for spline fit modified to handle sharp kinks better 2/5/11 */
/*  data precision in dump files increased 12/01/11   v.06         */
/*  add missing data treatment             12/12/11   v.06         */

/* Modifications by Lucas Macri, June-July 2012 */
/* user can specify root name for input data file & output files */
/* user can specify a single trial period or a range of periods */
/* best period appended to pdm2b.out */
/* increased # of bins to 100 to analyze CSTAR data */

/*  SEE - PDM2.doc - user's guide   */ 

#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <stdio.h>
#include <ctype.h>

/*--------------------------------------------------------------------------*/

#define TRUE  1
#define FALSE 0
#define MISSING -99999

#define LOGICAL    int

#define MAXDATP  100000                 /*  max data points  */
#define LPOINTS  10                     /*  frequency points to cover line  */
#define SC_MAX   1000                   /*  switch point to single cover    */  
#define MAXF     1000000                /*  max allowable points in scan    */
#define MAXSEGS  100                    /*  max allowable segments  */
#define MAXBINS  100                    /*  max allowable bins  */ 
#define THMAX    500                    /*  bins in theta distribution   */

#define sq(x)         ((x)*(x))
#define cub(x)        ((x)*(x)*(x))
#define fmax(x1,x2)   ((x1) > (x2) ? (x1) : (x2))
#define fmin(x1,x2)   ((x1) < (x2) ? (x1) : (x2))
#define x(i)          datx[i]
#define y(i)          daty[i]

#define MAXIT 200
#define EPS 1.2e-7
#define FPMIN 1.0e-30

/* old form, based on the F test                                                      */
/* #define signif2(th)    (370. * pow( th, 9.34 )/pow( 109.+(1.-th)*(ne-109), 1.26 )) */

/*  this form is the exact analytical result. Tested, looks ok              */
/*    reference: Schwarzenberg-Czerny, 1997, Ap.J. 489, 941-945, eq 11      */
#define signif2(th)    inc_beta( (ne-m0)/2., (m0-1)/2., (ne-m0)*th/(ne-1) )

/*---local data---*/

static int nseg, debug, m0, bign = SC_MAX;
static int seg_beg_i[MAXSEGS+1], seg_end_i[MAXSEGS+1], seg_pts[MAXSEGS+1];
static double sum_all, data_mean, av_sig, sig_var, beta;
static double dt_avg, dt_min, segmn[MAXSEGS+1], seg_var[MAXSEGS+1], seg_var0[MAXSEGS+1];
// static double delf, dt_avg, dt_min, segmn[MAXSEGS+1], seg_var[MAXSEGS+1], seg_var0[MAXSEGS+1];
static double sumy2, sig02, sig2, theta, rscale;
static double daty0[MAXDATP+1];
static long dummy;           /*  for ran1  */

/*--------------------------- prototypes ------------------------------------*/

int pdm2( int ne, double datx[], double daty[], double sig[], double f_min, double f_max, double delf, int nbins);
int p_sort( int n, double dat1[], double dat2[], double dat3[], int sgn );
int i_sort( int n, int dat1[], double dat2[], int sgn );
double inc_beta( double a, double b, double x );
double rnd( void );
double ran1(long *idum);

static int binner( double datx[], double daty[], double sig[], double f, int bin10, int first, int last, int nbins);
static double dophase( double tt, double t0, double tn, double f );
static int segset( int ne, double datx[], double daty[], double segdev );
static double sqrc( double x );
// static void error( char *str );
static double mcurve( double phase );
double bspl3( double z );
static double table_interp( double x0, int n, double xt[], double yt[] );
// static int read_data_file( char *file, int title_lines );
// static int fgetline(FILE *fp, char *buffer, int max_len);


/*----these globals are provided for optional external control of the process----*/

/*----run params - can be set externally----*/
/*    if zero, default values will be used  */
int invert_curve = FALSE;                  /*  plot negative of curve                 */
// int lpoints = 10;                          /*  number of points to cover line         */
double minf0, maxf0;                       /*  freq scan range                        */
double segdev = 200;                       /*  sensitivity for segments (big->1 seg)  */
double beta_min = -5, beta_max = 5;        /*  period change mode                     */
double beta_scale = 1./365.25e6;           /*  scale factor for beta                  */
double phase_shift = 0.;                   /*  shift phase of pdmcurve plot           */

int nb0 = 3;                               /*  points in beta scan                    */
int do_beta_scan = FALSE;                  /*  set to true for beta scan              */
int do_subharm = TRUE;                     /*  subharmonic averaging                  */
int do_dist = 0;                           /*  write theta distribution files         */
int bin_10 =  2;                           /*  0=5/2 bins,1=10/1,2=auto,>2=sw_#/auto  */
int do_linear_fit = FALSE;                 /*  linear curve fitting                   */
int do_spline_fit = TRUE;                  /*  Bspline curve fitting                  */
int pdm_verbose = 0;                       /*  generate screen output                 */
int pdm_debug = 0;                         /*  turn on local debugging                */
int do_non_par = FALSE;                    /*  select nonparametric sig test          */
int nb1 = 250;                             /*  points in Monte Carlo analysis         */
int do_sigmas = TRUE;                      /*  use sigmas in computation              */
char root[31] = "pdm";                     /* run prefix                              */

int Blazhko = FALSE;                       /* Blazhko flag                            */
int Nstep = 250;                           /* step length (points) for Blazhko        */
int Nrange = 1000;                         /* segment length for Blazhko              */

/*---results - can be used externally----*/
double trange;                                    /* range of time                    */
double fthmin[4], thmin[4], signf[4];             /* freq,theta,signif at 3 minima    */
double bin_mean[MAXBINS+1], bin_var[MAXBINS+1];   /* bin data, mean curve             */
int nbin[MAXBINS+1], nf;                          /* points per bin / scan            */ 
int theta_dist[THMAX+1], tot_points;              /* numerical theta distribution     */
double dtheta[THMAX+1];                           /* theta values for dist            */
double theta_dist2[THMAX+1];                      /* theta_min distribution           */
int tot_points2;                                  /* npoints for theta_min distr      */
// double f_min, f_max
double theta2[MAXDATP+1];           /* final theta scan result          */
int ran_array[MAXDATP+1];                         /* for Nemec significance test      */
double ratio;                                     /*  1/(S/N ratio) from sigmas       */
int nplot;                                        /* number of points in plot file    */
double ymin, ymax, yamp, ymean_mean, xmean;       /* mean curve parameters            */

/*---data arrays---*/
// int Numdat;
double times[MAXDATP+1];
double mags[MAXDATP+1];
double sigs[MAXDATP+1];

/*---output arrays---*/
double* f_array;
double* theta_array;

/*---number of bins---*/
// int nbins=100; // use this number as the number of bins, it should be lower than the MAXBINS

int main(){
    // double x[350], y[350], s[350];
    int i;
    printf("hello world\n");

    for(i=0; i<350; ++i){
        times[i] = 20. / 350. * i;
        mags[i] = sin(times[i]);
        sigs[i] = 0.0;
    }

    pdm2(350, times, mags, sigs, 0.01, 1, 0.001, 100);
    printf("nf is %d", nf);
    for(i=0; i<nf; ++i){
        printf(" f: %f, theta: %f", f_array[i], theta_array[i]);
    }

    // free(f_array);
    // free(theta_array);
    return 0;
}

int pdm2(int ne, double datx[], double daty[], double sig[], double f_min, double f_max, double delf, int nbins){
    int i, j, k, kk, bins, sub1_index, big_seg, big_pts, bin10, pts, ntot, pdv;
    // int line_points;
    int icurr, nb;
    double s2, s20, seg_x0[MAXSEGS+1], seg_xrange[MAXSEGS+1], sum, sum1, sum2, f, f1, th1, th2, phase, fslope;
    double tmp, tmp2, ymean, bin_mean0[MAXBINS+1], resid, sig0, signf0;
    double th0min[4], delbeta, beta0, ndof, theta_crit, hifact;
    // char *stp;
    // char file[31];
    // FILE *fo, *fp, *fp2, *fb;
    
    // printf("Received root as %s from main\n",root);
    
    nf = (int)((f_max - f_min)/delf + 1);

    f_array = realloc(f_array, nf * sizeof(double));
    theta_array = realloc(theta_array, nf * sizeof(double));

    // f_array = f_array_larger;
    // theta_array = theta_array_larger;

    if( ne <= 100 )  return -1; //pdm: too few points
    // if( pdm_debug )  debug = TRUE;
    
    // fo = stdout; 
    
    /*  estimate of significant theta level  */
    theta_crit = (bin_10 ? 1. - 11./pow((double)ne,0.8) : 1. - 8/pow((double)ne,0.8));
    
    /*  theta array for distributions  */
    // for( i = 0; i <= THMAX; i++ )  dtheta[i] = 1.2*i/THMAX;
    
    /*----------beta scan option---(changing period)-----------*/
    
    // nb = nb0;
    
    // if( do_beta_scan ) {
    //     pdm_verbose = FALSE;
    //     fprintf( fo, "\n Beta     Thetamin    F\n" );
    //     sprintf(file,"%s_beta.dat",root);
    //     fb = fopen(file, "w" );
    //     if( !fb )  error( "Could not open beta_scan.dat" );
    //     delbeta = (beta_max - beta_min) / (nb-1);
    //     if( !delbeta )  error( "Zero range of beta" );
    //     fprintf( fb, "#Beta Thetamin F\n" );
    // }
    // else {
    //     nb = 1;
    // }
    nb = 1;
    /*--------------------monte-carlo sig option----------------------------*/
    
    // if( do_non_par ) {
    //     printf("\n==========MONTE CARLO SIGNIFICANCE ANALYSIS============\n" );
        
    //     nb = nb1;
    //     for( i = 1; i <= ne; i++ ) {
    //         daty0[i] = daty[i];
    //     }
        
    //     for( i = 0; i <= THMAX; i++ )  theta_dist2[i] = 0.;
    //     tot_points2 = 0;
    //     pdv = pdm_verbose;
    //     pdm_verbose = FALSE;
    // }
    
    /*=============outer loop for beta or non-par runs=================*/
    
    for( kk = 1; kk <= nb; kk++ ) {
        // if( do_beta_scan ) {
        //     beta0 = beta_min - delbeta + kk * delbeta;
        //     beta = beta0 * beta_scale;  
        // }
        
        // else if( do_non_par ) {
        //     /*  set to original order  */
        //     for( i = 1; i <= ne; i++ ) {
        //         daty[i] = daty0[i];
        //     }
        //     /*  final pass  */
        //     if( kk == nb ) {
        //         pdm_verbose = pdv;
        //     }
        //     else {
        //         /*  scramble the data set  */
        //         for( i = 1; i <= ne; i++ ) {
        //             ran_array[i] = (int)(1.e8 * rnd());
        //         }
        //         i_sort( ne, ran_array, daty, 1 );
        //         printf( "   iteration %d / %d\r", kk, nb );
        //     }
        // }
        
        /*--------------------------------------------------------*/
        
        /* set bin structure  */
        bin10 = bin_10;
        if( bin_10 > 2 )  bign = bin_10;
        
        p_sort( ne, datx, daty, sig, 1 );           /*  sort the data  */
        
        // if( lpoints )  line_points = lpoints;
        // else           line_points = LPOINTS;
        
        /*  set segments  */ 
        if( !segdev )  segdev = 2.;
        if( !segset( ne, datx, daty, segdev ) ) return(0);
        
        // f_min = f_max = 0;
        
        /*  write data file  */
        // sprintf(file,"%s_data.dat",root);
        // fp = fopen(file, "w" );
        
        // if( invert_curve ) {
        //     fprintf( fp, "#Time Val Sig\n" );
        //     for( i = 1; i <= ne; i++ ) {
        //         if( daty[i] == MISSING ) {
        //             fprintf( fp, "%18.12g %18.12g %18.12g\n", datx[i], 99.9999999,fabs(sig[i]) );
        //         }
        //         else {
        //             fprintf( fp, "%18.12g %18.12g %18.12g\n", datx[i], -daty[i], fabs(sig[i]) );
        //         }
        //     }
        // }
        // else {
        //     fprintf( fp, "#Time Val Sig\n" );
        //     for( i = 1; i <= ne; i++ ) {
        //         if( daty[i] == MISSING ) {
        //             fprintf( fp, "%18.12g %18.12g %18.12g\n", datx[i], 99.9999999,fabs(sig[i]) );
        //         }
        //         else {
        //             fprintf( fp, "%18.12g %18.12g %18.12g\n", datx[i], daty[i], fabs(sig[i]) );
        //         }
        //     }
        // }
        // fclose( fp );
        
        /*  sigmas  */
        if( do_sigmas ) {
            sig_var = av_sig = 0.;
            for( i = 1, sum = 0, s2 = 0l; i <= ne; i++ ) {
                sig[i] = fabs( sig[i] );
                sig_var += sq( sig[i] ) / ne;
                av_sig += sig[i] / ne;
                sum += y(i);
                s2 += sq( y(i) );
            }
            sig2 = (s2 - sq( sum ) / ne ) / (ne - 1.);
            
            /*  sigma report  */
            if( sig2 )  ratio = sqrt(sig_var / sig2);
            if( ratio ) {
                rscale = fmin( 1., 0.2 / ratio );
            }
            else {
                rscale = 1.;
            }
            // if( sig_var && pdm_verbose )  fprintf( fo, "Variance=%.4f, Sig_var=%.4f, Rsc=%.4f\n", sig2, sig_var, rscale );
        }
        
        /*  print header  */ 
        // if( pdm_verbose ) {
        //     fprintf( fo, "\n*** PDM2 PERIOD ANALYSIS ***\n" );
        //     fprintf( fo, "\nN = %d     DOF = %d, %d\n", ne, ne - 1, ne - 10*nseg );
        // }
        
        /*  full data statistics  */
        seg_beg_i[nseg+1] = ne + 1;
        sum_all = s2 = s20 = trange = 0.;
        big_seg = 1;
        big_pts = 0;
        ntot = 0;
        for( i = 1; i <= nseg; i++ ) {
            seg_end_i[i] = seg_beg_i[i+1] - 1;
            seg_pts[i] = seg_beg_i[i+1] - seg_beg_i[i];
            if( seg_pts[i] > big_pts ) {
                big_seg = i;
                big_pts = seg_pts[i];
            }
            seg_x0[i] = x(seg_beg_i[i]);
            seg_xrange[i] = x(seg_end_i[i]) - x(seg_beg_i[i]);
            trange = fmax( trange, seg_xrange[i] );
            /*  loop to get means  */
            for( j = seg_beg_i[i], sum = 0., pts = 0; j <= seg_end_i[i]; j++ ) {
                /*  skip missing  */
                if( y(j) == MISSING )  continue;
                /*  skip outliers  */
                if( do_sigmas && sq( rscale*sig[j] ) > 3. * sig_var )  continue;
                sum_all += y(j);
                sum += y(j);
                pts++;
            }
            segmn[i] = sum / pts;
            ntot += pts;
            
            /*  full data variance, w/ sig correction  */
            for( j = seg_beg_i[i], sum1 = 0., sum2 = 0.; j <= seg_end_i[i]; j++ ) {
                /*  skip missing  */
                if( y(j) == MISSING )  continue;
                /*  skip outliers  */
                if( do_sigmas && sq( rscale*sig[j] ) > 3. * sig_var )  continue;
                tmp = fabs( y(j) - segmn[i] );
                sum1 += sq( tmp );
                
                if( do_sigmas ) tmp2 = rscale*sig[j];
                else            tmp2 = 0.;
                
                tmp = fmax(tmp - tmp2, 0);
                sum2 += sq( tmp );
            }
            if( seg_pts[i] > 1 ) {
                seg_var0[i] = sum1 / (seg_pts[i]-1.);
                seg_var[i] = sum2 / (seg_pts[i]-1.);
            }
            else {
                seg_var0[i] = 0.;
                seg_var[i] = 0.;
            }
            s20 += (seg_pts[i] - 1.) * seg_var0[i];
            s2 += (seg_pts[i] - 1.) * seg_var[i];
        }
        data_mean = sum_all / ntot;
        sig02 = s20 / (ntot - nseg);          /*  sig02 is the normal total variance  */
        sig2 = s2 / (ntot - nseg);            /*  sig2  is corrected for data sigmas  */
        
        /*  segment report  */
        // if( pdm_verbose ) {
        //     fprintf( fo, "\nSeg:  #    St     N     Tstart    Trange      Mean      s.d.     Bins\n" );
        //     for( i = 1; i <= nseg; i++ ) {
        //         if( bin10 >= 2 )  stp = ( seg_pts[i] > bign ) ? "100/1" : "50/2";
        //         else stp = bin10 ? "100/1" : "50/2";
        //         fprintf( fo, "    %3d   %3d    %3d  %8.2f  %8.2f  %8.2f  %8.2f     %s\n",
        //             i, seg_beg_i[i], seg_pts[i], seg_x0[i], seg_xrange[i], segmn[i], sqrt(seg_var[i]), stp );
        //     }
        //     fprintf( fo, "\nStandard Dev = %g, Variance = %g, Trange = %g, DTavg=%g\n", sqrc(sig2), sig2, trange, dt_avg );
        // }
        
        /*  pick scan range  */
        // delf = 1./(line_points * trange);
        // f_min = delf;
        // f_max = 1. / (2. * dt_avg); 
        /*  check and use input params  */
        // if( minf0 )  f_min = fmax( delf, minf0 );
        // if( maxf0 )  f_max = maxf0;
        
        nf = (int)((f_max - f_min)/delf + 1);
        hifact = f_max * 2. * trange / ne;
        
        fslope = (nf-1.) / (f_max-f_min);
        
        if( seg_pts[nseg] > SC_MAX ) m0 = 100;
        else                         m0 = 50;
        
        /*---significance distribution tester------*/
        // if( do_dist ) {
        //     sprintf(file,"%s_sig.dat",root);
        //     fp = fopen(file, "w" );
        //     fprintf( fp, "#Theta Beta Beta_bc\n" );
        //     for( i = 1; i <= 121; i++ ) {
        //         tmp = (i-1)/100.;
        //         sig0 = signif2( tmp );
        //         /*  apply bandwidth correction  */
        //         //tmp2 = 2.*(float)nf/fmin(5,line_points);
        //         //tmp2 = 10*dt_avg/dt_min;
        //         tmp2 = ne * hifact;
        //         if( bin10 )  tmp2 *= 2.;          /*  from monte carlo results  */
        //         signf0 = 1. - pow( 1-sig0, tmp2); 
        //         fprintf( fp, "%g %g %g\n", tmp, sig0, signf0 );
        //     }
        //     fclose( fp );
        //     if( pdm_verbose )  printf( "Beta distribution written to sig.dat\n" );
        // }
        
        // sprintf(file,"%s_pdmplot.dat",root);
        // fp = fopen(file, "w" );
        // fprintf( fp, "#Frequency Theta\n" );
        
        /*-----------------------------------------*/
        
        // if( pdm_verbose ) {
        //     fprintf( fo, "#f_min = %g, f_max = %g, delf = %g, nf = %d\n", f_min, f_max, delf, nf);
            
            
        //     if( nf > MAXF ) {
        //         printf( "     *** too many points, set to %d ***", MAXF );
        //         nf = MAXF;
        //         delf = (f_max - f_min) / (nf - 1);
        //     }
        //     if( nf < 1 ) {
        //         printf( ".*** frequency points=  %d abort scan ***", nf );
        //         return(0);
        //     }
            
        //     fprintf( fo, "\nTHETA SCAN:  F range [%g -> %g]", f_min, f_max );
            
        //     fprintf( fo, "    %d points/line\n", (int)(1. / (trange * delf )) );
            
        //     printf( "        ***SCANNING***\n" );
        // }
        
        
        /*--------------------- compute ---------------------------------------*/
        
        if( !do_non_par ) {
            for( i = 0; i <= THMAX; i++ )  theta_dist[i] = 0;
            tot_points = 0;
        }
        
        icurr = 0;
        for( f = f_min, k = 0; f <= f_max + delf; f += delf ) {
            icurr++;
            // if( pdm_verbose )  printf( "===scan %d / %d frequency points===\r", icurr, nf );
            if( k ) {
                th2 = th1;
                th1 = theta;
            }
            s2 = 0.;
            bins = 0;
            ndof = 0.;
            for( i = 1; i <= nseg; i++ ) {
                if( bin_10 >= 2 )  bin10 = seg_pts[i] > bign;
                
                /* fprintf( fo, "\nSEG %d   bin10 = %d\n", i, bin10 ); */
                
                binner( datx, daty, sig, f, bin10, seg_beg_i[i], seg_end_i[i], nbins);
                
                for( j = 0; j <= nbins; j++ ) {
                    if( nbin[j] > 1 ) {        /*  bin statistics  */
                        bins++;
                        ndof += nbin[j]-1;
                        s2 += (nbin[j] - 1) * bin_var[j];
                    }
                }
            }
            /*=============doit!!===============*/
            k++;
            theta = s2 / (ndof * sig2);
            theta2[k] = theta;                      /*  unaveraged theta array  */
            if( k == 1 ) {
                th1 = th2 = theta;
                thmin[1] = thmin[2] = thmin[3] = theta;
                th0min[1] = th0min[2] = th0min[3] = theta;
                fthmin[1] = fthmin[2] = fthmin[3] = 1.;
            }
            
            /*  subharmonic averaging when possible */
            if( do_subharm ) {
                if( f/2. >= f_min &&  theta < theta_crit ) {
                    sub1_index = (int)(1 + fslope * (f/2 - f_min)+0.5);
                    theta = (theta + theta2[sub1_index])/2.;
                }
            }
            /*  accumulate for distribution  */
            // if( do_dist ) {
            //     i = (int)floor(theta*THMAX/1.2);
            //     i = fmin( i, THMAX );
            //     theta_dist[i]++;
            //     tot_points++;
            // }
            
            // if( debug ) fprintf( fo, "f=%g, s2=%g, ndof=%g, theta=%g\n", f, s2, ndof, theta);
            
            /* fprintf( fo, "f=%g  theta=%g\n", f, theta);  */
            // printf("f=%g  theta=%g\n", f, theta);
            // printf("%d", icurr);
            f_array[icurr-1] = f;
            theta_array[icurr-1] = theta;
            // fprintf( fp, "%g %g\n", f, theta );
            
            /*  save 3 minima  */
            // if( th1 < theta && th1 < th2 ) {
            //     f1 = f - delf;
            //     if( debug )  printf( "min, f, th = %14.6e, %14.6e\n", f1, th1 );
            //     if( th1 < thmin[3] ) {
            //         if( th1 > thmin[2] ) {
            //             thmin[3] = th0min[3] = th1;
            //             fthmin[3] = f1;
            //         }
            //         else if( th1 > thmin[1] ) {
            //             thmin[3] = thmin[2];
            //             th0min[3] = th0min[2];
            //             fthmin[3] = fthmin[2];
            //             thmin[2] = th0min[2] = th1;
            //             fthmin[2] = f1;
            //         }
            //         else {
            //             thmin[3] = thmin[2];
            //             th0min[3] = th0min[2];
            //             fthmin[3] = fthmin[2];
            //             thmin[2] = thmin[1];
            //             th0min[2] = th0min[1];
            //             fthmin[2] = fthmin[1];
            //             thmin[1] = th0min[1] = th1;
            //             fthmin[1] = f1;
            //         }
            //     }
            // }
        }
        // fclose( fp );
        
        /*  accumulate for extreme value distribution  */
        // if( do_non_par ) {
        //     i = (int)floor(thmin[1]*THMAX/1.2);
        //     i = fmin( i, THMAX );
        //     theta_dist2[i] += 1.;
        //     tot_points2++;
        // }
        
        /*--------------- write the final data sets -----------------------------*/
        
        // if (write == 1) {
        //     fp = fopen("pdm2b.out","a");
        //     fprintf (fp,"%-8s %10.7f\n",root,1./fthmin[1]);
        //     fclose(fp);
        // }
        // sprintf(file,"%s_pdmcurve.dat",root);
        // fp = fopen(file, "w" );
        // sprintf(file,"%s_residuals.dat",root);
        // fp2 = fopen(file, "w" );
        // if( invert_curve )  fprintf( fp, "#Phase(F=%.5f)(P=%.7f)) Val- Mean- Sigma Num\n", fthmin[1], 1./fthmin[1] );
        // else  fprintf( fp, "#Phase(F=%.5f)(P=%.7f)) Val Mean Sigma Num\n", fthmin[1], 1./fthmin[1] );
        
        // fprintf( fp2, "#Resid X\n" );
        
        // binner( datx, daty, sig, fthmin[1], bin10, seg_beg_i[big_seg], seg_end_i[big_seg], nbins);
        
        // if( do_spline_fit && thmin[1] < theta_crit ) {
        //     /*  apply curvature corrections  - use for final mean curve only  */
        //     /*  will show on mean curve plot and improve residuals            */
        //     /*    equivalent to a converged Stobie iteration                  */
        //     for( j = 0; j <= nbins; j++ ) {
        //         bin_mean0[j] = bin_mean[j];
        //     }
        //     for( j = 1; j <= nbins-1; j++ ) {
        //         bin_mean[j] = (3.*bin_mean0[j] - 0.5*(bin_mean[j-1] + bin_mean[j+1])) / 2.;
        //     }
        //     bin_mean[0] = (3.*bin_mean0[0] - 0.5*(bin_mean[nbins-1] + bin_mean[1])) / 2.;
        //     bin_mean[nbins] = bin_mean[0];
        // }       
        /*------------- write data versus phase and residuals ------------------------*/
        /*  do for largest segment only to get meaningful data    */
        
        // ymean_mean = 0.;
        // for( i = seg_beg_i[big_seg]; i <= seg_end_i[big_seg]; i++ ) {
        //     if( y(i) == MISSING )  continue;
            
        //     phase = dophase( x(i), x(seg_beg_i[big_seg]), x(seg_end_i[big_seg]), fthmin[1] );

        //     /*  compute mean curve  */
        //     ymean = mcurve( phase );
        //     resid = y(i) - ymean;

        //     /* plot data at given phases  */
        //     if( invert_curve ) {
        //         tmp = -y(i);
        //         ymean *= -1.;
        //     }
        //     else  tmp = y(i);
            
		// 	/*  mean curve parameters for Blazhko file  */
        //     if( i == seg_beg_i[big_seg] )  ymin = ymax = ymean;
        //     ymax = fmax( ymax, ymean);
        //     ymin = fmin( ymin, ymean );
        //     ymean_mean += ymean;
            
        //     tmp2 = sig[i];
            
        //     // fprintf( fp, "%g %18.12g %g %g %d\n", phase, tmp, ymean, tmp2, i );
        //     // fprintf( fp2, "%18.12g %g\n", x(i), resid );
        // }
        // /*  parameters for mean curve  */
        // yamp = ymax - ymin;
        // xmean = (x(seg_beg_i[big_seg])+x(seg_end_i[big_seg])) / 2.;
        // nplot = seg_end_i[big_seg] - seg_beg_i[big_seg] + 1;
        // ymean_mean /= nplot;
        
        // fclose( fp );
        // if( pdm_verbose )  printf( "Phased data and mean curve written to pdmcurve.dat\n" );
        // fclose( fp2 );
        // if( pdm_verbose )  printf( "Residuals written to residuals.dat\n" );
        
        /*  distribution plot  */
        // if( do_dist ) {
        //     sprintf(file,"%s_theta_dist.dat",root);
        //     fp = fopen(file, "w" );
        //     if( !fp )  error( "Could not open theta_dist" );
        //     fprintf( fp, "#Theta Dist Dist_min\n" );
        //     tmp = 0.;
        //     tmp2 = 0.;
        //     for(i = 0; i <= THMAX; i++ ) {
        //         tmp += (double)theta_dist[i]/tot_points;
        //         if( do_non_par )  tmp2 += theta_dist2[i]/tot_points2;
        //         fprintf( fp, "%g %g %g\n", dtheta[i], tmp, tmp2 );
        //     }
        //     fclose( fp );
        //     if( pdm_verbose )  printf( "Theta distribution written to theta_dist.dat\n" );
        // }
        
        /*=======  compute significances here  =======*/
        // if( do_non_par && kk == nb ) {
        //     tmp2 = 0.;
        //     for(i = 0; i <= THMAX; i++ ) {
        //         tmp2 += theta_dist2[i]/tot_points2;
        //         theta_dist2[i] = tmp2;
        //     }
        //     for( i = 1; i <= 3; i++ ) {
        //         signf[i] = table_interp( th0min[i], THMAX, dtheta, theta_dist2 );
        //     }
        //     // if( pdm_verbose )  fprintf( fo, "\nMinima: #    Theta     Frequency      Period    Signif(MC)\n" );
        // }
        // else {
        //     for( i = 1; i <= 3; i++ ) {
        //         sig0 = signif2( th0min[i] );
        //         /*  apply bandwidth correction  */
        //         //tmp = 2.*(float)nf/fmin(5,line_points);
        //         //tmp = 10*dt_avg/dt_min;
        //         //tmp2 = ne;
        //         tmp = ne * hifact;
        //         if( bin10 )  tmp *= 2.;          /*  from monte carlo results  */
        //         signf[i] = 1. - pow( 1. - sig0, tmp); 
        //         if( thmin[i] <= 0. )  thmin[i] = 0.;
        //     }
        //     // if( pdm_verbose )  fprintf( fo, "\nMinima: #    Theta     Frequency      Period     Signif(beta)\n" );
        // }
        
        /*=======  screen summary  =======*/
        // if( pdm_verbose ) {
        //     //fprintf( fo, "\nMinima: #    Theta     Frequency      Period       Signif\n" );
        //     for( i = 1; i <= 3; i++ ) {
        //         fprintf( fo, "        %d   %6.6f   %10.8f   %10.8f   %6.6f\n",
        //             i, thmin[i], fthmin[i], 1./fthmin[i], signf[i] );
        //     }
        // }
        
        // if( do_beta_scan) {
        //     fprintf( fo, " %5.5f  %5.5f  %5.5f, \n", beta0, thmin[1], fthmin[1] );
        //     fprintf( fb, "%g %g %g\n", beta0, thmin[1], fthmin[1] );
        // }
        // if( fo != stdout ) {
        //     fclose(fo);
        // }
    }    // end of k block
    
    // if( do_beta_scan ) {
    //     fprintf( fo, "\nScan finished, data in beta_scan.dat\n" );
    //     fclose( fb );
    // }
    // if( do_non_par ) {
    //     printf( "\n %d data distributions analyzed\n   theta distributions are in theta_dist.dat\n", nb );
    // }

    // for(i=0; i<nf; ++i){
    //     printf(" f: %f, theta: %f", f_array[i], theta_array[i]);
    // }

    return(1);
}

/*------ compute value of mean curve fit to data at given phase ---------------------------*/
/*    pick frequency first and call binner()                */
/*    this is a linear interp between bin means             */
/*    optional: do a Bspline fit                            */

double mcurve( double phase )
{
    int i0, i1, i2, i3, j0, j1, j2, j3;
    double phase0, ymean;
    
    /*  i1 = left bin, i2 = right bin number   */
    /*  i's are for means, j's for the phases  */
    
    i1 = (int)floor(100.*phase)%100;
    if( i1 < 0 )  i1 += 100;
    i2 = i1 + 1;
    i0 = i1 - 1;
    if( i0 < 0 )  i0 += 100;
    i3 = i2 + 1;
    if( i3 > 100 )  i3 -= 100;
    
    j1 = i1;
    j0 = j1-1;
    j2 = j1+1;
    j3 = j2+1;
    
    /*  this takes care of the first 1/2 bin  */
    if( phase < 0.5 && i1 > 50 )  phase0 = phase + 1;
    else                          phase0 = phase;
    
    
    /*  4 splines contribute to each point inthe curve  */
    if( do_spline_fit ) {
        ymean = bin_mean[i1] * bspl3(100.*phase0 - j1) + bin_mean[i2] * bspl3(j2 - 100.*phase0);
        ymean += bin_mean[i0] * bspl3(100.*phase0 - j0) + bin_mean[i3] * bspl3(j3 - 100.*phase0);
    }
    /*  for linear fit, only two points contribute  */
    else if( do_linear_fit ) {
        ymean = bin_mean[i1] + 100. * (bin_mean[i2] - bin_mean[i1]) * (phase0 - i1/100.);
    }
    else {
        i1 = (int)(100.*phase+0.5) % 100;
        ymean = bin_mean[i1];
    }
    return( ymean );
}

/*  cubic B-Spline (sometimes called W4)  */
/*  unit width and unit area              */
/*  center at 0, zero beyond +/- 2        */
/*  peak value = 2/3 at origin            */

double bspl3( double z ) 
{
    double tmp;
    
    z = fabs( z );
    
    if( z < 1.0 ) {
        tmp = (2./3. - sq( z ) + cub( z )/2.);
    }
    else if( z < 2.0 ) {
        tmp = cub( 2. - z ) / 6.;
    }
    else {
        tmp = 0.;
    }
    //printf( "\nbspline, z = %g, tmp = %g\n", z, tmp );
    return( tmp );
}


/*------- compute bin sums and numbers ----------------------------------------------*/

int binner( double datx[], double daty[], double sig[], double f, int bin10, int first, int last, int nbins)
{
    double phase, t0, tn, sybin[MAXBINS+1], sumy2_dev[MAXBINS+1], ymean=0.; 
    double tmp, theta0, dev, dev0, bm1, bm2;
    int i, j, bin1, bin2, n0;
    int do_curve;
    
    /*  inits  */
    for( j = 0; j <= 99; j++ ) {
        nbin[j] = 0;
        sybin[j] = 0.;
        sumy2_dev[j] = 0.;
    }
    
    t0 = x(first);
    tn = x(last);
    
    /*  compute bin statistics  */
    for( j = first; j <= last; j++ ) {
        /*  skip missing values  */
        if( y(j) == MISSING )  continue;
        /*  skip outliers  */
        if( do_sigmas && sq( rscale*sig[j] ) > 3. * sig_var )  continue;
        
        phase = dophase( x(j), t0, tn, f );
        /*  (100,1) bins  */
        /* bins are 0->99, center of b0 at 0  */
        if( bin10 ) {
            bin1 = (int)(nbins*phase+0.5) % nbins;
            nbin[bin1] += 1;
            sybin[bin1] += y(j);
        }
        /*  (50,2) bins)  */
        else {
            // bin1 = 2 * ((int)(50.*phase) % 50) + 1;   /*  odd bins  */
            // bin2 = 2 * ((int)(50.*phase+0.5) % 50);   /*  even bins  */
            // bin1 = ceil(phase/(1./nbins)) - 1;
            // if (bin1&1) {
            //     bin2 = bin1 + 1;
            //     if(bin2==nbins) bin2=0;
            // }
            // else {
            //     bin2 = bin1;
            //     bin1--;
            //     if(bin1==-1) bin1=1;
            // }
            if (nbins&1){
                bin1 = 2 * ((int)(((nbins-1) / 2)*phase) % ((nbins-1) / 2)) + 1;
                bin2 = 2 * ((int)(((nbins-1) / 2)*phase+0.5) % ((nbins-1) / 2));
                if (bin1==nbins-2) bin2 = nbins - 1;
            }
            else{
                bin1 = 2 * ((int)((nbins / 2)*phase) % (nbins / 2)) + 1;
                bin2 = 2 * ((int)((nbins / 2)*phase+0.5) % (nbins / 2));
            }
            nbin[bin1] += 1;
            sybin[bin1] += y(j);
            nbin[bin2] += 1;
            sybin[bin2] += y(j);
        }
    }
    
    /*  compute bin means - bin 0 has at least 1 point  */
    for( j = 0; j <= nbins; j++ ) {
        if( nbin[j] )  bin_mean[j] = sybin[j] / nbin[j];
    }
    bin_mean[nbins] = bin_mean[0];
    /*  fix empty bins  */
    for( j = 1; j <= nbins; j++ ) {
        if( !nbin[j] ) {
            for( i = j+1; i <= nbins; i++ ) {
                if( nbin[i] ) {
                    bm1 = bin_mean[j-1];
                    bm2 = bin_mean[i];
                    bin_mean[j] = bm1 + (bm2-bm1)/(i-j+1);
                    /*     printf( "Bin %d empty, bin %d mean=%g, bin %d mean=%g, final=%g\n", j, j-1, bin_mean[j-1], i, bin_mean[i], bin_mean[j] );
                    getchar(); */
                    break;
                }
            }
        }
    }
    /*  well below noise result  */
    n0 = fmax( last - first, 2 );
    theta0 = (bin10 ? 1. - 11./pow((double)n0,0.8) : 1. - 8./pow((double)n0,0.8));
    
    do_curve = FALSE;
    /*  compute bin variances  */
    for( j = first; j <= last; j++ ) {
        /*  skip missing values  */
        if( y(j) == MISSING )  continue;
        
        phase = dophase( x(j), t0, tn, f );
        /*  compute mean curve at minima  */
        if( j > first && do_linear_fit  ) { 
            do_curve = TRUE;
            ymean = mcurve( phase );
        }
        /*  sig check  */
        if( do_sigmas ) tmp = rscale*sig[j];
        else            tmp = 0.;
        /*  (100,1) bins  */
        if( bin10 ) {
            bin1 = (int)(nbins*phase+0.5) % nbins;
            // bin1 = (int)(100.*phase+0.5) % 100;
            if( !do_curve )  ymean = bin_mean[bin1];
            dev0 = fabs(y(j) - ymean);
            /*  note sigma correction here  */
            dev = fmax( dev0 - tmp, 0. );
            sumy2_dev[bin1] += sq( dev  );
        }
        /*  (50,2) bins)  */
        else {
            // bin1 = 2 * ((int)(50.*phase) % 50) + 1;   /*  odd bins  */
            // bin2 = 2 * ((int)(50.*phase+0.5) % 50);   /*  even bins  */
            if (nbins&1){
                bin1 = 2 * ((int)(((nbins-1) / 2)*phase) % ((nbins-1) / 2)) + 1;
                bin2 = 2 * ((int)(((nbins-1) / 2)*phase+0.5) % ((nbins-1) / 2));
                if (bin1==nbins-2) bin2 = nbins - 1;
            }
            else{
                bin1 = 2 * ((int)((nbins / 2)*phase) % (nbins / 2)) + 1;
                bin2 = 2 * ((int)((nbins / 2)*phase+0.5) % (nbins / 2));
            }
            if( !do_curve )  ymean = bin_mean[bin2];
            dev0 = fabs(y(j) - ymean);
            /*  note sigma correction here  */
            dev = fmax( dev0 - tmp, 0. );
            sumy2_dev[bin2] += sq( dev );
            //printf( "BIN2: dev0=%g, dev=%g sig=%g\n", dev0, dev, sig[j] );
            if( !do_curve ) {
                ymean = bin_mean[bin1];
                dev0 = fabs(y(j) - ymean);
            }
            /*  note sigma correction here  */
            dev = fmax( dev0 - tmp, 0. );
            sumy2_dev[bin1] += sq( dev  );
        }
        
        //printf( "BIN1:  dev0=%g, dev=%g sig=%g\n", dev0, dev, sig[j] );
    }
    
    for( j = 0; j <= nbins-1; j++ ) {
        if( nbin[j] > 1 )  bin_var[j] = sumy2_dev[j] / (nbin[j] - 1);
        else               bin_var[j] = 0.;
    }
    return(1);
}


/*------------------ compute phase at time t --------------------------------*/

double dophase( double tt, double t0, double tn, double f )
{
    double t1;
    
    if( tt < t0 )  printf( "\ndophase:  t = %g < t0 = %g\n", tt, t0 );
    
    /*  this line aded 7/23/2012 - rfs  */
    t0 += phase_shift/f;
    
    t1 = (tt - t0) * f * (1. - beta * (tt - tn));   /*  beta is the period change rate  */
    return( t1 - floor( t1 ) );
}



/*----------------- set segments ---------------------------------------------*/

int segset( int ne, double datx[], double daty[], double segdev ) 
{
    int i, ndt=0;
    double dt, dtmin, dtavg, dtrange, dtsum=0;
    
    seg_beg_i[1] = 1;
    
    dtrange = x(ne) - x(1);
    dtavg = dtrange / ne;
    dtmin = x(2) - x(1);
    nseg = 1;
    for( i = 1; i <= ne-1; i++ ) {
        if( y(i) == MISSING )  continue;
        dt = x(i+1) - x(i);
        if( !dt )  continue;       //  drop dup pts
        if( dt > segdev * dtavg ) {
            nseg++;
            if( nseg > MAXSEGS ) {
                nseg--;
                printf( "***Warning - max segs (%d) reached, last %d points no seg'ed\n", MAXSEGS, ne-i );
                break;
            }
            seg_beg_i[nseg] = i+1;
        }
        else {
            ndt++;
            dtsum += dt;
            dtmin = fmin( dtmin, dt);
        }
    }
    if( ndt ) {
        dt_avg = dtsum / ndt;
        dt_min = dtmin;
    }
    // if( pdm_verbose ) {
        // fprintf( stdout, "\nAuto-Segmentation: segdev = %g,  nseg = %d,  dtrange = %g, dtavg = %g\n", segdev, nseg, dtrange, dtavg );
    // }
    return(1);
}


/*  heap sort:
sort on dat1, but also interchange dat2 use sgn to get reverse sort
if two dat1 elements are equal, use dat2 to order
dat3 also reordered, but not tested
n = number of elements to sort, so partial arrays can be handled
array index runs from 1 -> n                                       */

int p_sort( int n, double dat1[], double dat2[], double dat3[], int sgn ) 
{
    int ihire, j, iret, i;
    double tmp1, tmp2, tmp3;
    
    ihire = (n >> 1) + 1;
    iret = n;
    for( ;; ) {
        if( ihire > 1 ) {
            tmp1 = dat1[--ihire];       /*  first hire workers  */
            tmp2 = dat2[ihire];
            tmp3 = dat3[ihire];
        }
        else {
            tmp1 = dat1[iret];
            tmp2 = dat2[iret];
            tmp3 = dat3[iret];
            dat1[iret] = dat1[1];
            dat2[iret] = dat2[1];
            dat3[iret] = dat3[1];
            if( --iret == 1 ) {         /*  then retire them  */
                dat1[1] = tmp1;
                dat2[1] = tmp2;
                dat3[1] = tmp3;
                return(1);
            }
        }
        /*  bubble tmp1 down to its level  */
        i = ihire;                      /*  current level  */
        j = ihire << 1;                 /*  first underling   */
        while( j <= iret ) {
            if( (j < iret && sgn * dat1[j] < sgn * dat1[j+1]) ||
                (dat1[j] == dat1[j+1] && sgn * dat2[j] < sgn * dat2[j+1]) )  ++j;
            if( (sgn * tmp1 < sgn * dat1[j]) ||
                ((tmp1 == dat1[j]) && (sgn * tmp2 < sgn * dat2[j])) ) {
                dat1[i] = dat1[j];
                dat2[i] = dat2[j];
                dat3[i] = dat3[j];
                j += (i=j);
            }
            else j = iret + 1;
        }
        dat1[i] = tmp1;
        dat2[i] = tmp2;
        dat3[i] = tmp3;
    }
}

/*  same as above, but modified for data reordering  */

int i_sort( int n, int dat1[], double dat2[], int sgn ) 
{
    int ihire, j, iret, i;
    int tmp1; 
    double tmp2;
    
    ihire = (n >> 1) + 1;
    iret = n;
    for( ;; ) {
        if( ihire > 1 ) {
            tmp1 = dat1[--ihire];       /*  first hire workers  */
            tmp2 = dat2[ihire];
        }
        else {
            tmp1 = dat1[iret];
            tmp2 = dat2[iret];
            dat1[iret] = dat1[1];
            dat2[iret] = dat2[1];
            if( --iret == 1 ) {         /*  then retire them  */
                dat1[1] = tmp1;
                dat2[1] = tmp2;
                return(1);
            }
        }
        /*  bubble tmp1 down to its level  */
        i = ihire;                      /*  current level  */
        j = ihire << 1;                 /*  first underling   */
        while( j <= iret ) {
            if( (j < iret && sgn * dat1[j] < sgn * dat1[j+1]) ||
                (dat1[j] == dat1[j+1] && sgn * dat2[j] < sgn * dat2[j+1]) )  ++j;
            if( (sgn * tmp1 < sgn * dat1[j]) ||
                ((tmp1 == dat1[j]) && (sgn * tmp2 < sgn * dat2[j])) ) {
                dat1[i] = dat1[j];
                dat2[i] = dat2[j];
                j += (i=j);
            }
            else j = iret + 1;
        }
        dat1[i] = tmp1;
        dat2[i] = tmp2;
    }
}


/*  linear interpolate x0 in the table xt / yt  */
/*  xt assumed monotonic, n is the array size   */

double table_interp( double x0, int n, double xt[], double yt[] )
{
    double tmp;
    int i;
    
    if( x0 < xt[1] || x0 > xt[n] )  return -1; //error( "table_interp: value out of range" );
    for( i = 1; i <= n-1; i++ ) {
        if( x0 > xt[i] && x0 < xt[i+1] ) break;
    }
    tmp = yt[i] + (x0 - xt[i]) * (yt[i+1] - yt[i]) / (xt[i+1] - xt[i]);
    
    //printf( "interp x0 = %g, n=%d, i = %d, ret=%g\n", x0, n, i, tmp );
    //printf( "      xs = %g %g  ys = %g %g\n", xt[i], xt[i+1], yt[i], yt[i+1] );
    return( tmp );
}


double sqrc( double x )
{
    if( x <= 0. )  return( 0. );
    else           return( sqrt( x ) );
}



double inc_beta( double a, double b, double x )
{
    double beta2( double a, double b, double x );
    double gaml( double xx );
    double bt;
    
    if( x <=  0.0 || x >=  1.0 ) {
        bt = 0.0;
    }
    else {
        bt = exp(gaml(a + b) - gaml( a ) - gaml( b ) + a * log( x ) + b * log( 1.0 - x ));
    }
    if( x < (a + 1.0) / (a + b + 2.0) ) {
        return( bt * beta2( a, b, x) / a );
    }
    else {
        return( 1.0 - bt * beta2( b, a, 1.0 - x) / b );
    }
}


double beta2( double a, double b, double x )
{
    int m, m2;
    double aa, c, d, del, h, qab, qam, qap;
    
    qab = a + b;
    qap = a + 1.0;
    qam = a - 1.0;
    c = 1.0;
    d = 1.0 - qab * x / qap;
    if( fabs( d ) < FPMIN )  d = FPMIN;
    d = 1.0 / d;
    h = d;
    for( m = 1; m <= MAXIT; m++ ) {
        m2 = 2 * m;
        aa = m * (b - m) * x / ((qam + m2) * (a + m2));
        d = 1.0 + aa * d;
        if( fabs( d) < FPMIN ) d = FPMIN;
        c = 1.0 + aa / c;
        if( fabs( c ) < FPMIN ) c = FPMIN;
        d = 1.0 / d;
        h  *=  d * c;
        aa = -(a + m) * (qab + m) * x / ((a + m2) * (qap + m2));
        d = 1.0 + aa * d;
        if( fabs( d ) < FPMIN ) d = FPMIN;
        c = 1.0 + aa / c;
        if( fabs( c ) < FPMIN ) c = FPMIN;
        d = 1.0 / d;
        del = d * c;
        h  *=  del;
        if( fabs( del - 1.0 ) < EPS) break;
    }
    if( m > MAXIT ) {
        printf( "a or b too big, or MAXIT too small in beta2\n" );
    }
    return( h );
}

double gaml( double xx )
{
    double x, y, tmp, ser;
    static double cof[6] = {76.18009172947146, -86.50532032941677,
        24.01409824083091, -1.231739572450155,
        0.1208650973866179e-2, -0.5395239384953e-5};
    int j;
    
    y = x = xx;
    tmp = x + 5.5;
    tmp  -=  (x + 0.5) * log( tmp );
    ser = 1.000000000190015;
    for( j = 0; j <= 5; j++ ) ser += cof[j] / ++y;
    return( -tmp + log( (2.5066282746310005 * ser) / x ) );
}


double rnd()
{
    /* return( (double)rand() / (RAND_MAX+1.) ); */
    return( ran1( &dummy ) );
}



/*  NRC favorite generator, returns rans 0->1, period = 1.e8  */
/*  called here using rnd()  */

#define IA 16807
#define IM 2147483647
#define AM (1.0/IM)
#define IQ 127773
#define IR 2836
#define NTAB 32
#define NDIV (1+(IM-1)/NTAB)
#define RNMX (1.0-EPS)

double ran1(long *idum)
{
    int j;
    long k;
    static long iy=0;
    static long iv[NTAB];
    double temp;
    
    if (*idum <= 0 || !iy) {
        if (-(*idum) < 1) *idum=1;
        else *idum = -(*idum);
        for (j=NTAB+7;j>=0;j--) {
            k=(*idum)/IQ;
            *idum=IA*(*idum-k*IQ)-IR*k;
            if (*idum < 0) *idum += IM;
            if (j < NTAB) iv[j] = *idum;
        }
        iy=iv[0];
    }
    k=(*idum)/IQ;
    *idum=IA*(*idum-k*IQ)-IR*k;
    if (*idum < 0) *idum += IM;
    j=iy/NDIV;
    iy=iv[j];
    iv[j] = *idum;
    if ((temp=AM*iy) > RNMX) return RNMX;
    else return temp;
}


/*  old UNIX function to read a line of data  */

// int fgetline(FILE *fp, char *buffer, int max_len)
// {
//     int len, ch;
    
//     ch = getc(fp);
//     if (ch == EOF)
//         return(EOF);
    
//     len = 0;
//     while (ch != '\n' && ch != EOF ) {
//         buffer[len] = ch;
//         len ++;
//         if (len >= max_len)
//             return(len);
//         ch = getc(fp);
//     }
//     buffer[len] = '\0';
//     return(len);
// }
