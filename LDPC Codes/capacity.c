// ------------------------------------------------------------------------
// File: capacity.c
// 
// Compute the constellation-constrained capacity of an AWGN channel
//
// References:
//
// G. Ungerboeck, "Channel Coding with Multilevel/Phase Signalling,"
// IEEE Trans. Info. Theory, vol. IT-28, no. 1, pp. 55-67, Jan. 1982.
//
// P. E. McIllree, "Channel Capacity Calculations for M-ary N-dimensional
// Signal Sets," Master's Thesis, School of Electronic Engineering, The
// University of South Australia, February 1995.
//
// M. Abramowitz and I.A. Stegun, Eds., "Handbook of Mathematical Functions,"
// Dover Publications: New York, 1972.
// ------------------------------------------------------------------------
// This program is complementary material for the book:
//
// R.H. Morelos-Zaragoza, The Art of Error Correcting Coding, Wiley, 2002.
//
// ISBN 0471 49581 6
//
// This and other programs are available at http://the-art-of-ecc.com
//
// You may use this program for academic and personal purposes only.
// If this program is used to perform simulations whose results are
// published in a journal or book, please refer to the book above.
//
// The use of this program in a commercial product requires explicit
// written permission from the author. The author is not responsible or
// liable for damage or loss that may be caused by the use of this program.
//
// Copyright (c) 2002. Robert H. Morelos-Zaragoza. All rights reserved.
// ------------------------------------------------------------------------

#include <stdio.h>
#include <math.h>
#include <float.h>
#include <limits.h>

#define MAXRAND LONG_MAX    // for random number generation

double gaussian(void);

main(int argc, char *argv[])
{
int    i, j, k, l; 
long   seed;
int    mod;
int    dimension;
double nbits;
double snr, init_snr, final_snr, snr_increment;
double ebn0;
double sum, sum2;
FILE   *fp;
char   name2[40];
int    npoints;
double x[16], y[16];
double aux, aux2, N0;
double cs, sn;
double z, w, rate;

// Coefficients and evaluation points for Hermite integration
double wp[10] = {
7.64043285560789861e-06,
1.34364574684723431e-03,
3.38743944571450600e-02,
2.40138611094110666e-01,
6.10862633765332230e-01,
6.10862633765332230e-01,
2.40138611094110666e-01,
3.38743944571450600e-02,
1.34364574684723431e-03,
7.64043285560789861e-06 
};

// Used for one-dimensional constellations
double xp[10] = {
3.43615911883773739e+00,
2.53273167423278966e+00,
1.75668364929988163e+00,
1.03661082978951358e+00,
3.42901327223704588e-01,
-3.42901327223704588e-01,
-1.03661082978951358e+00,
-1.75668364929988163e+00,
-2.53273167423278966e+00,
-3.43615911883773739e+00
};

// Used for two-dimensional constellations
double yp[10] = {
3.43615911883773739e+00,
2.53273167423278966e+00,
1.75668364929988163e+00,
1.03661082978951358e+00,
3.42901327223704588e-01,
-3.42901327223704588e-01,
-1.03661082978951358e+00,
-1.75668364929988163e+00,
-2.53273167423278966e+00,
-3.43615911883773739e+00
};

  // Command line processing
  if (argc != 7)
    {
      printf("Usage: %s modulation (1=BPSK 2=QPSK 3=8PSK 4=16QAM 5=64QAM) init_snr final_snr snr_inc output_file seed\n", 
                      argv[0]);
      exit(0);
    }

  sscanf(argv[1],"%d", &mod);
  sscanf(argv[2],"%lf", &init_snr);
  sscanf(argv[3],"%lf", &final_snr);
  sscanf(argv[4],"%lf", &snr_increment);
  sscanf(argv[5],"%s", name2);
  sscanf(argv[6],"%ld",&seed);
  time(&seed);
  srandom(seed);
  fp = fopen(name2,"w");

  dimension = 2;
  switch (mod) {
    case 1:  npoints = 2; 
             nbits = 1.0;
             x[0] = -1.0; 
             x[1] = 1.0; 
             dimension = 1;
             break;
    case 2:  npoints = 4;
             nbits = 2.0;
             x[0] =  M_SQRT1_2; y[0] =  M_SQRT1_2;
             x[1] = -M_SQRT1_2; y[1] =  M_SQRT1_2;
             x[3] = -M_SQRT1_2; y[3] = -M_SQRT1_2;
             x[2] =  M_SQRT1_2; y[2] = -M_SQRT1_2;
             break;
    case 3:  npoints = 8;
             nbits = 3.0;
             cs = cos(M_PI/8.0); sn = sin(M_PI/8.0);
             x[0] =  cs; y[0] =  sn;
             x[1] =  sn; y[1] =  cs;
             x[3] = -sn; y[3] =  cs;
             x[2] = -cs; y[2] =  sn;
             x[4] = -cs; y[4] = -sn;
             x[5] = -sn; y[5] = -cs;
             x[7] =  sn; y[7] = -cs;
             x[6] =  cs; y[6] = -sn;
             break;
    case 4:  npoints = 16;
             nbits = 4.0;
             cs = 1.0/sqrt(10.0); // Normalized to unit energy
             x[0] =  3.0*cs; y[0] =  1.0*cs;
             x[1] =  3.0*cs; y[1] =  3.0*cs;
             x[2] =  1.0*cs; y[2] =  3.0*cs;
             x[3] =  1.0*cs; y[3] =  1.0*cs;
             x[4] = x[0]-4.0*cs; y[4] = y[0];
             x[5] = x[1]-4.0*cs; y[5] = y[1];
             x[6] = x[2]-4.0*cs; y[6] = y[2];
             x[7] = x[3]-4.0*cs; y[7] = y[3];
             x[8] = x[4]; y[8] = y[4]-4.0*cs;
             x[9] = x[5]; y[9] = y[5]-4.0*cs;
             x[10]= x[6]; y[10]= y[6]-4.0*cs;
             x[11]= x[7]; y[11]= y[7]-4.0*cs;
             x[12]= x[0]; y[12]= y[0]-4.0*cs;
             x[13]= x[1]; y[13]= y[1]-4.0*cs;
             x[14]= x[2]; y[14]= y[2]-4.0*cs;
             x[15]= x[3]; y[15]= y[3]-4.0*cs;
             break;
    default: printf("Modulation format = %d is not implemented... exiting\n",
             mod);
             exit(0); 
             break;
    }

  printf("dimension = %d\n", dimension);
  printf("    SNR      \t    EBN0      \t    RATE\n");

  snr = init_snr;
  while ( snr < (final_snr+0.001) ) {
    switch (dimension) {

      // ------------------------------------------------
      //            One-dimentional modulation
      // ------------------------------------------------
      case 1:  N0 = 2.0/pow(10.0,(snr/10.0));
               sum = 0.0;
               for (i=0; i<npoints; i++) {
                 // Numerical Hermite integration
                 for (k=0; k<10; k++) {
                   aux2 = 0;
                   for (j=0; j<npoints; j++)
                     aux2 += exp( -2.0*xp[k]*((x[i]-x[j])/sqrt(N0))
                                      -((x[i]-x[j])*(x[i]-x[j])/N0)  );
                   sum += wp[k]*log(aux2)/M_LN2;
                   }
                 }
               sum /= (npoints * sqrt(M_PI));
               rate = nbits - sum;
               ebn0 = snr - 10.0*log10(2.0*rate);
               break;

      // ------------------------------------------------
      //           Two-dimensional modulations
      // ------------------------------------------------
      case 2: N0 = 1.0/pow(10.0,(snr/10.0));
               sum2 = 0.0;
               for (i=0; i<npoints; i++) {
                 // Numerical Hermite integration
                 for (k=0; k<10; k++) {
                   sum = 0;
                   for (l=0; l<10; l++) {
                     aux2 = 0;
                     for (j=0; j<npoints; j++)
                       aux2 += exp(  -2.0*((xp[k]*(x[i]-x[j])
                                           +yp[l]*(y[i]-y[j]))/sqrt(N0))
                                             -((x[i]-x[j])*(x[i]-x[j])+
                                               (y[i]-y[j])*(y[i]-y[j]))/N0);
                     sum += wp[l]*log(aux2)/M_LN2;
                     }
                   sum2 += wp[k]*sum;
                   }
                 }
               sum2 /= (npoints*M_PI);
               rate = nbits - sum2;
               ebn0 = snr - 10.0*log10(rate);
               break;

      }

    printf("%lf \t %lf \t %lf\n", snr, ebn0, rate);
    fprintf(fp, "%lf \t %lf\n", ebn0, rate);
    fflush(stdout);
    snr += snr_increment;
  }
}

