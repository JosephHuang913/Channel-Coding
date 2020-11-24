// ------------------------------------------------------------------------
// File: shannon.c
//
// Compute Shannon's capacity
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

main(int argc, char *argv[])
{
int    i, j, k; 
long   seed;
int    mod;
double snr, init_snr, final_snr, snr_increment;
double ebn0;
double sim, num_sim;
FILE   *fp;
char   name2[40];
int    npoints;
double x[16], y[16];
double sum, aux, aux2, N0;
double z, rate;

  // Command line processing
  if (argc != 5)
    {
      printf("Usage: %s init_snr final_snr snr_inc output_file \n", 
                      argv[0]);
      exit(0);
    }

  sscanf(argv[1],"%lf", &init_snr);
  sscanf(argv[2],"%lf", &final_snr);
  sscanf(argv[3],"%lf", &snr_increment);
  sscanf(argv[4],"%s", name2);

  fp = fopen(name2,"w");

  snr = init_snr;
  while ( snr < (final_snr+0.001) ) {
    aux = pow(10.0,(snr/10.0));
    rate = 0.5 * log(1+aux)/M_LN2;
    printf("%lf \t %lf\n", snr, rate);
    fprintf(fp, "%lf \t %lf\n", snr, rate);
    fflush(stdout);
    snr += snr_increment;
    }
}

