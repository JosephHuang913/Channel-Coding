// ------------------------------------------------------------------------
// File: shannon2.c
//
// Compute Shannon's capacity in terms of Eb/N0 (dB)
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
double ebn0;
FILE   *fp;
char   name2[40];
double rate, rate_step;

  if (argc != 3)
    {
      printf("Usage: %s rate_step output_file \n", argv[0]);
      exit(0);
    }

  sscanf(argv[1],"%lf", &rate_step);
  sscanf(argv[2],"%s", name2);
  fp = fopen(name2,"w");

  rate = rate_step;
  while ( rate < 5 ) {
    ebn0 = (pow(2.0,(rate))-1.0)/(rate);
    ebn0 = 10.0 * log10(ebn0);
    printf("%lf \t %lf\n", ebn0, rate);
    fprintf(fp, "%lf \t %lf\n", ebn0, rate);
    fflush(stdout);
    rate += rate_step;
    }
}

