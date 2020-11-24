// ------------------------------------------------------------------------
// File:    bitf.c 
// Created: February 14, 2000
//
// Gallager's iterative BIT-FLIP decoding of linear block codes
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
#include <stdlib.h>

#define MAX_RANDOM     LONG_MAX   // Maximum value of random() 
#define NODES             16384   // Maximum of number of code/check nodes
#define J                    17   // Maximum number of checks per code bit
#define K                    17   // Maximum number of code bits per check

int max_size_M;
int max_size_N;

int threshold;

int n;                            // length
int k;                            // dimension
int nk;                           // redundancy
float  rate;                      // code rate

int M,N;                          // Size of parity-check matrix

// ---------------
// NODE STRUCTURES
// ---------------

struct parent_node {
    int size;
    int index[J];                     // indexes of children
    };

struct child_node {
    int size;
    int index[K];                     // indexes of parents 
    int syndrome;
    };

struct parent_node code_node[NODES];
struct child_node check_node[NODES];

double init_snr, final_snr, snr_increment;
double sim, num_sim, ber, amp;
long seed;
int error;
int data[NODES], codeword[NODES];
int data_int;
double snr;
float  transmited[NODES], received[NODES];
int hard[NODES], decoded[NODES];
int i,j, iter, max_iter;
char filename[40], name2[40];
FILE *fp, *fp2;

void initialize(void);
void awgn(void);
void encode(void);
void bitflip(void);

main(int argc, char *argv[])
{
  // Command line processing
  if (argc != 12)
    {
      printf("Usage: %s length(n) dimension(k) file_parity-check threshold max_iter init_snr final_snr snr_inc num_sim output_file seed\n", 
                      argv[0]);
      exit(0);
    }

  sscanf(argv[1],"%d", &n);
  sscanf(argv[2],"%d", &k);
  sscanf(argv[3],"%s", filename);
  sscanf(argv[4],"%d", &threshold);
  sscanf(argv[5],"%d", &max_iter);
  sscanf(argv[6],"%lf", &init_snr);
  sscanf(argv[7],"%lf", &final_snr);
  sscanf(argv[8],"%lf", &snr_increment);
  sscanf(argv[9],"%lf",&num_sim);
  sscanf(argv[10],"%s", name2);
  sscanf(argv[11],"%ld",&seed);

  nk = n-k;
  rate = (float) k / (float) n;

  printf("\nGALLAGER'S BIT-FLIP DECODING OF LINEAR BLOCK CODES\n");
  printf("threshold = %d \tmax_iter = %d\n", threshold, max_iter);
  printf("SNR from %lf to %lf in increments of %lf\n",
         init_snr, final_snr, snr_increment);
  printf("%.0f codewords transmitted per SNR\n", num_sim);
  printf("\nn=%d, k=%d, n-k=%d, and rate = %lf\n\n",n,k,nk,rate);

  if ((fp = fopen(filename,"r")) != NULL)
    {
      fscanf(fp, "%d %d", &N, &M);
      fscanf(fp, "%d %d", &max_size_N, &max_size_M);
      for (i=0; i<M; i++)
        fscanf(fp, "%d", &check_node[i].size);
      for (i=0; i<N; i++)
        fscanf(fp, "%d", &code_node[i].size);

      // Read index sets for check nodes
      for (i=0; i<M; i++)
      {
        for (j=0; j<check_node[i].size; j++)
          fscanf(fp, "%d", &check_node[i].index[j]);
      }

      // Read index sets for code nodes
      for (i=0; i<N; i++)
      {
        for (j=0; j<code_node[i].size; j++)
          fscanf(fp, "%d", &code_node[i].index[j]);
      }
    }
  else 
    { 
      printf("incorrect input file name ...\n"); 
      exit(0);
    }

  fclose(fp);

#ifdef PRINT_MATRIX
  printf("%d %d\n", N, M);
  printf("%d %d\n", max_size_N, max_size_M);
  for (i=0; i<M; i++)
    printf("%4d", check_node[i].size);
  printf("\n");
  for (i=0; i<N; i++)
    printf("%4d", code_node[i].size);
  printf("\n");
  for (i=0; i<M; i++)
  {
    for (j=0; j<check_node[i].size; j++)
      printf("%4d", check_node[i].index[j]);
    printf("\n");
  }
  for (i=0; i<N; i++)
  {
    for (j=0; j<code_node[i].size; j++)
      printf("%4d", code_node[i].index[j]);
    printf("\n");
  }
  printf("\n");
#endif

  fp2 = fopen(name2,"w");

  snr = init_snr;
  srandom(seed);

  // -------------------------------------------------------------------
  //                  S I M U L A T I O N   L O O P 
  // -------------------------------------------------------------------

  while ( snr < (final_snr+0.001) )
    {

      initialize();

      while (sim < num_sim)   //  <--- Fixed number of simulations
      // while (ber < 1000)      //  <-- Minimum number of errors
      { 

        // ----------        FOR CONVENIENCE, MAKE DATA EQUAL TO ZERO
        for (i=0; i<k; i++)
          data[i] = 0;

        // -----------       BPSK MAPPING: "0" --> +1,  "1" --> -1
        for (i=0; i<n; i++)
          transmited[i] = 1.0;

        // -----------       ADDITIVE WHITE GAUSSIAN NOISE CHANNEL
        awgn();

        // -----------       BPSK SYMBOL-BY-SYMBOL ESTIMATION
        for (i=0; i<n; i++)
          if (received[i] > 0.0)
            hard[i] = 0;
          else
            hard[i] = 1;

        // -----------       ITERATIVE BIT-FLIP DECODING
        bitflip();

        // -----------       COUNT THE NUMBER OF BIT ERRORS
        for (i=0; i<n; i++)
          if (decoded[i]) ber++;

        sim+=1.0;

      }

    printf("%lf \t%8.0lf %8.0lf \t%13.8e\n", snr, ber, (n*sim), (ber/(sim*n))); 
    fflush(stdout);
    fprintf(fp2, "%lf %13.8e\n", snr, (ber/(sim*n)) );
    fflush(fp2);

    snr += snr_increment;

  }

  fclose(fp2);

}

void bitflip()
{
// BIT-FLIP decoding

int i,j,l,iter;
int delt,m,aux;
int all_zero;                             // Flag for syndrome testing
int count;

  // -------------------
  // INITIALIZATION STEP
  // -------------------

  // Prior values (used to be probabilities in soft-decision)

  for (i=0;i<N;i++)
    {
    decoded[i] = hard[i];
    }

  iter = 0;                  // Counter of iterations

  do {

  // ---------------------------------------
  // HORIZONTAL STEP = BOTTOM-UP PROPAGATION
  // ---------------------------------------
  //
  // Run through the checks m and compute, for each n in N(m) the
  // probabilitiy of a check symbol when code symbol is 0 (or 1)
  // given that the other code symbols have values 0, 1
  //
  // Pearl:
  // Node x_m computes new "lambda" messages to be sent to its parents
  // u_1, u_2, ..., u_K

  // Flag to determine if syndrome is all zero
  all_zero = 1;

  for (i=0; i<M; i++)
  {
    delt = 0;
    for (j=0; j<check_node[i].size; j++)
    {
      aux = check_node[i].index[j];
      delt ^= decoded[aux-1];
    }
  check_node[i].syndrome = delt;

  // Check if anyone of the syndromes is not zero
  if (delt) all_zero = 0;
  }


  // CONTINUE IF A CODEWORD HAS NOT BEEN FOUND
  if (!all_zero) {

  // ------------------------------------
  // VERTICAL STEP = TOP-DOWN PROPAGATION
  // ------------------------------------
  //
  // MacKay:
  // Take the computed values of rm0, rm1 and update the values of
  // the probabilities qm0, qm1
  //
  // Pearl:
  // Each node u_l computes new "pi" messages to be send to its
  // children x_1, x_2, ..., x_J

  for (i=0; i<N; i++)
    {
    count = 0;
    for (j=0; j<code_node[i].size; j++)
      {
      aux = code_node[i].index[j]-1; 

      // Compute index "m" of message from children
      m = 0;
      while (  ( (check_node[aux].index[m]-1) != i )
                     && ( m < check_node[aux].size )  ) m++;

      if (check_node[aux].syndrome)
        count++;
      }

    // If more that 1/2 checks are unsatisfied, FLIP the BIT
    if (count > threshold)
    // if (count > (code_node[i].size-1)/2 )
      {
      decoded[i] ^= 1;
      }
    }
  }

  // Increment the number of iterations, and check if maximum reached
  iter++;

  } while (iter < max_iter);

}

void encode()
//
// Systematic encoding 
//
{
//int i,j;
//for (j=0; j<n; j++)
//  {
//    if (j<k)                     // information position
//      codeword[j] = data[j];
//    else                         // redundant position
//      {
//        codeword[j] = 0;
//        for (i=0; i<k; i++)
//          // codeword[j] ^= ( data[i] * H[j-k][i] ) & 0x01;
//          codeword[j] ^= ( data[i] * H[(j-k+i)%n] ) & 0x01;
//      }
//  }
}


void awgn()
//
// AWGN generation
// 
{
  double u1,u2,s,noise,randmum;
  int i;
  for (i=0; i<n; i++)
    {
      do {
            randmum = (double)(random())/MAX_RANDOM;
            u1 = randmum*2.0 - 1.0;
            randmum = (double)(random())/MAX_RANDOM;
            u2 = randmum*2.0 - 1.0;
            s = u1*u1 + u2*u2;
            } while( s >= 1);
      noise = u1 * sqrt( (-2.0*log(s))/s );
      received[i] = transmited[i] + noise/amp;
#ifdef NO_NOISE
      received[i] = transmited[i];
#endif
    }
}


void initialize()
{
  amp = sqrt(2.0*rate*pow(10.0,(snr/10.0)));
  ber = 0.0;
  sim = 0.0;
}


