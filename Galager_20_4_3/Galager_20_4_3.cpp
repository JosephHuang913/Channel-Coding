/*=================================================================*/
/* Author: Chao-wang Huang                                                                                        */
/* Date: Tuesday, June 16, 2009                                                                                 */
/* A (20,7) Galager Low-Density Parity-Check code (LDPC) is simulated                     */
/* Decoding algorithm: Message Passing Algorithm (Log-Domain Belief Propagation)  */
/*=================================================================*/

#include "stdafx.h"
#include <math.h>
#include <cmath>
#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include <ctime>
#include <conio.h>
#include <limits.h>
#include <float.h>
#include <cstdlib>
#include <iostream>
using namespace std;
#include <fstream>
#include <cstdio>
#include <assert.h>

void AWGN_noise(double, double, double *);
int Error_count(int, int);
void LDPC_Encode(int *, int *);
int Belief_Propagation(double *, int *, double, int);
double F(double);

#define NODES          20							// Maximum of number of code/check nodes
#define MAX_CHK      3   						// Maximum number of checks per code bit
#define MAX_BIT        4   						// Maximum number of code bits per check
#define Num_packet 	500000				// number of packets simulated
int Max_col_weight;                 					// Maximum column weight
int Max_row_weight;                 				// Maximum row weight
const int Iteration = 50;	   						// Iterations of message passing ( >= 1 )
const int n = 20;               						// Codeword length of a linear block code
const int k = 7;                							// Information bits of a linear block code
int n_k;                          							// Redundancy (Parity check bits)
int M, N;                          							// Size of parity-check matrix (Row, Column) = (M, N)
int **G, **H;												// Generator matrix & parity check matrix
int *coded_bit, *err_count;

// ---------------
// NODE STRUCTURES
// ---------------
struct bit_node {
	int Num_chks_bit;									// Number of check nodes of a bit node
    int index[MAX_CHK];							// Check nodes set of a bit node
    double pi1[MAX_CHK];							// messages "pi" to check node
    };

struct check_node {
    int Num_bits_chk;									// Number of bit nodes of a check node		
    int index[MAX_BIT];								// Bit nodes set of a check node
    double lambda1[MAX_BIT];					// messages "lambda" to bit node
    };

struct bit_node Bit_node[NODES];
struct check_node CHK_node[NODES];

int main(void)
{
	time_t t, start, end;
	int i, j, p, iter, *data_bit, *decoded;
	double *bpsk, snr, Eb_No, noise_pwr, noise[2], *Yk, *err_rate;
	FILE *parity_check, *ber, *records;

	start = time(NULL);
	printf("BER Performance of (20,7) LDPC Code in AWGN Channel\n");
	printf("Dimension of Parity Check matrix: (15,20)\n");
	printf("Column Weight = 3, Row Weight = 4\n");
	printf("Number of bits of simulation = %d\n\n", n*Num_packet);
	printf("This program is running. Don't close, please!\n\n");

	fopen_s(&records, "Records_20_7_BP_awgn.log", "a");
	fprintf(records, "BER Performance of Galager (20,7) LDPC Code in AWGN Channel\n");
	fprintf(records, "Dimension of Parity Check matrix: (15,20)\n");
	fprintf(records, "Column Weight = 3, Row Weight = 4\n");
	fprintf(records, "number of bits of simulation = %d\n\n", n*Num_packet);
	fprintf(records, "Eb/No     BER\n");
	fflush(records);

	data_bit = new int[k];
	coded_bit = new int[n];
	decoded = new int[n];
	bpsk = new double[n];
	Yk = new double[n];								// Received signal
	err_count = new int[Iteration];
	err_rate = new double[Iteration];
	H = new int*[15];
	for(i=0; i<15; i++)
		H[i] = new int[20];
	G = new int*[k];
	for(i=0; i<k; i++)
		G[i] = new int[n];

	srand((unsigned) time(&t));

	// Define Parity Check Matrix
	fopen_s(&parity_check, "gallager_20_4_3.log", "r");
	fscanf_s(parity_check, "%d %d", &N, &M);
	fscanf_s(parity_check, "%d %d", &Max_row_weight, &Max_col_weight);
	
	for (i=0; i<M; i++)
        fscanf_s(parity_check, "%d", &CHK_node[i].Num_bits_chk);
    for (i=0; i<N; i++)
        fscanf_s(parity_check, "%d", &Bit_node[i].Num_chks_bit);

	// Read index sets for check nodes
    for (i=0; i<M; i++)
		for (j=0; j<CHK_node[i].Num_bits_chk; j++)
			fscanf_s(parity_check, "%d", &CHK_node[i].index[j]);
    
    // Read index sets for bit nodes
    for (i=0; i<N; i++)
    {
		for (j=0; j<Bit_node[i].Num_chks_bit; j++)
			fscanf_s(parity_check, "%d", &Bit_node[i].index[j]);
    }

   fclose(parity_check);

	for(i=0; i<M; i++)
		for(j=0; j<N; j++)
			H[i][j] = 0;

   for(i=0; i<M; i++)
	   for(j=0; j<CHK_node[i].Num_bits_chk; j++)
		   H[i][CHK_node[i].index[j]-1] = 1;
/*
	fopen_s(&parity_check, "Parity_20_4_3.log", "w");
	for(i=0; i<M; i++)
	{
		for(j=0; j<N; j++)
			fprintf(parity_check, "%d ", H[i][j]);
		fprintf(parity_check, "\n");
	}
	fclose(parity_check);
*/
	fopen_s(&parity_check, "G.log", "r");
	for(i=0; i<k; i++)
		for(j=0; j<n; j++)
			fscanf_s(parity_check, "%d", &G[i][j]);
	fclose(parity_check);

/*===========================================*/
/*             M A I N    S I M U L A T I O N    L O O P           */
/*===========================================*/
	fopen_s(&ber, "ber_Galager_20_4_3.log", "w");
	for(snr=0; snr<=10; snr+=1)
	{
   		for(iter=0; iter<Iteration; iter++)
   			err_count[iter] = 0;
   		
		// noise power calculation
		Eb_No = (double)snr;
		noise_pwr = 1.0/(((double)k/(double)n)*pow(10.0, Eb_No/10.0));	// BPSK, Nyquist filter assumption
		printf("Eb_No = %f, ", Eb_No);

		for(p=0; p<Num_packet; p++)
		{
			// Generate random information bit stream
      		for (i=0; i<k; i++)
				if(rand()/(double)RAND_MAX>=0.5)
					data_bit[i] = 1;
				else
					data_bit[i] = 0;

			/* LDPC code Encoding */
			LDPC_Encode(data_bit, coded_bit);

			for(i=0; i<n; i++)		// BPSK mapping
   				if(coded_bit[i] == 1)
      				bpsk[i] = 1.0;
      			else
      		  		bpsk[i] = -1.0;

			/* AWGN channel */
			for(i=0; i<n; i++)
			{
         		AWGN_noise(0, noise_pwr, &noise[0]);
         		Yk[i] = bpsk[i] + noise[0];
				
				// Soft Decision
				//Dk[i] = tanh(Yk[i]/noise_pwr);
				//Dk[i] = Yk[i];
	        }

/* Message Passing Decoding (Belief Propagation) */
			
			Belief_Propagation(Yk,  decoded, snr, Iteration);
		}

		// Statistics and records
		cout << "Error Rate = ";
		for(iter=0; iter<Iteration; iter++)
		{
      		err_rate[iter] = err_count[iter] / (double)(N*Num_packet);
	      	printf("%e, ", err_rate[iter]);
		}
		cout << endl;

		fprintf(ber, "%f ", Eb_No);
		fprintf(records, "%f ", Eb_No);
		for(iter=0; iter<Iteration; iter++)
		{
			fprintf(ber, "%e ", err_rate[iter]);
			fprintf(records, "%e ", err_rate[iter]);
		}
		fprintf(ber, "\n");
		fprintf(records, "\n");
		fflush(records);
		fflush(ber);
	}

	delete data_bit;
	delete coded_bit;
	delete decoded;
	delete bpsk;
	delete Yk;
	delete err_count;
	delete err_rate;
	for(i=0; i<M; i++)
		delete H[i];
	delete H;
	for(i=0; i<k; i++)
		delete G[i];
	delete G;
	

	end = time(NULL);
	printf("Total elapsed time: %.0f(sec)\n", difftime(end,start));
	fprintf(records, "\nTotal elapsed time: %.0f(sec)\n\n", difftime(end,start));
	fclose(ber);
	fclose(records);
	printf("This program is ended. Press any key to continue.\n");
	getchar();

	return 0;
}

int Belief_Propagation(double *received, int *decoded, double snr, int max_iter)
{
// Iterative decoding by belief propagation in code's Bayesian network
// Based on Pearl's book and MacKay's paper

	int i, j, l, iter, m, aux, sign;
	double delt;
	double llrp1[NODES];									// Prior probabilities (channel)
	double q1[NODES];										// Pseudo-posterior probabilities
	double snr_rms = 2.0 * sqrt(2.0*(double)k/(double)n*(pow(10.0,(snr/10.0))));

	// -------------------
	// ***** STEP 0 *****
	// INITIALIZATION STEP
	// -------------------

	// Prior log-likelihood ratios (channel metrics)
	for (i=0;i<N;i++)
    {
		// LOOK-UP TABLE (LUT)
		llrp1[i] = received[i]*snr_rms;	// ln{P(1)/P(-1)}
    }

	// For every (m,l) such that there is a link between parents and
	// children, qm0[i][j] and qm1[i][j] are initialized to pl[j].
	// Notation: pi (Pearl) = q (MacKay)

	for (i=0; i<N; i++)                         // run over code nodes
    {
		//printf("received = %10.7lf\n", received[i]);
		for (j=0; j<Bit_node[i].Num_chks_bit; j++)       // run over check nodes
		{
			Bit_node[i].pi1[j] = llrp1[i];
			//printf("i, j,  pi1 = %d,%d   %10.7lf\n", i, j, Bit_node[i].pi1[j]);
		}
	}

	iter = 0;                  // Counter of iterations
	do {

	// ---------------------------------------
	//         ***** STEP 1 *****
	// HORIZONTAL STEP = BOTTOM-UP PROPAGATION
	// ---------------------------------------
	//
	// MacKay:
	// Run through the checks m and compute, for each n in N(m) the
	// probabilitiy of a check symbol when code symbol is 0 (or 1)
	// given that the other code symbols have distribution qm0, qm1
	//
	// Pearl:
	// Node x_m computes new "lambda" messages to be sent to its parents
	// u_1, u_2, ..., u_K

		for (i=0; i<M; i++)
			for (j=0; j<CHK_node[i].Num_bits_chk; j++)
			{	
				delt = 0.0;
				sign = 0;                           // Keep track of sign of delt, 0: positive, 1: negative

				for (l=0; l<CHK_node[i].Num_bits_chk; l++)
				{
					aux = CHK_node[i].index[l];
					if (aux != CHK_node[i].index[j])
					{
						// --------------------------------------------------------
						//  Compute the index "m" of the message from parent node
						// --------------------------------------------------------
						m = 0;
						while (  ( (Bit_node[aux-1].index[m]-1) != i ) && ( m < Bit_node[aux-1].Num_chks_bit)  ) 
							m++;

						if (Bit_node[aux-1].pi1[m] < 0.0) 
							sign ^= 1;
					
						delt += F(fabs(Bit_node[aux-1].pi1[m]));
						//printf("pi1, delt =  %lf, %lf \n", Bit_node[aux-1].pi1[m], delt);
					}
				}	
      
				if (sign == 0)
					CHK_node[i].lambda1[j] = F(delt);
				else
					CHK_node[i].lambda1[j] = -F(delt);

				// Normalization
				if (CHK_node[i].lambda1[j] < -30.0)
					CHK_node[i].lambda1[j] = -30.0;

				//printf("i, j, lambda1 = %d,%d  %10.7lf\n", i, j, CHK_node[i].lambda1[j]);
			}

		// ------------------------------------
		//         ***** STEP 2 *****
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
			for (j=0; j<Bit_node[i].Num_chks_bit; j++)
			{
				Bit_node[i].pi1[j] = 0.0;

				for (l=0; l<Bit_node[i].Num_chks_bit; l++)
				{
					aux = Bit_node[i].index[l] - 1; 

					if ( aux != (Bit_node[i].index[j]-1) )
					{
						// Compute index "m" of message from children
						m = 0;
						while (  ( (CHK_node[aux].index[m]-1) != i ) && ( m < CHK_node[aux].Num_bits_chk )  ) 
							m++;

						Bit_node[i].pi1[j] += CHK_node[aux].lambda1[m];
					}
				}

				Bit_node[i].pi1[j] += llrp1[i];
				//printf("---->  i, j,  pi1 = %d,%d    %10.7lf\n", i, j, Bit_node[i].pi1[j]);

				if (Bit_node[i].pi1[j] < -30.0)
					Bit_node[i].pi1[j] = -30.0;
			}

		// DECODING:
		// MacKay: At this step we also compute the (unconditional) pseudo-
		// posterior probalilities "q0, q1" to make tentative decisions

		for (i=0; i<N; i++)
	    {
			q1[i] = 0.0;

			for (j=0; j<Bit_node[i].Num_chks_bit; j++)
			{
				aux = Bit_node[i].index[j] - 1; 

				// Compute index "m" of message from children
				m = 0;
				while (  ( (CHK_node[aux].index[m]-1) != i ) && ( m < CHK_node[aux].Num_bits_chk )  ) 
					m++;

				q1[i] += CHK_node[aux].lambda1[m];
			}

			q1[i] += llrp1[i];
			//printf("iter = %d, q1 = %10.7lf\n", iter, q1[i]);

			if (q1[i] >= 0.0) 
				decoded[i] = 1;
			else 
				decoded[i] = 0;

			err_count[iter] += Error_count(coded_bit[i], decoded[i]);
		}

/*
		sum = 0;
		for(i=0; i<M; i++)
			for(j=0; j<N; j++)
				sum += decoded[j] * H[i][j]
		if(sum%2 = 0)
			return 0;
*/
		// Increment the number of iterations, and check if maximum reached
		iter++;
	} while (iter < max_iter);

	return 0;
}

double F(double x)
{
	double interm;
	
	if (x == 0.0) 
		return(1.0e+30);
	if (fabs(x) > 30.0)
		interm = 0.0;
	else
		interm = log ( (exp(x)+1.0)/(exp(x)-1.0) );
	
	return(interm);
}

void LDPC_Encode(int *data, int *codeword)
// Systematic encoding 
{
	int i, j;
	
	for (j=0; j<n; j++)
	{
		if (j>=n-k)                     // information bits
			codeword[j] = data[j-(n-k)];
		else								// parity bits
		{
			codeword[j] = 0;
			for (i=0; i<k; i++)
				codeword[j] ^= ( data[i] * G[i][j] ) & 0x01;
		}
	}
}

void AWGN_noise(double mu, double variance, double *noise)
{
	double Pi = 3.14159265358979;
   double u1, u2;
   do
   {
   	u1 = (double)rand()/(double)RAND_MAX;
      u2 = (double)rand()/(double)RAND_MAX;
   }
   while(u1 == 0.0 || u2 == 0.0);

   *(noise+0) = (sqrt(-2.0*log(u1))*cos(2*Pi*u2))*sqrt(variance/2.0)+mu/sqrt(2.0);
   *(noise+1) = (sqrt(-2.0*log(u1))*sin(2*Pi*u2))*sqrt(variance/2.0)+mu/sqrt(2.0);
}

int Error_count(int x, int y)
{
	if(x == y)
   	return 0;
   else
   	return 1;
}
