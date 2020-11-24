/**********************************************/
/* Author: Chao-wang Huang                    */
/* Date: Saturday, August 28, 2004            */
/* An (n,k,m) Convolutional code is simulated */
/* Decoding algorithm: MAP decoding           */
/**********************************************/

#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include <conio.h>
#include <iostream.h>
#include <limits.h>
#include <float.h>

void AWGN_noise(float, double, double *);
int Error_count(int, int);
void Trellis_diagram(int, int, int, int *);
void CC_Encoder(int *, int *, int, int, int *, int, int);
void MAP_dec_BCJR(double *, double *, double *, double, int, int, int);

#define Pi 	3.14159265358979
//#define n	2					// (n,k,m) Convolutional Code, # of output symbols
//#define k	1					// (n,k,m) Convolutional Code, # of input symbols
//#define m	2					// (n,k,m) Convolutional Code, # of memories
const int n = 2;
const int k = 1;
const int m = 2;
//#define K	1000		  		// packet length of information bit stream
//#define N	K*(n/k)	 		// packet length of coded bit stream
const int K = 1000;
const int N = 2000;
#define num_packet 10000		// number of packets simulated
#define Iteration 3
const int num_state = pow(2,m);	// number of states
const int num_in_sym = pow(2,k);	// number of input symbols
int gp[2] = {5,7};   // Generator polynomial of CC given in hexadecimal and in reverse order

struct trel_diag			// Data structure of trellis diagram for each branch
		{
      	int from;		// from_state of trellis diagram
         int to;			// to_state of trellis diagram
         int in;			// input data bit of trellis diagram
         int out[n];		// output codeword symbol of trellis diagram
      };
struct trel_diag Trellis[4][2];		// Trellis[num_state][num_in_sym]

int main(void)
{
	time_t  t, start, end;
	int i, p, s, *data_bit, *coded_bit, *Ak, *err_count;
   double *intrinsic, *bpsk, snr, Eb_No, noise_pwr, noise[2], *Yk, *err_rate, *LLR;
   FILE *ber, *records;

   start = time(NULL);

   printf("BER Performance of Convolutional Code (2,1,2) in AWGN Channel\n");
	printf("Generator polynomials are {5,7} in Octal\n");
   printf("This program is running. Don't close, please!\n\n");

   records = fopen("Records_Ite_CC.log", "a");
   fprintf(records, "(2,1,2) {5,7} convolutional code, MAP decoder, AWGN channel\n");
   fprintf(records, "number of bits of simulation = %d\n\n", K*num_packet);
   fprintf(records, "Iterations = %d\n\n", Iteration);
   fprintf(records, "Eb/No     BER\n");
   fflush(records);

   data_bit = new int[K];
   coded_bit = new int[N];
   bpsk = new double[N];
   Ak = new int[K];
   Yk = new double[N];
   LLR = new double[K];
   intrinsic = new double[K];
   err_count = new int[Iteration];
   err_rate = new double[Iteration];

	srand((unsigned) time(&t));

	Trellis_diagram(n, k, m, &gp[0]);

/************************/
/* main simulation loop */
/************************/
   ber=fopen("ber_ite_CC.log", "w");
   for(snr=0; snr<=10; snr++)
   {
   	for(s=0; s<Iteration; s++)
   		err_count[s] = 0;
   	// noise power calculation
      Eb_No = (double)snr;
      noise_pwr = 1.0/(((float)k/(float)n)*pow(10.0, Eb_No/10.0));	// BPSK, Nyquist filter assumption
      printf("Eb_No = %f, ", Eb_No);
      for(p=0; p<num_packet; p++)
      {
   		for(i=0; i<K-2; i++)
				data_bit[i] = random(2);		// Generate random information bit stream
   		data_bit[K-2] = data_bit[K-1] = 0;

         CC_Encoder(data_bit, coded_bit, n, m, &gp[0], K, num_state);

			for(i=0; i<N; i++)		// BPSK mapping
   			if(coded_bit[i]==1)
      			bpsk[i] = 1.0;
      		else
      			bpsk[i] = -1.0;

         /* AWGN channel */
         for(i=0; i<N; i++)
         {
         	AWGN_noise(0, noise_pwr, &noise[0]);
         	Yk[i] = bpsk[i] + noise[0];
         }

         for(i=0; i<K; i++)
				intrinsic[i] = 0.0;

         for(s=0; s<Iteration; s++)
         {
         	MAP_dec_BCJR(Yk, LLR, intrinsic, noise_pwr, K, num_state, num_in_sym);

				for(i=0; i<K; i++)	// data decision
         	{
         		if(LLR[i]>=0)
            		Ak[i] = 1;
            	else
            		Ak[i] = 0;

            	err_count[s] += Error_count(data_bit[i],Ak[i]);

               //intrinsic[i] = LLR[i] - intrinsic[i];
               intrinsic[i] = LLR[i];
         	}
         }
		}

      // Statistics and records
      cout << "Error Rate = ";
      for(s=0; s<Iteration; s++)
      {
      	err_rate[s] = err_count[s] / (double)(K*num_packet);
      	printf("%e, ", err_rate[s]);
      }
      cout << endl;

      fprintf(ber, "%f ", Eb_No);
      fprintf(records, "%f ", Eb_No);
		for(s=0; s<Iteration; s++)
      {
	      fprintf(ber, "%e ", err_rate[s]);
         fprintf(records, "%e ", err_rate[s]);
      }
      fprintf(ber, "\n");
      fprintf(records, "\n");
      fflush(records);
      fflush(ber);
   }

   delete data_bit;
   delete coded_bit;
   delete bpsk;
   delete Ak;
   delete Yk;
   delete LLR;
   delete intrinsic;
   delete err_count;
   delete err_rate;

   end = time(NULL);
   printf("Total elapsed time: %.0f(sec)\n", difftime(end,start));
   fprintf(records, "Total elapsed time: %.0f(sec)\n\n", difftime(end,start));
   fclose(ber);
   fclose(records);
   printf("This program is ended. Press any key to continue.\n");
   getch();

   return 0;
}

void Trellis_diagram(int n, int k, int m, int *gp)
{
/**************************************************************/
/* Generate TrellisDiagram for (n=2, k=1, m=2) {5,7} CC code */
/**************************************************************/
	int i, j, input, out_bit, out_sym, from_state, to_state, tmp_state;
   int num_state = pow(2,m), num_in_sym = pow(2,k);

   for(from_state=0; from_state<num_state; from_state++) // from_state of trellis diagram
   {
   	for(input=0; input<num_in_sym; input++)		// input of trellis diagram for (n, k, m)
      {
      	tmp_state = from_state;
         out_sym = 0;                      // output codeword symbol of trellis diagram
         tmp_state = (tmp_state << 1) ^ (input & 0x01);  // read input bit
         for(i=0; i<n; i++)
         {
         	out_bit = 0;						// output bit of trellis diagram
            for(j=m; j>=0; j--)
            	out_bit ^= ((tmp_state & gp[i]) >> j) & 1;  	// Calculate output bit

            out_sym = (out_sym << 1) ^ out_bit;				  	// Calculate output symbol
         }
         to_state = tmp_state & (num_state-1); 					// to_state of trellis diagram

         Trellis[from_state][input].from = from_state;
         Trellis[from_state][input].to = to_state;
         Trellis[from_state][input].in = input;
         Trellis[from_state][input].out[0] = ((out_sym>>1)&1);
         Trellis[from_state][input].out[1] = out_sym&1;
      }
   }
}

void CC_Encoder(int *Data_in, int *Code_bit, int n, int m, int *gp, int Packet_length, int Num_state)
{
/**********************************************************************/
/* Convolutional Encoder (n=2, k=1, m=2) Generator polynomial: {5, 7} */
/**********************************************************************/
   int i, j, l, from_state, tmp_state, out_bit;

	from_state = 0;
   for(i=0; i<Packet_length; i++)
   {
   	tmp_state = from_state;
//    out_sym = 0;                      // output codeword symbol of trellis diagram
		tmp_state = (tmp_state << 1) ^ (Data_in[i] & 0x01);  // read input bit
      for(j=0; j<n; j++)
      {
      	out_bit = 0;						// output bit of trellis diagram
         for(l=m; l>=0; l--)
         	out_bit ^= ((tmp_state & gp[j]) >> l) & 1;  	// Calculate output bit

         Code_bit[2*i+j] = out_bit;
//       out_sym = (out_sym << 1) ^ out_bit;				  	// Calculate output symbol
		}
      from_state = tmp_state & (Num_state-1); 					// to_state of trellis diagram
	}
}
/******************************************************************************/

void MAP_dec_BCJR(double *Yk, double *LLR, double *intrinsic,
						double noise_pwr, int Packet_length, int Num_state, int Num_in)
{
/**********************************************************************/
/* MAP decoder (modified BCJR Algorithm)                              */
/* Convolutional Decoder (n=2, k=1, m=2) Generator polynomial: {5, 7} */
/**********************************************************************/
	int i, j, l;
   double p1, p2, *normal, **alpha, **beta, ***gamma, *delta, min, **a_priori;

   a_priori = new double*[2];
   for(l=0; l<2; l++)
   	a_priori[l] = new double[Packet_length];
   alpha = new double*[Packet_length+1];			// alpha[time index][state]
   for(i=0; i<=Packet_length; i++)
      alpha[i] = new double[Num_state];
   beta = new double*[Packet_length+1];       	// beta[time index][state]
   for(i=0; i<=Packet_length; i++)
      beta[i] = new double[Num_state];
   gamma = new double**[Packet_length];			// gamma[time index][state][input]
   for(i=0; i<Packet_length; i++)
      gamma[i] = new double*[Num_state];
   for(i=0; i<Packet_length; i++)
   	for(j=0; j<Num_state; j++)
         gamma[i][j] = new double[Num_in];
   normal = new double[Packet_length+1];			// renormalization of BCJR
   delta = new double[Num_in];

	// Initialization of alpha and beta
   for(i=0; i<=Packet_length; i++)					// alpha[time index][state]
   	for(j=0; j<Num_state; j++)
      	alpha[i][j] = 0.0;
   alpha[0][0] = 1.0;

   for(i=0; i<=Packet_length; i++)           	// beta[time index][state]
   	for(j=0; j<Num_state; j++)
      	beta[i][j] = 0.0;
   beta[Packet_length][0] = 1.0;

   for(i=0; i<Packet_length; i++)
   	for(l=0; l<2; l++)
      	a_priori[l][i] = expl(l*intrinsic[i]) / (1 + expl(intrinsic[i]));

   // calculate gamma[time index][state][input]
   for(i=0; i<Packet_length; i++)					// time index
   	for(j=0; j<Num_state; j++)						// state index
      	for(l=0; l<Num_in; l++)						// input symbol
         {
            p1 = exp(-pow(Yk[2*i]-(2*Trellis[j][l].out[0]-1),2)/(2*noise_pwr))/sqrt(2*Pi*noise_pwr);
            p2 = exp(-pow(Yk[2*i+1]-(2*Trellis[j][l].out[1]-1),2)/(2*noise_pwr))/sqrt(2*Pi*noise_pwr);
            gamma[i][j][l] = a_priori[l][i] * p1 * p2;		// gamma[time index][state][input]
         }

   // calculate alpha[time index][state]
   for(i=1; i<=Packet_length; i++)					// time index
   {
   	for(j=0; j<Num_state; j++)						// from_state index
      	for(l=0; l<Num_in; l++)						// input bit
         	alpha[i][Trellis[j][l].to] += alpha[i-1][j] * gamma[i-1][j][l];

      normal[i] = 0.0;									// for renormalization
      for(j=0; j<Num_state; j++)						// to_state index
      	normal[i] += alpha[i][j];

      for(j=0; j<Num_state; j++)
      	alpha[i][j] = alpha[i][j] / normal[i];
   }
   normal[0] = 1.0;

   // calculate beta[time index][state]
   for(i=Packet_length-1; i>0; i--)					// time index
   {
   	for(j=0; j<Num_state; j++)						// from_state index
      	for(l=0; l<Num_in; l++)						// input bit
         	beta[i][j] += beta[i+1][Trellis[j][l].to] * gamma[i][j][l];

      for(j=0; j<Num_state; j++)
      	beta[i][j] = beta[i][j] / normal[i];
   }

   // Calculate conditional LLR of the information bit
   for(i=0; i<Packet_length; i++)					// time index
   {
   	min = 0.0;											// find the minimum product of alpha*gamma*beta
      for(j=0; j<Num_state; j++)
      	for(l=0; l<Num_in; l++)
         {
         	delta[0] = alpha[i][j] * gamma[i][j][l] * beta[i+1][Trellis[j][l].to];

            if((delta[0] < min && delta[0] != 0.0) || min == 0.0)
            	min = delta[0];
         }

      if(min == 0.0 || min > 1.0)					// if all else fails, make min real small
      	min = 1E-100;

      delta[0] = delta[1] = 0.0;
      for(j=0; j<Num_state; j++)						// from_state index
      	for(l=0; l<Num_in; l++)						// input bit
         	delta[l] += alpha[i][j] * gamma[i][j][l] * beta[i+1][Trellis[j][l].to];

      if(delta[1] == 0.0)
      	delta[1] = min;
      if(delta[0] == 0.0)
      	delta[0] = min;

      LLR[i] = logl(delta[1]/delta[0]);
   }

   for(l=0; l<2; l++)
   	delete a_priori[l];
   delete a_priori;
   for(i=0; i<=Packet_length; i++)
      delete alpha[i];
	delete alpha;
   for(i=0; i<=Packet_length; i++)
       delete beta[i];
   delete beta;
   for(i=0; i<Packet_length; i++)
   	for(j=0; j<Num_state; j++)
         delete gamma[i][j];
   for(i=0; i<Packet_length; i++)
      delete gamma[i];
   delete gamma;
   delete normal;
   delete delta;
}
/******************************************************************************/

void AWGN_noise(float mu, double variance, double *noise)
{
//	const  float Pi = 3.14159265358979;
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

