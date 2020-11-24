/*========================================================================================*/
/* Author: Chao-wang Huang                                                                */
/* Date: Monday, July 31, 2006                                                            */
/* An (2,1,6) Tail-Biting Convolutional encoder is simulated                              */
/* Initialize the encoder's memory with the last data bits of the FEC block being encoded */
/* Decoding algorithm: modified BCJR algorithm                                            */
/*========================================================================================*/

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
void MAP_dec_BCJR(double *, double *, double *, int, int, int, double);

#define Pi 	3.14159265358979
//#define n	2							// (n,k,m) Convolutional Code, # of output symbols
const int n = 2;
//#define k	1							// (n,k,m) Convolutional Code, # of input symbols
const int k = 1;
//#define m	2							// (n,k,m) Convolutional Code, # of memories
const int m = 6;
#define K	1000		  				// packet length of information bit stream
#define N	K*(n/k)	 				// packet length of coded bit stream
#define num_packet 100000				// number of packets simulated
const int num_state = pow(2,m);	// number of states
#define N_state  64
const int num_in_sym = pow(2,k);	// number of input symbols
//const int gp[2] = {171,133};	// Generator polynomial of CC given in Octal
int gp[2] = {79,109};           	// Generator polynomial of CC given in Decimal and in reverse order

struct trel_diag						// Data structure of trellis diagram for each branch
		{
      	int from;					// from_state of trellis diagram
         int to;						// to_state of trellis diagram
         int in;						// input data bit of trellis diagram
         int out[n];					// output codeword symbol of trellis diagram
      };
struct trel_diag Trellis[64][2];	// Trellis[num_state][num_in_sym]

int main(void)
{
	time_t  t, start, end;
	int i, /*j,*/ p, *data_bit, *coded_bit, err_count, *Ak;
   double *bpsk, snr, Eb_No, noise_pwr, noise[2], *Yk, err_rate, *Dk, *LLR, *Hk;
   FILE /*trelis,*/ *ber, *records;

   start = time(NULL);
   printf("BER Performance of Convolutional Code (2,1,6) in AWGN Channel\n");
	printf("Generator polynomials are {171,133} in Octal\n");
	cout << "Minimum Free Distance: 10" << endl;
   printf("Maximum number of bits of simulation = %d\n\n", K*num_packet);
   printf("This program is running. Don't close, please!\n\n");

   records = fopen("Records_BCJR_awgn.log", "a");
   fprintf(records, "(2,1,6) {171,133} Octal convolutional code, BCJR decoder, AWGN channel\n");
   fprintf(records, "Minimum Free Distance: 10\n");
   fprintf(records, "Maximum number of bits of simulation = %d\n\n", K*num_packet);
   fprintf(records, "Eb/No     BER\n");
   fflush(records);

   data_bit = new int[K];
   Ak = new int[K];
   Hk = new double[N];
   coded_bit = new int[N];
   bpsk = new double[N];
   Yk = new double[N];			// Received signal
   Dk = new double[N];			// Soft Decision
   LLR = new double[K];			// Log-likelihood Ratio

	srand((unsigned) time(&t));

   Trellis_diagram(n, k, m, &gp[0]);
/*
	trelis = fopen("Trellis_diagram.log", "w");
   fprintf(trelis, "Trellis diagram of (%d,%d,%d) convolutional code\n", n, k, m);
   fprintf(trelis, "Generator polynomials are {171,133} in Octal\n");
   fprintf(trelis, "Generator polynomials are {%d,%d} in Decimal reverse order\n\n", gp[0], gp[1]);
   fprintf(trelis, "s(k-1) s(k) input out171 out133\n");
   for(i=0; i<num_state; i++) // from_state of trellis diagram
   	for(j=0; j<num_in_sym; j++)		// input of trellis diagram for (n, k, m)
         fprintf(trelis, "%4d %4d %5d %5d %4d\n", Trellis[i][j].from, Trellis[i][j].to, Trellis[i][j].in, Trellis[i][j].out[0], Trellis[i][j].out[1]);
   fclose(trelis);
  */
/************************/
/* main simulation loop */
/************************/
   ber=fopen("ber.log","w");
   for(snr=0; snr<=5; snr++)
   {
   	err_count = 0;
   	// noise power calculation
      Eb_No = (double)snr;
      noise_pwr = 1.0/(((float)k/(float)n)*pow(10.0, Eb_No/10.0));	// BPSK, Nyquist filter assumption
      printf("Eb_No = %f, ", Eb_No);

      p = 0;
      do
      //for(p=0; p<num_packet; p++)
      {
   		for(i=0; i<K-m; i++)
				data_bit[i] = random(2);		// Generate random information bit stream
         for(i=K-m; i<K; i++)
         	data_bit[i] = 0;

         // Tail-Biting Convolutional Encoder
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
            // Soft Decision (LLR)
            //Dk[i] = sqrt(2.0)*Yk[i]/noise_pwr;
				//Dk[i] = sqrt(2.0)*Yk[i];
            Dk[i] = Yk[i];
         }

         MAP_dec_BCJR(Dk, LLR, Hk, K, N_state, num_in_sym, noise_pwr);

         // Hard Decision
         for(i=0; i<K; i++)
         {
         	if(LLR[i]>=0)
            	Ak[i] = 1;
            else
            	Ak[i] = 0;
         }

         // Bit error count
         for(i=0; i<K-m; i++)
            err_count += Error_count(data_bit[i], Ak[i]);

         p++;
		} while(err_count <= 1000 && p < num_packet);

      // Statistics and records
      err_rate = err_count / (double)((K-m)*p);
      printf("Error rate = %e\n", err_rate);
      fprintf(ber, "%f %e\n", Eb_No, err_rate);
      fprintf(records, "%f %e\n", Eb_No, err_rate);
      fflush(records);
      fflush(ber);
   }

   fclose(ber);
   delete data_bit;
   delete Ak;
   delete Hk;
   delete coded_bit;
   delete bpsk;
   delete Yk;
   delete Dk;
   delete LLR;

   end = time(NULL);
   printf("Total elapsed time: %.0f(sec)\n", difftime(end,start));
   fprintf(records, "\nTotal elapsed time: %.0f(sec)\n\n", difftime(end,start));
   fclose(records);
   printf("This program is ended. Press any key to continue.\n");
   getch();

   return 0;
}

void Trellis_diagram(int n, int k, int m, int *gp)
{
/******************************************************************/
/* Generate TrellisDiagram for (n=2, k=1, m=6) {171,133} CC code */
/******************************************************************/
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
/******************************************************************************************/
/* Tail-Biting Convolutional Encoder (n=2, k=1, m=6) Generator polynomial: {171, 133}     */
/* Initialize the encoder's memory with the last data bits of the FEC block being encoded */
/******************************************************************************************/
   int i, j, l, from_state, tmp_state, out_bit;

   from_state = 0;
   //for(i=0; i<m; i++)
	  //	from_state += Data_in[Packet_length-i-1] * pow(2,i);

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

void MAP_dec_BCJR(double *Received/*intrinsic*/, double *LLR, double *extrinsic, int Packet_length, int Num_state, int Num_in, double noise_pwr)
{
/*========================================================================================*/
/* A Tail-Biting Convolutional encoder is simulated                                       */
/* Initialize the encoder's memory with the last data bits of the FEC block being encoded */
/* MAP decoder (modified BCJR Algorithm)                                                  */
/* Convolutional Decoder (n=2, k=1, m=6) Generator polynomial: {171, 133}                 */
/*========================================================================================*/
	int i, j, l;
   double p1, p2;
   double a_priori=0.5, *normal, **alpha, **beta, ***gamma, *delta, min, **code_bit_prob;
   	/* Information Bit */

   code_bit_prob = new double*[2*Packet_length];// a-priori probability, a-priori[time index][input]
   for(i=0; i<2*Packet_length; i++)
   	code_bit_prob[i] = new double[Num_in];
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
/*
   // Calculate a-priori probability of the code bit
   for(i=0; i<2*Packet_length; i++)					// time index
   	for(l=0; l<Num_in; l++)							// code bit
      	code_bit_prob[i][l] = exp(l*intrinsic[i]) / (1 + exp(intrinsic[i]));
  */
   // calculate gamma[time index][state][input]
   for(i=0; i<Packet_length; i++)					// time index
   	for(j=0; j<Num_state; j++)						// state index
      	for(l=0; l<Num_in; l++)						// input symbol
         {
         	p1 = exp(-pow(Received[2*i]-(2*Trellis[j][l].out[0]-1),2)/(2*noise_pwr))/sqrt(2*Pi*noise_pwr);
            p2 = exp(-pow(Received[2*i+1]-(2*Trellis[j][l].out[1]-1),2)/(2*noise_pwr))/sqrt(2*Pi*noise_pwr);
            gamma[i][j][l] = a_priori * p1 * p2;		// gamma[time index][state][input]
         	//gamma[i][j][l] = a_priori * code_bit_prob[2*i][Trellis[j][l].out[0]]
              //       				     * code_bit_prob[2*i+1][Trellis[j][l].out[1]];
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

      LLR[i] = log(delta[1]/delta[0]);

	   // Calculate Extrinsic Information of the 1st code bit
		delta[0] = delta[1] = 0.0;
      for(j=0; j<Num_state; j++)						// from_state index
      	for(l=0; l<Num_in; l++)                // input bit
         	delta[Trellis[j][l].out[0]] += alpha[i][j] * gamma[i][j][l]
                     										 	 * beta[i+1][Trellis[j][l].to];

      if(delta[1] == 0.0)
      	delta[1] = min;
      if(delta[0] == 0.0)
      	delta[0] = min;

      extrinsic[2*i] = log(delta[1]/delta[0]);	// calculate extrinsic information

      // Calculate Extrinsic Information of the 2nd code bit
      delta[0] = delta[1] = 0.0;
      for(j=0; j<Num_state; j++)						// from_state index
      	for(l=0; l<Num_in; l++)                // input bit
         	delta[Trellis[j][l].out[1]] += alpha[i][j] * gamma[i][j][l]
                     										 	 * beta[i+1][Trellis[j][l].to];

      if(delta[1] == 0.0)
      	delta[1] = min;
      if(delta[0] == 0.0)
      	delta[0] = min;

      extrinsic[2*i+1] = log(delta[1]/delta[0]);// calculate extrinsic information
   }

	for(i=0; i<2*Packet_length; i++)
   	delete code_bit_prob[i];
   delete code_bit_prob;
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

