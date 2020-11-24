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

void AWGN_noise(float, double, double *);
int Error_count(int, int);

#define Pi 	3.14159265358979
#define n	2					// (n,k,m) Convolutional Code, # of output symbols
#define k	1					// (n,k,m) Convolutional Code, # of input symbols
#define m	2					// (n,k,m) Convolutional Code, # of memories
#define K	1000		  		// packet length of information bit stream
#define N	K*(n/k)	 		// packet length of coded bit stream
#define num_packet 10000		// number of packets simulated
const int num_state = pow(2,m);	// number of states
const int num_in_sym = pow(2,k);	// number of input symbols
const int gp[2] = {5,7};   // Generator polynomial of CC given in hexadecimal and in reverse order

int main(void)
{
	time_t  t, start, end;
	int i, j, l, p, input, out_sym, out_bit, *data_bit, *coded_bit;
   int from_state, to_state, tmp_state, *Ak, err_count;
   double a_priori=0.5, *bpsk, snr, Eb_No, noise_pwr, noise[2], *Yk, err_rate;
   double p1, p2, /*sum,*/ **alpha, **beta, ***gamma, *normal, *delta, *LLR, min;
   FILE *trelis,*ber, *records;
   struct trel_diag			// Data structure of trellis diagram for each branch
   		{
         	int from;		// from_state of trellis diagram
            int to;			// to_state of trellis diagram
            int in;			// input data bit of trellis diagram
            int out[n];		// output codeword symbol of trellis diagram
         };
   struct trel_diag Trellis[4][2];		// Trellis[num_state][num_in_sym]

   start = time(NULL);
   printf("BER Performance of Convolutional Code (2,1,2) in AWGN Channel\n");
	printf("Generator polynomials are {5,7} in Octal\n");
   printf("This program is running. Don't close, please!\n\n");

   data_bit = new int[K];
   coded_bit = new int[N];
   bpsk = new double[N];
   Ak = new int[K];
   Yk = new double[N];
   delta = new double[num_in_sym];
   LLR = new double[K];
   alpha = new double*[K+1];			  		// alpha[time index][state]
   for(i=0; i<=K; i++)
      alpha[i] = new double[num_state];
   beta = new double*[K+1];       			// beta[time index][state]
   for(i=0; i<=K; i++)
      beta[i] = new double[num_state];
   gamma = new double**[K];					// gamma[time index][state][input]
   for(i=0; i<K; i++)
      gamma[i] = new double*[num_state];
   for(i=0; i<K; i++)
   	for(j=0; j<num_state; j++)
         gamma[i][j] = new double[num_in_sym];
   normal = new double[K+1];					// for normalization of alpha & beta

	srand((unsigned) time(&t));

/**************************************************************/
/* Generate TrellisDiagram for (n=2, k=1, m=2) {5,7} CC code */
/**************************************************************/
	trelis = fopen("Trellis_diagram.log", "w");
   fprintf(trelis, "Trellis diagram of (%d,%d,%d) convolutional code\n", n, k, m);
   fprintf(trelis, "Generator polynomials are {%d,%d}\n\n", gp[0], gp[1]);
   fprintf(trelis, "s(k-1) s(k) input out5 out7\n");
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
         fprintf(trelis, "%4d %4d %5d %5d %4d\n", from_state, to_state, input, ((out_sym>>1)&1), out_sym&1);
      }
   }
   fclose(trelis);

/************************/
/* main simulation loop */
/************************/
   for(snr=0; snr<=10; snr++)
   {
   	err_count = 0;
   	// noise power calculation
      Eb_No = (double)snr;
      noise_pwr = 1.0/(((float)k/(float)n)*pow(10.0, Eb_No/10.0));	// BPSK, Nyquist filter assumption
      printf("Eb_No = %f, ", Eb_No);
      for(p=0; p<num_packet; p++)
      {
/*   		for(i=0; i<K-2; i++)		// Generate random information bit stream
   		{
   			if(rand()/RAND_MAX>=0.5)
      			data_bit[i] = 1;
      		else
      			data_bit[i] = 0;
   		}
*/
   		for(i=0; i<K-2; i++)
				data_bit[i] = random(2);		// Generate random information bit stream
   		data_bit[K-2] = data_bit[K-1] = 0;

/**********************************************************************/
/* Convolutional Encoder (n=2, k=1, m=2) Generator polynomial: {5, 7} */
/**********************************************************************/
   		from_state = 0;
   		for(i=0; i<K; i++)
   		{
				tmp_state = from_state;
//      		out_sym = 0;                      // output codeword symbol of trellis diagram
      		tmp_state = (tmp_state << 1) ^ (data_bit[i] & 0x01);  // read input bit
      		for(j=0; j<n; j++)
      		{
      			out_bit = 0;						// output bit of trellis diagram
         		for(l=m; l>=0; l--)
         			out_bit ^= ((tmp_state & gp[j]) >> l) & 1;  	// Calculate output bit

         		coded_bit[2*i+j] = out_bit;
//         		out_sym = (out_sym << 1) ^ out_bit;				  	// Calculate output symbol
      		}
      		from_state = tmp_state & (num_state-1); 					// to_state of trellis diagram
   		}

/**********************************************************************/

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

/*******************************************/
/* MAP decoder (modified BCJR Algorithm) */
/*******************************************/
			// Initialization of alpha and beta
         for(i=0; i<=K; i++)				// alpha[time index][state]
         	for(j=0; j<num_state; j++)
         		alpha[i][j] = 0.0;
         alpha[0][0] = 1.0;

         for(i=0; i<=K; i++)           // beta[time index][state]
         	for(j=0; j<num_state; j++)
         		beta[i][j] = 0.0;
         beta[K][0] = 1.0;

			// calculate gamma[time index][state][input]
         for(i=0; i<K; i++)	// time index
         	for(j=0; j<num_state; j++)		// state index
            	for(l=0; l<num_in_sym; l++)	// input symbol
               {
						p1 = exp(-pow(Yk[2*i]-(2*Trellis[j][l].out[0]-1),2)/(2*noise_pwr))/sqrt(2*Pi*noise_pwr);
                  p2 = exp(-pow(Yk[2*i+1]-(2*Trellis[j][l].out[1]-1),2)/(2*noise_pwr))/sqrt(2*Pi*noise_pwr);
         			gamma[i][j][l] = a_priori * p1 * p2;		// gamma[time index][state][input]
               }

         // calculate alpha[time index][state]
         for(i=1; i<=K; i++)		// time index
         {
         	for(j=0; j<num_state; j++)		// from_state index
            	for(l=0; l<num_in_sym; l++)	// input bit
         			alpha[i][Trellis[j][l].to] += alpha[i-1][j] * gamma[i-1][j][l];

//            sum = 0.0;		// for renormalization
				normal[i] = 0.0;
            for(j=0; j<num_state; j++)		// to_state index
            	normal[i] += alpha[i][j];

            for(j=0; j<num_state; j++)
            	alpha[i][j] = alpha[i][j] / normal[i];
         }
         normal[0] = 1.0;

         // calculate beta[time index][state]
         for(i=K-1; i>0; i--)		// time index
         {
         	for(j=0; j<num_state; j++)		// from_state index
            	for(l=0; l<num_in_sym; l++)	// input bit
         			beta[i][j] += beta[i+1][Trellis[j][l].to] * gamma[i][j][l];

/*            sum = 0.0;		// for renormalization
            for(j=0; j<num_state; j++)		// from_state index
            	sum += beta[i][j];
*/
            for(j=0; j<num_state; j++)
            	beta[i][j] = beta[i][j] / normal[i];
         }

         // calculate conditional LLR
         for(i=0; i<K; i++)		// time index
         {
         	min = 0.0;		// find the minimum product of alpha*gamma*beta
            for(j=0; j<num_state; j++)
            	for(l=0; l<num_in_sym; l++)
               {
                  delta[0] = alpha[i][j] * gamma[i][j][l] * beta[i+1][Trellis[j][l].to];

                  if((delta[0] < min && delta[0] != 0.0) || min == 0.0)
                  	min = delta[0];
					}

            if(min == 0.0 || min > 1.0)	// if all else fails, make min real small
            	min = 1E-100;

         	delta[0] = delta[1] = 0.0;
         	for(j=0; j<num_state; j++)		// from_state index
            	for(l=0; l<num_in_sym; l++)	// input bit
               	delta[l] += alpha[i][j] * gamma[i][j][l] * beta[i+1][Trellis[j][l].to];

            if(delta[1] == 0.0)
            	delta[1] = min;
//            else if(delta[0] == 0.0)
            if(delta[0] == 0.0)
            	delta[0] = min;

            LLR[i] = log(delta[1]/delta[0]);
         }

/*******************************************/

			for(i=0; i<K; i++)	// data decision
         {
         	if(LLR[i]>=0)
            	Ak[i] = 1;
            else
            	Ak[i] = 0;

            err_count += Error_count(data_bit[i],Ak[i]);
         }
		}

      // Statistics and records
      err_rate = err_count / (double)(K*num_packet);
      printf("Error rate = %e\n", err_rate);
      ber=fopen("ber.log","a");
   	records = fopen("cc_MAP_awgn.log", "a");
      fprintf(ber, "%f %e\n", Eb_No, err_rate);
      fprintf(records, "(2,1,2) {5,7} convolutional code, MAP decoder, AWGN channel with Eb/N0 = %f dB\n", Eb_No);
   	fprintf(records, "Average bit error rate = %e\n", err_rate);
   	fprintf(records, "number of bits of simulation = %d\n\n", K*num_packet);

      fclose(ber);
	   fclose(records);
   }

   delete data_bit;
   delete coded_bit;
   delete bpsk;
   delete Ak;
   delete Yk;
   delete delta;
   delete LLR;
   delete alpha;
   delete beta;
   delete gamma;
   delete normal;

   records = fopen("cc_MAP_awgn.log", "a");
   end = time(NULL);
   printf("Total elapsed time: %.0f(sec)\n", difftime(end,start));
   fprintf(records, "Total elapsed time: %.0f(sec)\n\n", difftime(end,start));
   fclose(records);
   printf("This program is ended. Press any key to continue.\n");
   getch();

   return 0;
}

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

