/*********************************************************/
/* Author: Chao-wang Huang                               */
/* Date: Thursday, November 24, 2005                     */
/* An (n,k,m) Convolutional code is simulated            */
/* Decoding algorithm: Viterbi Algorithm (Soft Decision) */
/* Modulation Scheme: 16-QAM                             */
/*********************************************************/

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

#define Pi 	3.14159265358979
#define n	2					// (n,k,m) Convolutional Code, # of output symbols
#define k	1					// (n,k,m) Convolutional Code, # of input symbols
#define m	2					// (n,k,m) Convolutional Code, # of memories
#define K	1000		  		// packet length of information bit stream
#define N	K*(n/k)	 		// packet length of coded bit stream
#define num_packet 10000		// number of packets simulated
const int num_state = pow(2,m);	// number of states
#define N_state  4
const int num_in_sym = pow(2,k);	// number of input symbols
const int gp[2] = {5,7};   // Generator polynomial of CC given in hexadecimal and in reverse order

int main(void)
{
	time_t  t, start, end;
	int i, j, l, p, q, input, out_sym, out_bit, *data_bit, *coded_bit;
   int from_state, to_state, tmp_state, err_count;
   double snr, Eb_No, noise_pwr, noise[2], *Yk_I, *Yk_Q, err_rate;
   double *survival_metric, metric, *sym_I, *sym_Q, *Lamda, *Dk;
   FILE *trelis,*ber, *records;
   struct trel_diag			// Data structure of trellis diagram for each branch
   		{
         	int from;		// from_state of trellis diagram
            int to;			// to_state of trellis diagram
            int in;			// input data bit of trellis diagram
            int out[n];		// output codeword symbol of trellis diagram
         };
   struct trel_diag Trellis[4][2];		// Trellis[num_state][num_in_sym]

   struct surv        			// Data structure of survival path
   		{
         	double metric;		// Path metric
            int data_in[K];	// input bit stream
            int state[K];		// state transition sequence
         };
   struct surv survival_path[N_state], survival_temp[N_state];	// Survival_path[num_state]

   start = time(NULL);
   printf("BER Performance of Convolutional Code (2,1,2) in AWGN Channel\n");
	printf("Generator polynomials are {5,7} in Octal\n");
	cout << "Minimum Free Distance: 5" << endl;
   printf("This program is running. Don't close, please!\n\n");

   records = fopen("Records_VA_awgn.log", "a");
   fprintf(records, "(2,1,2) {5,7} convolutional code, Viterbi decoder, AWGN channel\n");
   fprintf(records, "number of bits of simulation = %d\n\n", K*num_packet);
   fprintf(records, "Eb/No     BER\n");
   fflush(records);

   data_bit = new int[K];
   coded_bit = new int[N];
   sym_I = new double[N/4];
   sym_Q = new double[N/4];
   Yk_I = new double[N/4];
   Yk_Q = new double[N/4];
   Lamda = new double[N];
   Dk = new double[N];
   survival_metric = new double[num_state];

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
   ber=fopen("ber_SDVA_16QAM.log","w");
   for(snr=0; snr<=15; snr++)
   {
   	err_count = 0;
   	// noise power calculation
      Eb_No = (double)snr;
      noise_pwr = 0.25/(((float)k/(float)n)*pow(10.0, Eb_No/10.0));	// 16-QAM, Nyquist filter assumption
      printf("Eb_No = %f, ", Eb_No);
      for(p=0; p<num_packet; p++)
      {
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

         // 16-QAM Mapping
         for(i=0; i<N/4; i++)
         {
	   	   if(coded_bit[4*i] == 0 && coded_bit[4*i+2] == 0)
   	   		sym_I[i] = 1.0/sqrt(10.0);
	   	   else if(coded_bit[4*i] == 0 && coded_bit[4*i+2] == 1)
   	   		sym_I[i] = 3.0/sqrt(10.0);
	      	else if(coded_bit[4*i] == 1 && coded_bit[4*i+2] == 1)
   	   		sym_I[i] = -3.0/sqrt(10.0);
	   	   else
   	   		sym_I[i] = -1.0/sqrt(10.0);

	      	if(coded_bit[4*i+1] == 0 && coded_bit[4*i+3] == 0)
   	   		sym_Q[i] = 1.0/sqrt(10.0);
	   	   else if(coded_bit[4*i+1] == 0 && coded_bit[4*i+3] == 1)
   	   		sym_Q[i] = 3.0/sqrt(10.0);
	      	else if(coded_bit[4*i+1] == 1 && coded_bit[4*i+3] == 1)
   	   		sym_Q[i] = -3.0/sqrt(10.0);
	   	   else
   	   		sym_Q[i] = -1.0/sqrt(10.0);
         }

         /* AWGN channel */
         for(i=0; i<N/4; i++)
         {
         	AWGN_noise(0, noise_pwr, &noise[0]);
         	Yk_I[i] = sym_I[i] + noise[0];
            Yk_Q[i] = sym_Q[i] + noise[1];

            // Soft Decision
            Lamda[4*i] = -log((exp(2*(1.0/sqrt(10.0))*Yk_I[i]/noise_pwr)+exp(((6*(1.0/sqrt(10.0))*Yk_I[i])-(8*0.1))/noise_pwr))
                         /(exp(-2*(1.0/sqrt(10.0))*Yk_I[i]/noise_pwr)+exp((-6*(1.0/sqrt(10.0))*Yk_I[i]-(8*0.1))/noise_pwr)));
   			Lamda[4*i+1] = -log((exp(2*(1.0/sqrt(10.0))*Yk_Q[i]/noise_pwr)+exp(((6*(1.0/sqrt(10.0))*Yk_Q[i])-(8*0.1))/noise_pwr))
                           /(exp(-2*(1.0/sqrt(10.0))*Yk_Q[i]/noise_pwr)+exp((-6*(1.0/sqrt(10.0))*Yk_Q[i]-(8*0.1))/noise_pwr)));
            Lamda[4*i+2] = -log((exp(2*(1.0/sqrt(10.0))*Yk_I[i]/noise_pwr)+exp(-2*(1.0/sqrt(10.0))*Yk_I[i]/noise_pwr))
                           /(exp(((6*(1.0/sqrt(10.0))*Yk_I[i])-(8*0.1))/noise_pwr)+exp((-6*(1.0/sqrt(10.0))*Yk_I[i]-(8*0.1))/noise_pwr)));
				Lamda[4*i+3] = -log((exp(2*(1.0/sqrt(10.0))*Yk_Q[i]/noise_pwr)+exp(-2*(1.0/sqrt(10.0))*Yk_Q[i]/noise_pwr))
                           /(exp(((6*(1.0/sqrt(10.0))*Yk_Q[i])-(8*0.1))/noise_pwr)+exp((-6*(1.0/sqrt(10.0))*Yk_Q[i]-(8*0.1))/noise_pwr)));

            //printf("%2.2f %2.2f %2.2f %2.2f\n", Lamda[4*i], Lamda[4*i+1], Lamda[4*i+2], Lamda[4*i+3]);
            //getch();

         }
/*
         // 16-QAM De-mapping (Hard Decision)
         for(i=0; i<N/4; i++)
         {
   	   	if(Yk_I[i] < -2.0/sqrt(10.0))
      		{
      			Dk[4*i] = 1;
		         Dk[4*i+2] = 1;
   		   }
      		else if(Yk_I[i] >= -2.0/sqrt(10.0) && Yk_I[i] < 0.0)
		      {
   		   	Dk[4*i] = 1;
      		   Dk[4*i+2] = -1;
	      	}
	   	   else if(Yk_I[i] >= 0.0 && Yk_I[i] < 2.0/sqrt(10.0))
   	   	{
      			Dk[4*i] = -1;
	      	   Dk[4*i+2] = -1;
	   	   }
   	   	else
	   	   {
   	   		Dk[4*i] = -1;
      	   	Dk[4*i+2] = 1;
		      }

   		   if(Yk_Q[i] < -2.0/sqrt(10.0))
      		{
      			Dk[4*i+1] = 1;
		         Dk[4*i+3] = 1;
   		   }
      		else if(Yk_Q[i] >= -2.0/sqrt(10.0) && Yk_Q[i] < 0.0)
	      	{
   	   		Dk[4*i+1] = 1;
	      	   Dk[4*i+3] = -1;
		      }
   		   else if(Yk_Q[i] >= 0.0 && Yk_Q[i] < 2.0/sqrt(10.0))
      		{
      			Dk[4*i+1] = -1;
		         Dk[4*i+3] = -1;
   		   }
      		else
	      	{
   	   		Dk[4*i+1] = -1;
	      	   Dk[4*i+3] = 1;
		      }
         }
*/
         // Soft Bit
         for(i=0; i<N; i++)
         {
				Dk[i] = tanh(0.5*Lamda[i]);
            //printf("%f\n", Dk[i]);
            //getch();
         }


/*********************************************/
/* Viterbi Algorithm decoder (Soft Decision) */
/*********************************************/
         // Initialize survival path
         for(i=0; i<num_state; i++)
   		{
//         	survival_path[i].metric = DBL_MAX;	// Initial maximum value for Euclidean distance
   			survival_path[i].metric = -DBL_MAX;	// Initial minimum value for cross-correlation
      		for(j=0; j<K; j++)
      		{
      			survival_path[i].data_in[j] = 0;
         		survival_path[i].state[j] = 0;
      		}
   		}
         survival_path[0].metric = 0.0;

         for(i=0; i<K; i++)
         {
         	for(j=0; j<num_state; j++)		// Initialize the survival path metric
//         		survival_metric[j] = DBL_MAX;		// Initial maximum value for Euclidean distance
               survival_metric[j] = -DBL_MAX;	// Initial minimum value for cross-correlation

         	for(j=0; j<num_state; j++)		// from_state index
            	for(l=0; l<num_in_sym; l++)	// input bit
               {
/*               	metric = 0.0;		// branch metric, Euclidean Distance
                  metric += pow(Dk[2*i]-(2*Trellis[j][l].out[0]-1),2);		// brahch metric
                  metric += pow(Dk[2*i+1]-(2*Trellis[j][l].out[1]-1),2);	// branch metric
                  metric += survival_path[j].metric;
*/
                  metric = 0.0;		// branch metric, Cross-correlation
                  metric += (Dk[2*i] * (2*Trellis[j][l].out[0]-1));		// brahch metric
                  metric += (Dk[2*i+1] * (2*Trellis[j][l].out[1]-1));	// branch metric
                  metric += survival_path[j].metric;

                  // find the survival path metric
//                  if(metric < survival_metric[Trellis[j][l].to])	//	Euclidean distance (Minimize)
                  if(metric > survival_metric[Trellis[j][l].to])	// Cross-correlation (Maximize)
                  {
                  	survival_metric[Trellis[j][l].to] = metric;

                  	// Record and refresh the survival path
                  	for(q=0; q<i; q++)
                  	{
                  		survival_temp[Trellis[j][l].to].data_in[q] = survival_path[j].data_in[q];
                     	survival_temp[Trellis[j][l].to].state[q] = survival_path[j].state[q];
                  	}
                  	survival_temp[Trellis[j][l].to].data_in[i] = l;
                  	survival_temp[Trellis[j][l].to].state[i] = Trellis[j][l].to;
                  }
               }

            // Record and refresh the survival path
            for(j=0; j<num_state; j++)		// to_state index
            {
            	survival_path[j].metric = survival_metric[j];
               for(q=0; q<=i; q++)
               {
               	survival_path[j].data_in[q] = survival_temp[j].data_in[q];
                  survival_path[j].state[q] = survival_temp[j].state[q];
               }
            }
         }

/*******************************************/

			for(i=0; i<K; i++)	// Bit error count
            err_count += Error_count(data_bit[i], survival_path[0].data_in[i]);
		}

      // Statistics and records
      err_rate = err_count / (double)(K*num_packet);
      printf("Error rate = %e\n", err_rate);
      fprintf(ber, "%f %e\n", Eb_No, err_rate);
      fprintf(records, "%f %e\n", Eb_No, err_rate);
      fflush(records);
   }

   fclose(ber);
   delete data_bit;
   delete coded_bit;
   delete sym_I;
   delete sym_Q;
   delete Yk_I;
   delete Yk_Q;
   delete Lamda;
   delete Dk;
   delete survival_metric;

   end = time(NULL);
   printf("Total elapsed time: %.0f(sec)\n", difftime(end,start));
   fprintf(records, "\nTotal elapsed time: %.0f(sec)\n\n", difftime(end,start));
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

