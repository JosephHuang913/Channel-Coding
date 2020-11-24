/************************************************************/
/* Author: Chao-wang Huang                                  */
/* Date: Monday, July 10, 2006                              */
/* An (2,1,4) Convolutional code is simulated               */
/* Decoding algorithm: Soft Output Viterbi Algorithm (SOVA) */
/************************************************************/

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
void Viterbi_dec_SOVA(double *, int *, double *, int, int, int);

#define Pi 	3.14159265358979
//#define n	2					// (n,k,m) Convolutional Code, # of output symbols
const int n = 2;
//#define k	1					// (n,k,m) Convolutional Code, # of input symbols
const int k = 1;
//#define m	2					// (n,k,m) Convolutional Code, # of memories
const int m = 4;
#define K	1000		  		// packet length of information bit stream
#define N	K*(n/k)	 		// packet length of coded bit stream
#define num_packet 100000		// number of packets simulated
const int num_state = pow(2,m);	// number of states
#define N_state  16
const int num_in_sym = pow(2,k);	// number of input symbols
//const int gp[2] = {23,35};		// Generator polynomial of CC given in Octal
int gp[2] = {25,23};           // Generator polynomial of CC given in Decimal and in reverse order

struct trel_diag			// Data structure of trellis diagram for each branch
		{
      	int from;		// from_state of trellis diagram
         int to;			// to_state of trellis diagram
         int in;			// input data bit of trellis diagram
         int out[n];		// output codeword symbol of trellis diagram
      };
struct trel_diag Trellis[16][2];		// Trellis[num_state][num_in_sym]

int main(void)
{
	time_t  t, start, end;
	int i, j, p, *data_bit, *coded_bit, err_count[2], *Ak, *Hk;
   double *bpsk, snr, Eb_No, noise_pwr, noise[2], *Yk, err_rate[2], *Dk, *LLR;
   FILE *trelis, *ber, *records;

   start = time(NULL);
   printf("BER Performance of Convolutional Code (2,1,4) in AWGN Channel\n");
	printf("Generator polynomials are {23,35} in Octal\n");
	cout << "Minimum Free Distance: 7" << endl;
   printf("Maximum number of bits of simulation = %d\n\n", K*num_packet);
   printf("This program is running. Don't close, please!\n\n");

   records = fopen("Records_SOVA_awgn.log", "a");
   fprintf(records, "(2,1,4) {23,35} Octal convolutional code, soft output Viterbi decoder, AWGN channel\n");
   fprintf(records, "Minimum Free Distance: 7\n");
   fprintf(records, "Maximum number of bits of simulation = %d\n\n", K*num_packet);
   fprintf(records, "Eb/No     BER\n");
   fflush(records);

   data_bit = new int[K];
   Ak = new int[K];
   Hk = new int[K];
   coded_bit = new int[N];
   bpsk = new double[N];
   Yk = new double[N];		// Received signal
   Dk = new double[N];			// Soft Decision
   LLR = new double[K];		// Log-likelihood Ratio

	srand((unsigned) time(&t));

   Trellis_diagram(n, k, m, &gp[0]);

	trelis = fopen("Trellis_diagram.log", "w");
   fprintf(trelis, "Trellis diagram of (%d,%d,%d) convolutional code\n", n, k, m);
   fprintf(trelis, "Generator polynomials are {23,35} in Octal\n");
   fprintf(trelis, "Generator polynomials are {%d,%d} in Decimal reverse order\n\n", gp[0], gp[1]);
   fprintf(trelis, "s(k-1) s(k) input out23 out35\n");
   for(i=0; i<num_state; i++) // from_state of trellis diagram
   	for(j=0; j<num_in_sym; j++)		// input of trellis diagram for (n, k, m)
         fprintf(trelis, "%4d %4d %5d %5d %4d\n", Trellis[i][j].from, Trellis[i][j].to, Trellis[i][j].in, Trellis[i][j].out[0], Trellis[i][j].out[1]);
   fclose(trelis);

/************************/
/* main simulation loop */
/************************/
   ber=fopen("ber.log","w");
   for(snr=0; snr<=10; snr++)
   {
   	err_count[0] = err_count[1] = 0;
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
            // Soft Decision
//            Dk[i] = tanh(Yk[i]/noise_pwr);
            Dk[i] = Yk[i];
         }

         Viterbi_dec_SOVA(Dk, Hk, LLR, K, num_state, num_in_sym);

         // Hard Decision
         for(i=0; i<K; i++)
         {
         	if(LLR[i]>=0)
            	Ak[i] = 1;
            else
            	Ak[i] = 0;

			  	// Bit error count
            err_count[0] += Error_count(data_bit[i], Ak[i]);
            err_count[1] += Error_count(data_bit[i], Hk[i]);
         }

         p++;
		} while(err_count[1] <= 1000 && p < num_packet);

      // Statistics and records
      err_rate[0] = err_count[0] / (double)(K*p);
      err_rate[1] = err_count[1] / (double)(K*p);
      printf("Error rate = %e %e\n", err_rate[0], err_rate[1]);
      fprintf(ber, "%f %e %e\n", Eb_No, err_rate[0], err_rate[1]);
      fprintf(records, "%f %e %e\n", Eb_No, err_rate[0], err_rate[1]);
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
/****************************************************************/
/* Generate TrellisDiagram for (n=2, k=1, m=4) {23,35} CC code */
/****************************************************************/
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
/************************************************************************/
/* Convolutional Encoder (n=2, k=1, m=4) Generator polynomial: {23, 35} */
/************************************************************************/
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

void Viterbi_dec_SOVA(double *Data_in, int *Data_out, double *Soft_out, int Packet_length, int Num_state, int Num_in)
{
/**********************************************************************/
/* Soft Output Viterbi Algorithm decoder (SOVA)                       */
/* Convolutional Decoder (n=2, k=1, m=2) Generator polynomial: {5, 7} */
/**********************************************************************/
	int i, j, l, q, pre_state;
   double **mju_f, *mju, *survival_metric, metric, ***branch_metric, **mju_b, mju_tmp;

   struct surv        			// Data structure of survival path
   		{
         	double metric;		// Path metric
            int data_in[K];	// input bit stream, K: packet length, global
            int code_bit[N];	// output code bit stream, N: packet length, global
            int state[K];		// state transition sequence, K: packet length, global
         };
   struct surv survival_path[N_state], survival_temp[N_state];
   									// Survival_path[N_state], N_state: number of states of CC code, global

   survival_metric = new double[Num_state];
   mju = new double[Num_in];							// minimum path metric
   mju_f = new double*[Packet_length+1];			// forward path-metric[time index][state]
   for(i=0; i<=Packet_length; i++)
      mju_f[i] = new double[Num_state];
   mju_b = new double*[Packet_length+1];       	// backward path-metric[time index][state]
   for(i=0; i<=Packet_length; i++)
      mju_b[i] = new double[Num_state];
   branch_metric = new double**[Packet_length];	// branch[time index][state][input]
   for(i=0; i<Packet_length; i++)
      branch_metric[i] = new double*[Num_state];
   for(i=0; i<Packet_length; i++)
   	for(j=0; j<Num_state; j++)
         branch_metric[i][j] = new double[Num_in];

   // Initialize survival path
   for(i=0; i<Num_state; i++)
   {
// 	survival_path[i].metric = DBL_MAX;	// Initial maximum value for Euclidean distance
		survival_path[i].metric = -DBL_MAX;	// Initial minimum value for cross-correlation
//    mju_f[0][i] = DBL_MAX;					// Initial maximum value for Euclidean distance
//    mju_b[K][i] = DBL_MAX;					// Initial maximum value for Euclidean distance
		mju_f[0][i] = -DBL_MAX;					// Initial minimum value for cross-correlation
      mju_b[K][i] = -DBL_MAX;					// Initial minimum value for cross-correlation
	}
   survival_path[0].metric = 0.0;
   mju_f[0][0] = 0.0;
   mju_b[K][0] = 0.0;

/*********************/
/* Forward Recursion */
/*********************/
	for(i=0; i<Packet_length; i++)
   {
   	for(j=0; j<Num_state; j++)					// Initialize the survival path metric
      	survival_metric[j] = -DBL_MAX;

      for(j=0; j<Num_state; j++)					// from_state index
      	for(l=0; l<Num_in; l++)					// input bit
         {
         	// branch metric, Euclidean Distance
/*          branch_metric[i][j][l] = 0.0;
				branch_metric[i][j][l] += pow(Data_in[2*i]-(2*Trellis[j][l].out[0]-1),2);		// brahch metric
            branch_metric[i][j][l] += pow(Data_in[2*i+1]-(2*Trellis[j][l].out[1]-1),2);	// branch metric
            metric = survival_path[j].metric + branch_metric[i][j][l];
*/
				branch_metric[i][j][l] = 0.0;		// branch metric, Cross-correlation
            branch_metric[i][j][l] += (Data_in[2*i] * (2*Trellis[j][l].out[0]-1));		// brahch metric
            branch_metric[i][j][l] += (Data_in[2*i+1] * (2*Trellis[j][l].out[1]-1));	// branch metric
            metric = survival_path[j].metric + branch_metric[i][j][l];

            // find the survival path metric
//          if(metric < survival_metric[Trellis[j][l].to])	//	Euclidean distance (Minimize)
				if(metric > survival_metric[Trellis[j][l].to])	// Cross-correlation (Maximize)
            {
            	survival_metric[Trellis[j][l].to] = metric;

               // Record and refresh the survival path
               /*for(q=0; q<i; q++)
               {
               	survival_temp[Trellis[j][l].to].data_in[q] = survival_path[j].data_in[q];
                  survival_temp[Trellis[j][l].to].code_bit[2*q] = survival_path[j].code_bit[2*q];
                  survival_temp[Trellis[j][l].to].code_bit[2*q+1] = survival_path[j].code_bit[2*q+1];
                  survival_temp[Trellis[j][l].to].state[q] = survival_path[j].state[q];
               } */
               survival_temp[Trellis[j][l].to].data_in[i] = l;
               survival_temp[Trellis[j][l].to].code_bit[2*i] = Trellis[j][l].out[0];
               survival_temp[Trellis[j][l].to].code_bit[2*i+1] = Trellis[j][l].out[1];
               //survival_temp[Trellis[j][l].to].state[i] = Trellis[j][l].to;
               survival_temp[Trellis[j][l].to].state[i] = Trellis[j][l].from;
            }
			}

		// Record and refresh the survival path
      for(j=0; j<Num_state; j++)		// to_state index
      {
      	survival_path[j].metric = survival_metric[j];
         mju_f[i+1][j] = survival_metric[j];
         /*for(q=0; q<=i; q++)
         {
         	survival_path[j].data_in[q] = survival_temp[j].data_in[q];
            survival_path[j].code_bit[2*q] = survival_temp[j].code_bit[2*q];
            survival_path[j].code_bit[2*q+1] = survival_temp[j].code_bit[2*q+1];
            survival_path[j].state[q] = survival_temp[j].state[q];
         } */
      }
	}

   for(j=0; j<Num_state; j++)		// to_state index   
   {
      survival_path[j].data_in[Packet_length-1] = 
		          survival_temp[j].data_in[(Packet_length-1)];
      survival_path[j].code_bit[2*(Packet_length-1)] = 
		          survival_temp[j].code_bit[2*(Packet_length-1)];
      survival_path[j].code_bit[2*(Packet_length-1)+1] = 
		          survival_temp[j].code_bit[2*(Packet_length-1)+1];
  //    survival_path[j].state[Packet_length-1] = 
	//	          survival_temp[j].state[Packet_length-1];
   }

   for(j=0; j<Num_state; j++)		// to_state index
   {
   	pre_state = survival_temp[j].state[Packet_length-1];  // from state

   	for( q=Packet_length-2; q>=0; q--)
   	{
      	survival_path[j].data_in[q] = survival_temp[pre_state].data_in[q];
	      survival_path[j].code_bit[2*q] = survival_temp[pre_state].code_bit[2*q];
   	   survival_path[j].code_bit[2*q+1] = survival_temp[pre_state].code_bit[2*q+1];
      	// survival_path[j].state[q] = survival_temp[pre_state].state[q];
	      pre_state = survival_temp[pre_state].state[q];  // from state
   	}
  	}

/****************************************/
/* Backward Recursion and Soft Decision */
/****************************************/
	for(i=Packet_length-1; i>=0; i--)
   {
   	for(j=0; j<Num_state; j++)				// Initialize the survival path metric
//    	survival_metric[j] = DBL_MAX;		// Initial maximum value for Euclidean distance
			survival_metric[j] = -DBL_MAX;	// Initial minimum value for cross-correlation

		for(j=0; j<Num_state; j++)		// from_state index
      	for(l=0; l<Num_in; l++)		// input bit
         {
         	metric = mju_b[i+1][Trellis[j][l].to] + branch_metric[i][j][l];

            // find the survival path metric
//          if(metric < survival_metric[j])	//	Euclidean distance (Minimize)
				if(metric > survival_metric[j])	// Cross-correlation (Maximize)
            	survival_metric[j] = metric;
			}

		// Record the survival path metric
      for(j=0; j<Num_state; j++)		// from_state index
      	mju_b[i][j] = survival_metric[j];

      // LLR Calculation for the information bit
      mju[survival_path[0].data_in[i]] = mju_f[Packet_length][0];

//    mju[(survival_path[0].data_in[i]+1)%2] = DBL_MAX;		//	Euclidean distance (Minimize)
		mju[(survival_path[0].data_in[i]+1)%2] = -DBL_MAX;    // Cross-correlation (Maximize)
      for(j=0; j<Num_state; j++)		// from_state index
      {
      	mju_tmp = mju_f[i][j] + branch_metric[i][j][(survival_path[0].data_in[i]+1)%2]
         			 + mju_b[i+1][Trellis[j][(survival_path[0].data_in[i]+1)%2].to];

//       if(mju_tmp < mju[(survival_path[0].data_in[i]+1)%2])	//	Euclidean distance (Minimize)
			if(mju_tmp > mju[(survival_path[0].data_in[i]+1)%2])	// Cross-correlation (Maximize)
         	mju[(survival_path[0].data_in[i]+1)%2] = mju_tmp;
		}

//    Soft_out[i] = mju[0] - mju[1];		// Euclidean Distance
		Soft_out[i] = mju[1] - mju[0];		// Cross-correlation
/*
      // LLR Calculation (for the 1st code bit)
      mju[survival_path[0].code_bit[2*i]] = mju_f[Packet_length][0];

//    mju[(survival_path[0].code_bit[2*i]+1)%2] = DBL_MAX;			//	Euclidean distance (Minimize)
		mju[(survival_path[0].code_bit[2*i]+1)%2] = -DBL_MAX;    	// Cross-correlation (Maximize)
      for(j=0; j<Num_state; j++)		// from_state index
      	for(l=0; l<Num_in; l++)
         	if(Trellis[j][l].out[0] != survival_path[0].code_bit[2*i])
            {
            	mju_tmp = mju_f[i][j] + branch_metric[i][j][l] + mju_b[i+1][Trellis[j][l].to];
//             if(mju_tmp < mju[(survival_path[0].code_bit[2*i]+1)%2])	//	Euclidean distance (Minimize)
					if(mju_tmp > mju[(survival_path[0].code_bit[2*i]+1)%2])	// Cross-correlation (Maximize)
               	mju[(survival_path[0].code_bit[2*i]+1)%2] = mju_tmp;
				}

//    Soft_out[2*i] = mju[0] - mju[1];		// Euclidean Distance
		Soft_out[2*i] = mju[1] - mju[0];		// Cross-correlation

		// LLR Calculation (for 2nd code bit)
      mju[survival_path[0].code_bit[2*i+1]] = mju_f[Packet_length][0];

//    mju[(survival_path[0].code_bit[2*i+1]+1)%2] = DBL_MAX;		//	Euclidean distance (Minimize)
		mju[(survival_path[0].code_bit[2*i+1]+1)%2] = -DBL_MAX;		// Cross-correlation (Maximize)
      for(j=0; j<Num_state; j++)		// from_state index
      	for(l=0; l<Num_in; l++)
         	if(Trellis[j][l].out[1] != survival_path[0].code_bit[2*i+1])
            {
            	mju_tmp = mju_f[i][j] + branch_metric[i][j][l] + mju_b[i+1][Trellis[j][l].to];
//   	         if(mju_tmp < mju[(survival_path[0].code_bit[2*i+1]+1)%2])	//	Euclidean distance (Minimize)
					if(mju_tmp > mju[(survival_path[0].code_bit[2*i+1]+1)%2])	// Cross-correlation (Maximize)
               	mju[(survival_path[0].code_bit[2*i+1]+1)%2] = mju_tmp;
				}

//    Soft_out[2*i+1] = mju[0] - mju[1];		// Euclidean Distance
		Soft_out[2*i+1] = mju[1] - mju[0];		// Cross-correlation
*/
      Data_out[i] = survival_path[0].data_in[i];
   }

	delete survival_metric;
   delete mju;
   for(i=0; i<=Packet_length; i++)
      delete mju_f[i];
   delete mju_f;
   for(i=0; i<=Packet_length; i++)
       delete mju_b[i];
   delete mju_b;
   for(i=0; i<Packet_length; i++)
   	for(j=0; j<Num_state; j++)
         delete branch_metric[i][j];
   for(i=0; i<Packet_length; i++)
      delete branch_metric[i];
   delete branch_metric;
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

