/*=====================================================================*/
/* Author: Chao-wang Huang                                             */
/* Date: Monday, June 19, 2006                                         */
/* An (20,5) Low-Density Parity-Check code (LDPC) is simulated         */
/* Decoding algorithm: Message Passing Algorithm (Belief Propagation)  */
/*=====================================================================*/

#include <stdio.h>
#include <math.h>
#include <float.h>
#include <limits.h>
#include <stdlib.h>
#include <time.h>
#include <conio.h>
#include <iostream.h>

void AWGN_noise(float, double, double *);
int Error_count(int, int);
//void Belief_Propagation(double *, double ***, double, int);
void Belief_Propagation(double *, double, int);

//#define Pi 				3.14159265358979
#define MAX_RANDOM   LONG_MAX   			// Maximum value of random()
#define MAX_NODES    20   					// Maximum of number of code/check nodes
#define MAX_CHK      3   					// Maximum number of checks per code bit
#define MAX_BIT      4   					// Maximum number of code bits per check
#define num_packet 	500000				// number of packets simulated
int max_weight_M;                 		// Maximum column weight
int max_weight_N;                 		// Maximum row weight
int num_check_M[MAX_NODES];           	// Number of check nodes of a bit node
int num_check_N[MAX_NODES];           	// Number of bit nodes of a check node
int check_set_M[MAX_NODES][MAX_CHK];	// Check nodes set of a bit node
int check_set_N[MAX_NODES][MAX_BIT];	// Bit nodes set of a check node
const int Iteration = 10;	   			// Iterations of Turbo Equalization ( >= 1 )
const int n = 20;               			// Code word length of a linear block code
const int k = 5;                			// Information bits of a linear block code
int n_k;                          		// Redundancy (Parity check bits)
float rate;                      		// code rate
int M, N;                          		// Size of parity-check matrix (Row,Column) = (M,N)
double q[20][10][2], snr_rms;

int main(void)
{
	time_t t, start, end;
	int i, j, p, iter, *data_bit, /*coded_bit,*/ **Dk;
   int *err_count;
   double *bpsk, snr, Eb_No, noise_pwr, noise[2], *Yk;
   double *err_rate/*, ***out_prob*/;
   FILE *parity_check, *ber, *records;

   start = time(NULL);
   printf("BER Performance of Galager (20,5) LDPC Code in AWGN Channel\n");
   printf("Dimension of Parity Check matrix: (20,15)\n");
   printf("Column Weight = 3, Row Weight = 4\n");
   printf("This program is running. Don't close, please!\n\n");

   records = fopen("Records_LDPC_BP_awgn.log", "a");
   fprintf(records, "BER Performance of Galager (20,5) LDPC Code in AWGN Channel\n");
   fprintf(records, "Dimension of Parity Check matrix: (20,15)\n");
   fprintf(records, "Column Weight = 3, Row Weight = 4\n");
   fprintf(records, "number of bits of simulation = %d\n\n", k*num_packet);
   fprintf(records, "Eb/No     BER\n");
   fflush(records);

   data_bit = new int[k];
   //coded_bit = new int[n];
   bpsk = new double[n];
	Dk = new int*[n];
   for(i=0; i<n; i++)
   	Dk[i] = new int[Iteration];
   Yk = new double[n];								// Received signal
/*   out_prob = new double**[MAX_NODES];			// estimated bit probability
   for(i=0; i<MAX_NODES; i++)          		// out_prob[bit node][iteration][bit]
   	out_prob[i] = new double*[Iteration];
   for(i=0; i<MAX_NODES; i++)
   	for(j=0; j<Iteration; j++)
   		out_prob[i][j] = new double[2];*/
   err_count = new int[Iteration];
   err_rate = new double[Iteration];

	srand((unsigned) time(&t));
	n_k = n-k;
  	rate = (float) k / (float) n;

	// Define Parity Check Matrix
  	parity_check = fopen("gallager_20_4_3.log", "r");
   fscanf(parity_check, "%d %d", &N, &M);
   fscanf(parity_check, "%d %d", &max_weight_N, &max_weight_M);
   for (i=0; i<M; i++)
   	fscanf(parity_check, "%d", &num_check_N[i]);

   for (i=0; i<N; i++)
   	fscanf(parity_check, "%d", &num_check_M[i]);

   for (i=0; i<M; i++)
   	for (j=0; j<num_check_N[i]; j++)
      	fscanf(parity_check, "%d", &check_set_N[i][j]);

   for (i=0; i<N; i++)
   	for (j=0; j<num_check_M[i]; j++)
      	fscanf(parity_check, "%d", &check_set_M[i][j]);

   fclose(parity_check);

/*===========================================*/
/* M A I N    S I M U L A T I O N    L O O P */
/*===========================================*/
   ber=fopen("ber_LDPC_20_5.log", "w");
   for(snr=0; snr<=10; snr++)
   {
   	for(iter=0; iter<Iteration; iter++)
   		err_count[iter] = 0;
   	// noise power calculation
      Eb_No = (double)snr;
      noise_pwr = 1.0/(((float)k/(float)n)*pow(10.0, Eb_No/10.0));	// BPSK, Nyquist filter assumption
      snr_rms = 2.0 * sqrt(2.0*(float)k/(float)n*(pow(10.0,(snr/10.0))));
      printf("Eb_No = %f, ", Eb_No);

      for(p=0; p<num_packet; p++)
      {
      	for (i=0; i<k; i++)
          	data_bit[i] = 0;

   		/*for(i=0; i<k; i++)
				data_bit[i] = random(2);		// Generate random information bit stream
           */
         /* LDPC code Encoding */

			for(i=0; i<n; i++)		// BPSK mapping
   			//if(coded_bit[i]==1)
      			bpsk[i] = 1.0;
      		//else
      		  //	bpsk[i] = -1.0;

         /* AWGN channel */
         for(i=0; i<n; i++)
         {
         	AWGN_noise(0, noise_pwr, &noise[0]);
         	Yk[i] = bpsk[i] + noise[0];
            // Soft Decision
//            Dk[i] = tanh(Yk[i]/noise_pwr);
//            Dk[i] = Yk[i];
         }

/* Message Passing Decoding (Belief Propagation) */
			//Belief_Propagation(Yk, out_prob, noise_pwr, Iteration);
         Belief_Propagation(Yk, noise_pwr, Iteration);

         // Hard Decision
         for(i=0; i<N; i++)
         	for(iter=0; iter<Iteration; iter++)
            {
	         	//if(out_prob[i][iter][0] >= 0.5)
               if(q[i][iter][0] >= 0.5)
            		Dk[i][iter] = 0;
            	else
            		Dk[i][iter] = 1;

         		// Bit error count
               if(Dk[i][iter] == 1)
         			err_count[iter]++;
            }
		}

      // Statistics and records
      cout << "Error Rate = ";
      for(iter=0; iter<Iteration; iter++)
      {
      	err_rate[iter] = err_count[iter] / (double)(N*num_packet);
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
   //delete coded_bit;
   delete bpsk;
   delete Yk;
   for(i=0; i<N; i++)
   	delete Dk[i];
   delete Dk;
/*   for(i=0; i<MAX_NODES; i++)
      for(j=0; j<Iteration; j++)
   		delete out_prob[i][j];
   for(i=0; i<MAX_NODES; i++)
      delete out_prob[i];
   delete out_prob;*/
   delete err_count;
   delete err_rate;

   end = time(NULL);
   printf("Total elapsed time: %.0f(sec)\n", difftime(end,start));
   fprintf(records, "\nTotal elapsed time: %.0f(sec)\n\n", difftime(end,start));
	fclose(ber);
   fclose(records);
   printf("This program is ended. Press any key to continue.\n");
   getch();

   return 0;
}

void Belief_Propagation(double *received/*, double ***q/* output probability */, double noise_pwr, int Iteration)
{
// Iterative decoding by belief propagation in code's Bayesian network
// Based on Pearl's book and MacKay's paper

	int i, j, l, iter;
   double alpha, **p1, ***qm, ***rm, **delta;


   p1 = new double*[N];						// Probability of "1" and "0" were sent
	for(i=0; i<N; i++)						// p1[bit node][bit]
   	p1[i] = new double[2];
   qm = new double**[M];					// Probability from bit node to check node
	for(i=0; i<M; i++)						// qm[check node][bit node][bit]
      qm[i] = new double*[N];
   for(i=0; i<M; i++)
      for(j=0; j<N; j++)
   		qm[i][j] = new double[2];
   rm = new double**[M];					// Probability from check node to bit node
	for(i=0; i<M; i++)          			// rm[check node][bit node][bit]
      rm[i] = new double*[N];
   for(i=0; i<M; i++)
      for(j=0; j<N; j++)
   		rm[i][j] = new double[2];
   delta = new double*[M];					// qm[i][j][0] - qm[i][j][1]
   for(i=0; i<M; i++)          			// delta[check node][bit node]
   	delta[i] = new double[N];

  	// -------------------
  	// INITIALIZATION STEP
  	// -------------------
  	// a-priori probabilities
	for (i=0; i<N; i++)		// for every bit node
   {
    	// p1[i] = 1.0 / ( 1.0 + exp(fabs(received[i])*snr_rms) );
    	// p1[i] = 1.0 - exp(-snr_rms*fabs(received[i])) /
    	//         (1.0 + exp(-snr_rms*fabs(received[i])));
    	//p1[i] = 1.0 / ( 1.0 + exp(received[i]*snr_rms) );
      //p1[i][1] = 1.0 / ( 1.0 + exp(received[i]*2.0/noise_pwr) );
      				// probability of a "-1" was transmitted
      p1[i][1] = 1.0 / ( 1.0 + exp(received[i]*snr_rms) );
      p1[i][0] = 1.0 - p1[i][1];
   }

  	// For every (m,l) such that there is a link between parents and
  	// children, qm0[i][j] and qm1[i][j] are initialized to pl[j].
  	// Notation: pi (Pearl) = q (MacKay)

  	for (i=0; i<M; i++)		// for every check node
    	for (j=0; j<num_check_N[i]; j++)		// �o�@�B�J�n���S�N�q???
      {
      	qm[i][j][0] = 0.0;
      	qm[i][j][1] = 0.0;
      }

  	for (i=0; i<N; i++)		// for every bit node
    	for (j=0; j<num_check_M[i]; j++)	// for check nodes of every bit node
      {
      	qm[check_set_M[i][j]-1][i][0] = p1[i][0];
      	qm[check_set_M[i][j]-1][i][1] = p1[i][1];
      }

   for(iter=0; iter<Iteration; iter++)		// iterations of massage passing
   {
  		// ---------------------------------------
  		// HORIZONTAL STEP = BOTTOM-UP PROPAGATION
	  	// ---------------------------------------
  		for (i=0; i<M; i++)		// for every check node
    		for (j=0; j<num_check_N[i]; j++)		// for every bit of a check node
      	{
      		delta[i][check_set_N[i][j]-1] = 1.0;	// initial value for �֭�

      		for (l=0; l<num_check_N[i]; l++)
        			if (check_set_N[i][l] != check_set_N[i][j])
        				delta[i][check_set_N[i][j]-1] *= ( qm[i][check_set_N[i][l]-1][0]
                  	                                - qm[i][check_set_N[i][l]-1][1] );

      		rm[i][check_set_N[i][j]-1][0] = 0.5 * ( 1.0 + delta[i][check_set_N[i][j]-1] );
      		rm[i][check_set_N[i][j]-1][1] = 0.5 * ( 1.0 - delta[i][check_set_N[i][j]-1] );

      		// prevent from under-flow
      		if (rm[i][check_set_N[i][j]-1][0] == 0.0)
         	{
        			rm[i][check_set_N[i][j]-1][0] = 1.0e-10;
        			rm[i][check_set_N[i][j]-1][1] = 1.0 - 1.0e-10;
        		}
      	}

  		// ------------------------------------
  		// VERTICAL STEP = TOP-DOWN PROPAGATION
  		// ------------------------------------
  		for (i=0; i<N; i++)		// for every bit node
      	for (j=0; j<num_check_M[i]; j++)
      	{
		   	qm[check_set_M[i][j]-1][i][0] = 1.0;	// initial value for �֭�
      		qm[check_set_M[i][j]-1][i][1] = 1.0;

      		for (l=0; l<num_check_M[i]; l++)
        			if (check_set_M[i][l] != check_set_M[i][j])
          		{
          			qm[check_set_M[i][j]-1][i][0] *= rm[check_set_M[i][l]-1][i][0];
          			qm[check_set_M[i][j]-1][i][1] *= rm[check_set_M[i][l]-1][i][1];
	         	}

      		qm[check_set_M[i][j]-1][i][0] *= p1[i][0];
      		qm[check_set_M[i][j]-1][i][1] *= p1[i][1];

      		alpha = 1.0 / (qm[check_set_M[i][j]-1][i][0] + qm[check_set_M[i][j]-1][i][1]);

      		qm[check_set_M[i][j]-1][i][0] *= alpha;
      		qm[check_set_M[i][j]-1][i][1] *= alpha;

            // prevent from numerical under-flow
      		if (qm[check_set_M[i][j]-1][i][0] == 0.0)
        		{
            	qm[check_set_M[i][j]-1][i][0] = 1.0e-10;
               qm[check_set_M[i][j]-1][i][1] = 1.0 - 1.0e-10;
            }
      	}

		// MacKay: At this step we also compute the (unconditional) pseudo-
  		// posterior probalilities "q0, q1" to make tentative decisions

  		for (i=0; i<N; i++)		// for every bit node
    	{
      	q[i][iter][0] = 1.0;
    		q[i][iter][1] = 1.0;

    		for (j=0; j<num_check_M[i]; j++) // for every check node of the bit
      	{
      		q[i][iter][0] *= rm[check_set_M[i][j]-1][i][0];
      		q[i][iter][1] *= rm[check_set_M[i][j]-1][i][1];
      	}

    		q[i][iter][0] *= p1[i][0];
    		q[i][iter][1] *= p1[i][1];

    		alpha = 1.0 / (q[i][iter][0] + q[i][iter][1]);

    		q[i][iter][0] *= alpha;
    		q[i][iter][1] *= alpha;
    	}
   }

   for(i=0; i<N; i++)
   	delete p1[i];
   delete p1;
   for(i=0; i<M; i++)
      for(j=0; j<N; j++)
   		delete qm[i][j];
   for(i=0; i<M; i++)
      delete qm[i];
   delete qm;
   for(i=0; i<M; i++)
      for(j=0; j<N; j++)
   		delete rm[i][j];
   for(i=0; i<M; i++)
      delete rm[i];
   delete rm;
   for(i=0; i<M; i++)
   	delete delta[i];
   delete delta;
}





//void encode()
//
// Systematic encoding
//
//{
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
//}

void AWGN_noise(float mu, double variance, double *noise)
{
	const float Pi = 3.14159265358979;
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

