/**********************************************/
/* Author: Chao-wang Huang                  */
/* Date: Thursday, September 23, 2004 */
/* S-random interleaver generator          */
/**********************************************/

#include "stdafx.h"
#include <stdio.h>
#include <math.h>
//#include <time.h>
#include <ctime>
#include <conio.h>
#include <stdlib.h>
#include <iostream>
using namespace std;

#define N 8192		// Interleaver size
#define S 16		// S parameter of S-random interleaver

int main(void)
{
	time_t t, start, end;
	int i, j, k, sol, position, flag, count, *interleaver, *check;
	FILE *fp1;

	start = time(NULL);
	printf("S-random Interleaver Generator\n");
	printf("This program is running. Don't close, please!\n\n");
	srand((unsigned) time(&t));

	interleaver = new int[N];
	check = new int[N];
	fopen_s(&fp1, "s_random.log", "w");

	count = 0;
	do
	{
		count ++;
  		for(i=0; i<N; i++)
   		check[i] = 0;

  		for(i=0; i<N; i++)
   		{
   			do
      		{
      			flag = 0;
      			//position = random(N);
				position = rand() % N;

         		if(check[position]==1)
         			flag = 1;
         		else
         		{
         			for(j=i-1; (j>i-S && j>=0); j--)
         				if(abs(position - interleaver[j]) < S)
               			{
               				flag = 1;
                  			break;
               			}

            		if(flag==1)
         			{
            			for(k=0; k<N; k++)
            				if(check[k]==0)
                  			{
                  				sol = 1;
               					for(j=i-1; (j>i-S && j>=0); j--)
                     				if(abs(k - interleaver[j]) < S)
                     				{
                     					sol = 0;
                        				break;
                     				}

                     			if(sol==1)
                        			break;
                  			}

            			if(sol==0)
            			{
            				flag = 0;
            			}
         			}
         		}
      		}
      		while(flag);

//         printf("%d %d\n", i, position);
//         cout << i << " " << position << endl;
      		check[position] = 1;
      		interleaver[i] = position;
   		}

		if(count %10000 == 0)
			cout << "Iteration: " << count << endl;
	}
	while(sol==0);

	cout << endl << "Iteration: " << count << endl;
	for(i=0; i<N; i++)
   		fprintf(fp1, "%d %d\n", i, interleaver[i]);

	delete interleaver;
	delete check;
  	fclose(fp1);

	end = time(NULL);
	printf("Total elapsed time: %.0f(sec)\n", difftime(end,start));
	printf("This program is ended. Press any key to continue.\n");
	getchar();

	return 0;
}
