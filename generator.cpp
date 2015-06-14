#include<stdio.h>
#include<math.h>
#include<iostream>
using namespace std;
extern int trial;
extern int N;
extern int K;
extern int q;
/*---------------------------------------------------------------
FUNCTION: 
	void int2bin(int intstat, int *tempstat, int length)
	
DESCRIPTION:
	This function transfers a decimal integer to a binary sequence.

PARAMETERS:
	INPUT:
		intstat - Decimal integer needed to be changed.
		length - Length of binary sequence.
	OUTPUT:
		tempstat - Contains pointer to the binary sequence.

RETURN VALUE:
	None
---------------------------------------------------------------*/
void int2bin(int intstat, int *tempstat, int length)
{
	int i, temp;
	temp = intstat;
	for (i=length-1; i>=0; i--)
	{
		tempstat[i] = temp%2;
		temp = temp/2;
	}
}
/*-------------------------------------------
FUNCTION:
	int bin2int(int *binseq, int length)
DESCRIPTION:
	This function transfers a binary sequence to a decimal integer.
PARAMETERS:
	INPUT:
		binseq-Contains pointer to the binary sequence.
		length-Length of binary sequence.
RETURN VALUE:
	The decimal integer corresponding to the binary sequence.
--------------------------------------------*/
long double bin2int(int *binseq, int length)
{
	int i, j, temp;
	long double sum = 0;

	for(i = 0; i < length; i++){	
		temp = 1;
		for(j = 1; j <=i; j++){
			temp = temp*2;
		}
		sum = sum + temp*(*(binseq + length - 1 - i));
	}
	return sum;
}
/*---------------------------------------------
FUNCTION:
	void generator(int **generator, int n)
DESCRPTION:
	This function generators the polar generator according the Proposition 16 of the paper 
	"Channel Polarization: A method for constructing capacity-achieving codes for symmetric
	binary-input memoryless channel" by Erdal Arikan.
PARAMETERS:
	INPUT:
		generator-contains the generator matrix
		n-let N = 2**n, and the size of the generator matrix is N by N.
	OUTPUT:
		generator
RETURN VALUE:
	The generator for Polar code
----------------------------------------------*/
void generator(int **generator, int n)
{	
	int i, j, k;
	int indexForX[n], indexForY[n];
	N = (int)pow(2.0, n);
	for(i = 0; i < N; i++){
		int2bin(i, indexForX, n);
		for(j = 0; j < N; j++){
			generator[i][j] = 1;
			int2bin(j, indexForY, n);
			for(k = 0; k < n; k++){
				generator[i][j] *= (1 + indexForY[k] + indexForX[n-1-k]*indexForY[k])%2;						
			}			

		}
	}
	/*cout<<"The generator for polar is "<<endl;
	for(i = 0; i < N; i++){
	  for(j = 0; j < N; j++){
	    cout<<generator[i][j]<<" ";
	  }
	cout<<endl;
	}
	cout<<endl;*/
}
/*----------------------------------------
FUNCTION:
       void polarFMatrix(int **generator, int n)
DESCRIPTION:
       Generator F matrix for polar code based on equation 71 of channel polarizaiton
PARAMETERS:
       INPUT:
            generator
	    n        -- The row number of the matrix is 2^n
       OUTPUT:
            generator
RETURN VALUES:
       NONE
LAST UPDATE DATE
     06/12/2014
  ----------------------------------------*/
void polarFMatrix(int **generator, int n)
{	
	int i, j, k;
	int indexForX[n], indexForY[n];
	N = (int)pow(2.0, n);
	for(i = 0; i < N; i++){
		int2bin(i, indexForX, n);
		for(j = 0; j < N; j++){
			generator[i][j] = 1;
			int2bin(j, indexForY, n);
			for(k = 0; k < n; k++){
				generator[i][j] *= (1 + indexForY[k] + indexForX[k]*indexForY[k])%2;						
			}			

		}
	}
	/*cout<<"The generator for polar is "<<endl;
	for(i = 0; i < N; i++){
	  for(j = 0; j < N; j++){
	    cout<<generator[i][j]<<" ";
	  }
	cout<<endl;
	}
	cout<<endl;*/
}

/*----------------------------------
obtain the j-th column of polar generator matrix
--------------------------------------*/
void generatorC(int *generator, int j){
	int i, k;
	int n = (int)log2((double)N);
	int indexForX[n];
	int indexForY[n];	
	int2bin(j, indexForY, n);
	for(i = 0; i < N; i++){
		int2bin(i, indexForX, n);
		generator[i] = 1;
		for(k = 0; k < n; k++){
			generator[i] *= (1 + indexForY[k] + indexForX[n-1-k]*indexForY[k])%2;						
		}		
	}
}
/*-----------------------------------
obtain the i-th row of polar generator matrix
  -------------------------------------*/
void generatorR(int *generator, int i){
	int  j, k;
	int n = (int)log2((double)N);
	int indexForX[n], indexForY[n];

	int2bin(i, indexForX, n);
	for(j = 0; j < N; j++){
		generator[j] = 1;
		int2bin(j, indexForY, n);
		for(k = 0; k < n; k++){
			generator[j] *= (1 + indexForY[k] + indexForX[n-1-k]*indexForY[k])%2;						
		}			
	}
}
