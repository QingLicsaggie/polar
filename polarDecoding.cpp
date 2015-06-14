#include<cstdio>
#include<iostream>
#include<fstream>
#include<set>
#include<math.h>
#include<cstdlib>
#include"generator.h"
#include"polarDecoding.h"
#include"polarEncoding.h"
#include"channelPolarization.h"
using namespace std;
extern int N;
extern int K;

/*--------------------------------------
FUNCTION:
	void decodingSC(int *y, int *decodedcodeword, set <int>frozenIndex, int *frozenBitArray, double **W)
DESCRIPTION:
	This function implements the Successive Decoding of polar codes
PARAMETERS:
	INPUT:
		y--received message
		decodedcodeword--storing the decoded codeword
		frozenIndex--
		frozenBitArray--
		W--channel information
	OUTPUT:
		HIdecodedcodeword
RETURN VALUES:
	None
-------------------------------------------*/
void decodingSC(int *y, int *decodedcodeword, set<int> frozenIndex, int *frozenBitArray, double **W){
	
	int i, FIndex = 0;
	double hi;
	for(i = 0; i < N; i++){	
		//frozen bit case
		if(frozenIndex.find(i) != frozenIndex.end()){
			//cout<<"frozen bit case ";
			decodedcodeword[i] = frozenBitArray[FIndex++];
		}
		//information set case
		else{ 
 			hi = recursiveChannelTransformations(N, i, y,  0, N - 1, decodedcodeword, 0, W)/recursiveChannelTransformations(N, i, y, 0, N - 1, decodedcodeword, 1, W);
			decodedcodeword[i] = (hi > 1.0)? 0: 1;							
		}
	}
}		
