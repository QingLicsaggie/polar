#include<cstdio>
#include<set>
#include<iterator>
#include<iostream>
#include"polarEncoding.h"
#include"generator.h"
#include"channelPolarization.h"
using namespace std;
extern int K;
extern int N;
/*-------------------------------------
FUNCTION:
	void setInput(int *information, int *frozenIndex, int frozenBit, int *Message)
DESCRIPTION:
	This function sets input by filling the frozen bit with fronzenBit 
	non-frozen bit with random numbers
PARAMETERS:
	INPUT:
		information--stores the message
		frozenIndex--frozenbit index
		frozenBit--binary array
		Message
	OUTPUT:
		Message
RETURN VALUES:
	None 
--------------------------------------*/
void setInput(int *Information, set<int>frozenIndex, int *frozenBit, int *Message){	
	int i, Mindex = 0, Findex = 0;	
	for(i = 0; i < N; i++){
		if (frozenIndex.find(i) != frozenIndex.end()){
			Message[i] = frozenBit[Findex++];
		}
		else		
			Message[i] = Information[Mindex++];
	}
}

/*---------------------------------------
FUNCTION:
	void polarEncode(int *codeword, int *information, int *frorzenBit,  int *frozenIndex)

DESCRIPTION:
	This function implements the Polar Encoding process.
PARAMETERS:
	INPUT:
		codeword--
		information--original message(binary form)
		frozenBit--frozenBit(binary form)
		frozenIndex--indicates which indexes are frozen
	OUTPUT:
		codeword
RETURN VALUES:
	None
----------------------------------------*/
void polarEncode(int *codeword, int *information, int *frozenBit, set<int>frozenIndex){
	
	int messageBinary[N];
	setInput(information, frozenIndex, frozenBit, messageBinary);
	encode(codeword, messageBinary);
}
