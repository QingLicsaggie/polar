#include<stdio.h>
#include<math.h>
#include<stdlib.h>
#include<string>
#include<vector>
#include<set>
#include<iterator>
#include<time.h>
#include<sys/time.h>
#include<sys/resource.h>
#include<iostream>
#include<cstring>
#include<sstream>
#include<fstream>
#include<climits>
#include"analysis.h"
#include"generator.h"
#include"channelPolarization.h"
#include"polarEncoding.h"
#include"noise.h"
#include"polarDecoding.h"
#include"listdecoding.h"
#include"SC.h"
#include"rewritting.h"
#include"noisywem.h"
#include"crc16.h"
using namespace std;
int Dimension;
int N;
int q = 2;
int K;
int OUTPUTSIZE;
int trial = 1;
long double prob;
double average_error_rate = 0.0;
double theoretical_upper_bound = 1.0;
double flipCount1 = 0.0;
double flipCount2 = 0.0;
double distribution = 0.0;
/*-------------------------------------
FUNCTION: 
	void readMatrixes(double **W)
DESCRIPTION:
	This function obtains the channel information from file.
PARAMETERS:
	INPUT:
		W-channel transition probability
	OUTPUT
		The channel transition probability.
RETURN VALUES:
	None
--------------------------------------*/
void readMatrixes(double **W, double **M, long double prob){
	
	W[0][0] = W[1][1] = 1 - prob;
        W[0][1] = W[1][0] = prob;

	//cost is measured by hamming distance
	M[0][0] = M[1][1] = 0;
         M[0][1] = 1;  M[1][0] = 1;
	
}
/*-----------------------------------
FUNCTION:
	void channelCoding_experiment_simulation(set<int>frozenIndex, double **W, int **G)
DESCRIPTION:
	INPUT:
		frozenIndex--
		W--
		G--Generator
	OUTPUT:
		None
RETURN VALUES:
	None
-------------------------------------*/
void channelCoding_experiment_simulation(set<int>frozenIndex, double **W, const vector<long double>&channelInfo){
  average_error_rate = 0;
  int codeword[N], information[K], frozenBit[N - K], corruptedcodeword[N], decodedresult[N], Message[N];
  int flip = 0;
  for(int i = 0; i < trial; i++){
    //crc 16
    for(int i = 0; i < K-16; i++)
      information[i] =  rand()%2;
   
    int crc = crc16(information, K-16);
    int crcBit[16];
    int2bin(crc, crcBit, 16);
    for(int i = 0; i < 16; i++)
      information[K-16+i] = crcBit[i];

    for(int i = 0; i < N - K; i++)
      frozenBit[i] = 0;
	   	
    setInput(information, frozenIndex, frozenBit, Message);
    polarEncode(codeword, information, frozenBit, frozenIndex);
    BSC_corrupted(codeword, corruptedcodeword, W);
   

    flip += countFlip(codeword, corruptedcodeword);
    
    listDecoder ld(2);
    ld.SCLDecoder(corruptedcodeword, decodedresult, frozenIndex, frozenBit, W);

   

    average_error_rate += (double)blockErrorRate(Message, decodedresult);
    if(i%100==0){
      cout<<"Now, trial is "<<i<<endl;
      cout<<"so far average block error rate is "<<average_error_rate/(i+1)<<endl;
      cout<<"on average flip "<<(double)flip/(N*(1+i))<<endl;
    }
  }
   average_error_rate = average_error_rate/trial;
   cout<<"average error rate is "<<average_error_rate<<endl;
}
/*---------------------------------------------------------*/
int main(int argc, char*argv[]){
	/*------Get parameters from input--------*/
        OUTPUTSIZE = atoi(argv[1]);
	Dimension = atoi(argv[2]);
	K = atoi(argv[3]);
	long double prob = atof(argv[4]);
	char *fileP = argv[5];
	N = (int)pow(2.0, Dimension);
	double **W  = (double **)malloc(q*sizeof(double*));
	double **M = (double **)malloc(q*sizeof(double*));
  	for(int i = 0 ; i < q; i++){
	  W[i]  = (double*)malloc(OUTPUTSIZE*sizeof(double));
		M[i] = (double*)malloc(OUTPUTSIZE*sizeof(double));
	}		
	srand(time(NULL));
	readMatrixes(W, M, prob);
	

        /*--------Construct polar codes---------*/
	//DFSCalPforAllBitChannels(256, W);
	vector<long double> P;
	readVP(P, fileP);
	set<int> frozenIndex;
	
	
	/*-------Source coding ----------*/
	/*K = (compute_I(W))*N;
	cout<<K<<endl;
        selectChannels(P, frozenIndex, K); 
	sourceCodingVerify(W, M, frozenIndex);*/
	
	/*--------Polar WEM-----------------------*/
        /*K =  (compute_I(W))*N;
	selectChannels(P, frozenIndex, K);
	experiment_for_rewriting(W, M, frozenIndex);*/

	
	/*-------Noise WEM--------------------------*/
	/*int KForP = K;
	int KForW = (compute_I(WForW))*N;
	while(theoretical_upper_bound >= 0.00005 && KForP > 0){
		theoretical_upper_bound = 0;
		KForP--;
	        //K--;
	 	selectChannels(P, frozenIndex, KForP);
	
	//selectChannels(P, frozenIndex, K);
	}
	noisyWEM nw(KForP, KForW, L, prob, probforWEM);
        nw.experiment_of_noisyWEM(fileP, fileW, trial, W, WForW, M);*/
	
	
	/*--------channel coding part ------------*/
	channelCoding_experiment_simulation(frozenIndex, W, P);
	//print_info(frozenIndex, db, prob);
	
       /*---------clean up-------------------------*/
	for(int i = 0 ; i < q; i++){
		free (W[i]); 
		free (M[i]);
	}
	free(W);
	free(M);
	return 1;
}
