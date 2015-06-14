#include<cstdio>
#include<cstdlib>
#include<iostream>
#include<fstream>
#include<cstring>
#include<sstream>
#include<vector>
#include<time.h>
#include<math.h>
#include<set>
#include"analysis.h"
#include"channelPolarization.h"
#include"SC.h"
#include"generator.h"
#include"polarEncoding.h"
#include"polarDecoding.h"
#include"listdecoding.h"
using namespace std;
extern int trial;
extern int N;
extern int K;
extern int q;
extern long double prob;
int compareBit(int c1, int c2){
     int count = 0;
     int temp1[(int)log2(q)];
     int temp2[(int)log2(q)];
     int2bin(c1, temp1, (int)log2(q));
     int2bin(c2, temp2, (int)log2(q));
     for(int i = 0; i < (int)log2(q); i++)
	if (temp1[i] == temp2[i])
		count++;
    return count; 	
}
/*---------------------------
FUNCTION:
	double compute_D(double **W, double **M)
DESCRIPTION:
	This function computes d based on sokada 
PARAMETERS:
	INPUT:
		W--channel information
	        M--distoration
	OUTPUT:
		d
RETURN VALUES:
	d
-----------------------------------------------*/
double compute_D(double **W, double **M){
	
	double D = 0;
	/*for(int i = 0; i < q; i++){
		for(int j = 0; j < q; j++){
		int number = compareBit(i, j);
		W[i][j] = pow(prob, number)*pow(1 - prob, (int)log2(q) - number);		
		}
	}*/
	for(int y = 0; y < q; y++)
		for(int x = 0; x < q; x++){
			D += 0.5*W[x][y]*M[x][y];
                  
	}			
	return D;
}
/*--------------------------------
FUNCTION:
	long double measureDistoration(int *message, int *decodedcodeword, double **M)
DESCRIPTION:
	This function measures the distoration among message and decodedcodeword.
PARAMETERS:
	INPUT:
		message--
		decodedcodeword--
		M--Measure matrix
	OUTPUT:
		distoration
RETURN VALUES:
	distoration
---------------------------------*/
long double measureDistortion(int *message, int *decodedcodeword, double **M){
	
	int i;
	int distortion = 0;

	for(i = 0; i < N; i++){	
		distortion += M[message[i]][decodedcodeword[i]];
	}
	
	return (double)distortion/N; //((double)distortion/N < 0.5?  (double)distortion/N:  1- (double)distortion/N);
}
/*--------------------------------------------
Two small functions
----------------------------------------------*/
int measureBit(long double number){
	long double temp = number;
	int bit = 0;
	while(temp<1){
		bit++;
		temp = temp*10;
	}
	return bit;
}
int bin2dec(int *binseq, int length){
	long double sum = 0;
	for(int i = 0; i < length; i++){	
		int temp = 1;
		for(int j = 1; j <=i; j++){
			temp = temp*10;
		}
		sum = sum + temp*binseq[length - 1 - i];
	}
	return sum;
}
/*---------------------------------
FUNCTION:
	void SCEncoding(int *y, int *output, set<int>frozenIndex, int *frozenBit, double **W)
DESCRIPTION:
	This functin implements the SC Encoding according to Korada and Urbanke 's IEEE Trans paper.
PARAMETERS:
	INPUT:
		y--received bits.
		output--
		frozenIndex--
		frozenBit--
		W--channel information
	OUTPUT:
		output--
RETURN VALUES:
	None	
----------------------------------*/
void SCEncoding(int *received, int *decoded, set<int>FrozenIndex, int*FrozenBit, double**W){
   int FIndex = 0;
   int m = (int)log2((double)N);
   vector<vector<long double> >P0, P1;
   vector<vector<int> >C0, C1;

   for(int lambda = 0; lambda < m + 1; lambda++){
	vector<int>row1;
	vector<long double>row2;
	for(int j = 0; j < (int)pow(2.0, m - lambda ); j++){
		row1.push_back(0);
		row2.push_back(0.0);
	}
     P0.push_back(row2);
     P1.push_back(row2);
     C0.push_back(row1);
     C1.push_back(row1);
   }

   for(int beta = 0; beta < N; beta ++){//Initialization
     P0[0][beta] = W[0][received[beta]];
     P1[0][beta] = W[1][received[beta]];
   }
  
   for(int varphi = 0; varphi < N; varphi ++){//Main loop
     //cout<<"now, varphi is "<<varphi<<endl;
     recursivelyCalcPSE(m, varphi, P0, P1, C0, C1);
     if(FrozenIndex.find(varphi) != FrozenIndex.end()){//frozen bit
       //cout<<varphi<< " is Frozen bit"<<endl;
       if (varphi%2 == 0)
          C0[m][0] = FrozenBit[FIndex++];
       else
	  C1[m][0] = FrozenBit[FIndex++];
     }
     else{
       srand(time(NULL));
       //cout<<"Now being 0 is "<<P0[m][0]<<"; being 1 is "<<P1[m][0]<<" "<<endl;
	long double threshold = P0[m][0]/(P0[m][0] + P1[m][0]);
	//cout<<"Threshold is "<<threshold<<endl;
	long double ran =  (long double)rand()/RAND_MAX ;
	//	cout<<"random number is "<<ran<<endl;
       if(ran <= threshold){
	 //cout<<"set to 0"<<endl;
	 if (varphi%2 == 0)
	   C0[m][0] = 0;
	 else
           C1[m][0] = 0;
       }
       else {
	 // cout<<"set to 1"<<endl;
	 if (varphi %2 == 0)
	    C0[m][0] = 1;
	 else
	   C1[m][0] = 1;
       }
     }
     // cout<<"-------------------"<<endl;
     if(varphi%2 == 1){
       recursivelyUpdateC(m, varphi, C0, C1);
     }
   }

   for(int beta = 0; beta < N; beta++){
     decoded[beta] = C0[0][beta];
   }
 }
/*---------------------------------------
FUNCTION: 
	void SCDecoding(int *u, int *decodedresult)
DESCRIPTION:
	This functin implements SC Decoding.
PARAMETERS:
	INPUT:
		u--
		decodedresult--
	OUTPUT:
		decodedresult--
RETURN VALUES:
	None
---------------------------------------*/
void SCDecoding(int *u, int *decodedresult){
	encode(decodedresult, u);
}
/*------------------------------------
FUNCTION:
	double sourceCodingVerify(double **W)
DESCRIPTION:
	This function verifies source coding process described by Korada and Urbanke in IEEE Trans.
PARAMETERS:
	INPUT:
		W--channel information
		M--Measuring distoration
		frozenIndex--
	OUTPUT:
		average distortation
RETURN VALUES:
	average distoration	
------------------------------------*/
long double sourceCodingVerify(double **W, double **M, set<int>frozenIndex){
	
	
	double distortion = 0;
	int frozenBit[N-K], Message[N], reproductioncodeword[N];

	for(int i = 0; i < trial; i++){
		/*-------Generate random message-----*/	
		for(int j = 0; j < N; j++)
			Message[j] = rand()%2;
		int frozenbit = 0; 
		int2bin(frozenbit, frozenBit, N - K); 
		/*-------Map the message to a codeword of polar codes ------------*/
		SCEncoding(Message, reproductioncodeword, frozenIndex, frozenBit, W);
		/*-------Measuring distoration-------------*/
		distortion += measureDistortion(Message, reproductioncodeword, M);
		if(i%100 == 0)
		cout<<"so far distortion is "<<distortion/(i+1)<<endl;		 
	}	
	distortion /= trial;
	stringstream ss;
        ss<<N;
        string nstr = ss.str(); 
        string name = nstr + " Source Coding Result.txt";
	char *filename = new char[30];
        strcpy(filename, name.c_str());
        ofstream out(filename);
	out<<"-------------------------\n";
	out<<"Polar Code for Source Coding Experiment information\n";
	out<<"q is :"<<q<<endl;
	out<<"N is :"<<N<<endl;	
	out<<"K is :"<<K<<endl;
	out<<"I(W):" <<compute_I(W)<<endl;
	out<<"Rate (K/N):"<< (double) K/ (double)N<<endl;
	out<<"Average Distortion: "<<distortion<<endl;
	out<<"Theoretical bound for distortion: "<<compute_D(W,M)<<endl;
	out<<"Trial number: "<<trial<<endl;
	out<<"\n\n\n\n                             Qing Li(qingli@cse.tamu.edu)";
	return distortion;
}
