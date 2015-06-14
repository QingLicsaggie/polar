#include<cstdio>
#include<iostream>
#include<math.h>
#include"analysis.h"
#include<iostream>
#include<fstream>
#include<cstring>
#include<sstream>
#include<set>
using namespace std;
extern int q;
extern int N;
extern int K;
extern int trial;
extern double average_error_rate;
extern double theoretical_upper_bound;
extern double flipCount1;
/*---------------------------
FUNCTION:
	double compute_I(double **W)
DESCRIPTION:
	This function computes I(W) according to Arikan's paper
PARAMETERS:
	INPUT:
		W--channel information
	OUTPUT:
		I(W)
RETURN VALUES:
	I(W)
----------------------------*/
double compute_I(double **W){
	
	int y, x;
	double I = 0;

	for(y = 0; y < 2; y++)
		for(x = 0; x < 2; x++)
			I += 0.5*W[x][y]*log2(W[x][y]/(0.5*W[0][y] + 0.5*W[1][y]))/log2(2.0);			
	return I;
}
/*----------------------------
FUNCTION:
	double errorBitRate(int *originalMessage, int *decodedMessage)
PARAMETERS:
	INPUT:
		originalMessage--
		decodedMessage--
	OUTPUT:
		error bit
RETURN VALUES:
	error bit

------------------------------*/

double errorBitRate(int *originalMessage, int *decodedMessage){
	int errorNum = 0, i;
	for(i = 0; i < N; i++)
		if(originalMessage[i] != decodedMessage[i])			
			errorNum ++;
	return (double)errorNum/(double)N; 
}

double blockErrorRate(int *originalMessage, int *decodedMessage){
	if (N*errorBitRate(originalMessage, decodedMessage) > 0)
		return 1;
	else
		return 0;
}


/*---------------------------------------------------------
FUNCTION:
long double computeTheoreticalBER(set<int>frozenIndex, vector<long double>Probablity)
DESCRIPTION:
      Computes theoretical BER
PARAMETERS:
       INPUT:
         frozenIndex --
         Probability --
       OUTPUT
          upper bound on the BLOK ERROR RATE
RETURN VALUES:
         Theoretical BLOCK ERROR RATE
-----------------------------------------------------------*/
long double computeTheoreticalBER(set<int>frozenIndex, vector<long double>Probability){
  set<int>::iterator it;
  long double BER = 0.0;
  for(int i = 0; i < N; i++){
    it = frozenIndex.find(i);
    if(it == frozenIndex.end())
      BER += Probability[i];
  }
  return BER;
}
/*--------------------------
FUNCTION:
	void print_info(set<int>frozenSet, vector<long double>Probability, double **W)
DESCRIPTION:
	This function prints some information about Polar code
PARAMETERS:
	INPUT:
	  frozenSet--
	  Probability--
	  W--
	OUTPUT:
RETURN VALUES:
	None
----------------------------*/

void print_info(set<int>frozenSet, double prob, double prob1 ){
	stringstream ss;
        ss<<N<<" "<<prob<<" ";   
        string nstr = ss.str(); 
        string name = nstr;
	char *filename = new char[100];
        strcpy(filename, name.c_str());
        ofstream out(filename);
	out<<"-------------------------\n";
	out<<"Polar Code Experiment information\n";
	out<<"Original flip rate is "<<prob<<endl;
	out<<"Actual flip rate is "<<prob1<<endl;
	out<<"Actual flip rate in experiment is "<<flipCount1;
	out<<"Rate (K/N):"<< (double) K/ (double)N<<endl;
	out<<"BER obtained by experiment:"<<average_error_rate<<endl;
	out<<"Theoretical uppper bounded is "<<theoretical_upper_bound<<endl;
	out<<"trial is: "<<trial<<endl;
	delete []filename;
}


void concatePrint_info(set<int>frozenSet, double prob, double prob1 ){
	stringstream ss;
        ss<<N<<" "<<prob<<" concate";   
        string nstr = ss.str(); 
        string name = nstr;
	char *filename = new char[100];
        strcpy(filename, name.c_str());
        ofstream out(filename);
	//out<<"-------------------------\n";
	//out<<"Polar Code Experiment information\n";
	//out<<"Original flip rate is "<<prob<<endl;
	//out<<"Actual flip rate is "<<prob1<<endl;
	//out<<"Actual flip rate in experiment is "<<flipCount1;
	out<<"Rate (K/N):"<< (double) K/ (double)N<<endl;
	out<<"BER obtained by experiment:"<<average_error_rate<<endl;
	//out<<"Theoretical uppper bounded is "<<theoretical_upper_bound<<endl;
	//out<<"trial is: "<<trial<<endl;
	delete []filename;
}
