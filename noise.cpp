#include<stdio.h>
#include<time.h>
#include<stdlib.h>
#include<iostream>
#include<sstream>
#include<cstring>
#include<fstream>
#include<cmath>
#include<chrono>
#include<random>
#include<set>
using namespace std;
extern int N;

/*--------------------------
FUNCTION:
	void BSC_corrupted(int *codeword, int *corruptedcodeword, double **W)
DESCRIPTION:
	This function corruptes the codeword according to BSC channel.
PARAMETERS:
	INPUT:
		codeword--
		corruptedcodeword--
		W--Channel information
	OUTPUT:
		corruptedcodeword
RETURN VALUES:
	None
---------------------------*/
void BSC_corrupted(int *codeword, int *corruptedcodeword, double **W){
	for(int i = 0; i < N; i++){
		double p = (double)rand()/RAND_MAX;
		if(p < W[0][1]){
			corruptedcodeword[i] = (1 + codeword[i])%2;
			//cout<<"flip location "<<i<<endl;
		}
		else
			corruptedcodeword[i] = codeword[i];	
	}
	//cout<<endl;
}

/*---------------------------------------------
FUNCTION:
      flipBits(int *codeword, int *corruptedcodeword, int FlipNum)
DESCRIPTION:
       Given codeword, we flip FlipNum bits, and the result is stored in corruptedcodeword.
PARAMETERS:
       INPUT:
           codeword--received codeword.
	   corruptedcodeword--noisy codeword.
	   FlipNum          --number of bits we flip
       OUTPUT:
           corruptedcodeword
RETURN VALUES:
       NONE
LAST UPDATE DATE:
      05/28/2014
  ----------------------------------------------*/
void flipBits(int *codeword, int *corruptedcodeword, int FlipNum){
  for(int i = 0; i < N; i++)
    corruptedcodeword[i] = codeword[i];
  int count = 0;
  set<int>flip;
  while(count != FlipNum){
    int location = rand()%N;;
    if(flip.find(location) == flip.end()){
      count++;
      //cout<<"we have flipped the "<<location<<" -th bits "<<endl;
      corruptedcodeword[location] = (1+corruptedcodeword[location])%2;
    }    
  }
}

int countFlip(int *codeword, int *noisycodeword){
  int count = 0;
  for(int i = 0; i < N; i++){
    if(codeword[i] != noisycodeword[i]){
      count++;
      //cout<<"location "<<i<<endl;
     }
  }
  return count;
}

/*---------------------------------------------
FUNCTION: 
       void mgrns(double mean, double sigma, double seed, int n, double *a)
DESCRIPTION:
       Given gaussian distribution parameters mean and sigma,
       generate a sequence of numbers with gaussian distribution.
PARAMETERS:
       INPUT:
           mean--
	   sigma--
	   seed--
	   n--codeword length
	   a--storing result
       OUTPUT:
           a
RETURN VALUE:
       NONE
LAST UPDATE DATE:
      06/17/2014
  -----------------------------------------------*/
void mgrns(double mean, double sigma, int n, double *a){
  //construct a trivial random generator engine from a time-based seed
  unsigned seed = chrono::system_clock::now().time_since_epoch().count();
  default_random_engine generator(seed);
  
  normal_distribution<double> distribution(mean, sigma);
  for(int i  = 0; i < n; i++)
    a[i] = distribution(generator);
}
/*------------------------------------------------------
FUNCTION:
       void AWGN(double *original, double *noisy, double sigma, int length)
DESCRIPTION:
        Given the original codeword sequence and parameters for the gaussian
	parameters, generate noisy codeword following with gaussian distribution.
PARAMETERS:
        INPUT:
	     original--received codeword
	     noisy   --noisy codeword
	     sigma   --mean for gaussian is 0
	     length  --codeword length
        OUTPUT:
	     noisy
  -------------------------------------------------------*/
void AWGN(double *original, double *noisy, double sigma, int length){
  
  double gaussianNoise[length];
  mgrns(0, sigma, length, gaussianNoise);
 
  for(int i = 0; i < length; i++){
    //cout<<"noise "<<gaussianNoise[i]<<" and codeword "<<original[i];
    noisy[i] = original[i] + gaussianNoise[i];
    //cout<<" noise codeword is "<<noisy[i]<<endl;
  }
}
/*------------------------------------------------------------
FUNCTION:
        void AWGNModel(double sigma, int num)
DESCRIPTION:
        generate a approximate awgn model
PARAMETERS:
        INPUT:
	   sigma--
	   num--quantize number
	OUTPUT:
	   result in a txt file
RETURN VALUES:
       NONE
LAST UPDATE DATE:
       19/06/2014
  ----------------------------------------------------------------*/
void AWGNModel(double sigma, int num){
  const int nrolls = 10000; //number of experiments
  stringstream ss;
  ss<<sigma;
  string nstr = ss.str();
  string name = "Gaussian " + nstr+ ".txt";
  char *filename = new char[300];
  strcpy(filename, name.c_str());
  ofstream out(filename);

  default_random_engine generator;

  normal_distribution<double>distribution0(1, sigma);
  
  int p0[2] = {0};
  for(int i = 0; i < nrolls; i++){
    double number0 = distribution0(generator);
    if(number0>=0)
      p0[0]++;
    else 
      p0[1]++;
  }
  out<<(double)p0[0]/nrolls<<endl;
  out<<(double)p0[1]/nrolls<<endl;
  out<<(double)p0[1]/nrolls<<endl;
  out<<(double)p0[0]/nrolls<<endl;
  delete []filename;
}

/*-------------------------------------------------
FUNCTION:
       void AWGNMatrix(double **matrix, double sigma)
DESCRIPTION:
        Createa AWGN 2*2 matrix
PARAMETERS:
        INPUT
	    matrix--storing result
	    sigma --gaussian parameter
        OUTPUT:
	    matrix
RETURN VALUES;
        NONE
LAST UPDATE DATE:
       06/19/214
-------------------------------------------------*/
void AWGNMatrix(double **matrix, double sigma){
  const int nrolls = 10000; //number of experiments
  default_random_engine generator;
  normal_distribution<double>distribution0(1, sigma);
  
  int p0[2] = {0};
  for(int i = 0; i < nrolls; i++){
    double number0 = distribution0(generator);
    if(number0>=0)
      p0[0]++;
    else 
      p0[1]++;
  }

  matrix[0][0] = matrix[1][1] = (double)p0[0]/nrolls;
  matrix[0][1] = matrix[1][0] = 1 - matrix[0][0];
}
/*----------------------------------------------------
FUNCTION:
         int quantize(double receive, double sigma, int length)
DESCRIPTION:
         Given received value, receive, the quantize number and
	 the step value, return the quantize value.
PARAMETERS:
         INPUT:
	     receive -- received value
	     sigma   -- step value
	     length  -- quantize length
	 OUTPUT:
	      quantize number
RETURN VALUES:
          quantize number
LAST UPDATE DATE:
          06/19/2014
  -----------------------------------------------------*/

int quantize(double receive, double sigma, int length){
  double upper = sigma*length/2;
  double lower = -sigma*length/2;
  if(receive>upper)
    return length - 1;
  else if(receive<lower)
    return 0;
  else
    return ceil((receive - lower)/sigma);      
}

/*------------------------------------------
FUNCTION:
        double quantizeProbBeing1(int num, int length)
DESCRIPTION:
        Compute the probability of being 1.
  ------------------------------------------*/
double quantizeProbBeing1(int num, int length){
  return (double)(num+1)/length;
}

void quantizeMatrix(int **matrix, int length){
  for(int i = 0; i < length; i++)
    matrix[1][length-i] = matrix[0][i] = quantizeProbBeing1(i, length); 
}
