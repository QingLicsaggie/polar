#include"analysis.h"
#include"noisywem.h"
#include"channelPolarization.h"
#include"generator.h"
#include"listdecoding.h"
#include"rewritting.h"
#include"noise.h"
#include"SC.h"
#include<iostream>
#include<cstdio>
#include<fstream>
#include<sstream>
#include<cstring>
#include<cmath>
#include<set>
#include<vector>
#include<algorithm>
using namespace std;
extern int N;
extern int q;
extern long double theoretical_upper_bound;
extern int trial;
/*---------------------------
DESCRIPTION:
           Find out the Frozen index
           for noisyWEM: suppose
           the frozen index for channel
           is F1, the frozen index for
           WEM is F2, thus the frozen
           index for noisyWEM is F1 - F2
  ---------------------------*/
void noisyWEM::determineFrozenSet( char *fileforP, char *fileforW){
  vector<long double>probForP, probForW;
  readVP(probForP, fileforP);
  readVP(probForW, fileforW);
  selectChannels(probForP, frozenIndexForP, kForP);
  /*while(theoretical_upper_bound >= 0.00005 & KForP > 0){
    theoretical_upper_bound = 0;
    KForP--;
    selectChannels(PP, frozenIndexForP, KForP);
    }*/
  selectChannels(probForW, frozenIndexForW, kForW);
  set_difference(frozenIndexForW.begin(), frozenIndexForW.end(), frozenIndexForP.begin(), frozenIndexForP.end(), inserter(frozenIndex, frozenIndex.end()));
}
/*----------------------------------------------------------
DESCRIPTION:
           SCEncodingForNoisyWEM is a bit different than the
           SCEncoding as it involves to frozen sets. One is
           for WEM - channel, for which we set to the desired value, and
           the other is for channel, for which we set to 0. For
           the remaining parts, we set based on posterior probability.
	   For more details, please refer to algorithm 2 of my draft.
  -----------------------------------------------------------*/
void noisyWEM::SCEncodingForNoisyWEM(int *received, int *result, int *frozenBit, double **WForW){
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
     P0[0][beta] = WForW[0][received[beta]];
     P1[0][beta] = WForW[1][received[beta]];
   }
  
   for(int varphi = 0; varphi < N; varphi ++){//Main loop
     recursivelyCalcPSE(m, varphi, P0, P1, C0, C1);
     if(frozenIndex.find(varphi) != frozenIndex.end()){//frozen bit
       if (varphi%2 == 0)
          C0[m][0] = frozenBit[FIndex++];
       else
	  C1[m][0] = frozenBit[FIndex++];
     }
     else if(frozenIndexForP.find(varphi) != frozenIndexForP.end()){//frozen bit for channel. This part is different from the traditional one
	 if(varphi%2 == 0)
	   C0[m][0] = 0;
	 else
	   C1[m][0] = 0;
      }
     else{
	long double threshold = P0[m][0]/(P0[m][0] + P1[m][0]);
	long double ran = (long double)rand()/RAND_MAX ;
       if(ran <= threshold){
	 if (varphi%2 == 0)
	   C0[m][0] = 0;
	 else
           C1[m][0] = 0;
       }
       else {
	 if (varphi %2 == 0)
	    C0[m][0] = 1;
	 else
	   C1[m][0] = 1;
       }
     }
     if(varphi%2 == 1){
       recursivelyUpdateC(m, varphi, C0, C1);
     }
   }
	
   for(int beta = 0; beta < N; beta++){
     result[beta] = C0[0][beta];
   }
} 
/*-----------------------------------------------------------
DESCRIPTION:
            Please refer to algorithm 2 of my draft for details.
  ------------------------------------------------------------*/
void noisyWEM:: rewritterForNoisyWEM(int *current, int *dither, int *updated, int *frozenBit, double **WForW){

  int datavector[N];
  int tempUpdated[N];

  for(int i = 0; i < N; i++){
    datavector[i] = rand()%2;
    dither[i] = (current[i] + datavector[i])%2;
  }
  
  SCEncodingForNoisyWEM(datavector, tempUpdated, frozenBit, WForW);

  for(int i = 0; i < N; i++)
    updated[i] = (tempUpdated[i] + dither[i])%2;
}

/*------------------------------------------------------
DESCRIPTION:
     decoding operation.
     Please refer to algorithm 3 of my draft
  -----------------------------------------------------*/
void noisyWEM:: decoderForNoisyWEM(int *received, int *dither, int *result,  double **WForP){
  int data[N];
  for(int i = 0; i < N; i++)
    data[i] = (received[i] + dither[i])%2;
  int frozen[N - kForP];
  for(int i = 0; i < N - kForP; i++)
    frozen[i] = 0;

  listDecoder ld(getL()); //I use list decoder
  ld.SCLDecoder(data, result, frozenIndexForP, frozen, WForP);
}
/*-------------------------------------------------
DESCRIPTION:
          Compute decoding error. It has nothing to do with the information set
 --------------------------------------------------*/
int noisyWEM:: computeError(int *data, int *decoded){
  int messageDecoded[frozenIndex.size()];
  findFrozenBits(decoded, messageDecoded, frozenIndex);
  for(int i = 0; i < (int)frozenIndex.size(); i++){
	if(messageDecoded[i] != data[i])
		return 1;
  } 
    return 0;
}
/*-----------------------------------------------
DESCRIPTION:
          Print experiment results
  -------------------------------------------------*/
void noisyWEM:: print_info_for_noisyWEM(double **WForW, double **WForP, double **costM, int count){
  stringstream ss;
  ss<<N;
  stringstream bb;
  bb<<wemChannelP;
  string bstr = bb.str();
  string nstr = ss.str();
  string name = nstr + " " + bstr + "noisyWEM result(max).txt";
  char *filename = new char [30];
  strcpy(filename, name.c_str());
  ofstream out(filename);
  out<<"polar codes for NoisyWEM experiment information"<<endl;
  out<<"q: "<< q<<endl;
  out<<"N is :"<<N<<endl;
  out<<"The probability for WEM   Channel is "<<wemChannelP<<endl;
  out<<"The probability for Noise Channel is "<<channelP<<endl;
  out<<"Decoding method is list decoding with L = "<<getL()<<endl;
  out<<"Frozen size: "<<frozenIndex.size()<<endl;
  out<<"KForP "<<kForP<<endl;
  out<<"kForW "<<kForW<<endl;
  out<<"kForP - KForW: "<<kForP - kForW<<endl;
  out<<"Theoretical lower bound: "<< (compute_I(WForP) - compute_I(WForW))<<endl;
  out<<"Theoretical rewriting cost: "<<compute_D(WForW, costM)<<endl;
  out<<"Rewriting rate: "<<(double)frozenIndex.size()/N<<endl;
  out<<"Block error rate: "<<error<<endl;
  out<<"Rewriting cost "<<cost<<endl;
  out<<"Percentage is "<<(double)count/trial<<endl;
  out<<"                           Qing Li(qingli@cse.tamu.edu)"<<endl;
  out<<"------------------------"<<endl;
  delete[] filename;
}
/*------------------------------------------------
DESCRIPTION:
          Experiment of noisyWEM
PARAMETERS:
          WForP - probability information for channel
          WForW - probability information for WEM
  -------------------------------------------------*/
void noisyWEM::experiment_of_noisyWEM(char *fileforP, char *fileforW, int trial, double **WForP, double **WForW,  double **costM ){

  determineFrozenSet(fileforP, fileforW);
  int size = frozenIndex.size();
  int current[N], dither[N], updated[N], corruptedcodeword[N], decodedresult[N];
  for(int i = 0; i < N; i++)
    current[i] = rand()%2; 
  int count = 0;
  for(int i = 0; i < trial; i++){
    int data[size];
    for(int j = 0; j < size; j++)
	data[j] = rand()%2;
    rewritterForNoisyWEM(current, dither, updated, data, WForW);
    cost += measureDistortion(current, updated, costM); // measureL1(current, updated); 
    if(measureDistortion(current, updated, costM) > 1.1*compute_D(WForW, costM))
	count++;
    BSC_corrupted(updated, corruptedcodeword, WForP);

   //update current    
    for(int j = 0; j < N; j++)
    	current[j] = corruptedcodeword[j];

    //decoder
    decoderForNoisyWEM(corruptedcodeword, dither, decodedresult, WForP);  
    error += computeError(data, decodedresult); 
    if(i%500==0){
	cout<<"Now, i is "<<i<<endl;
	cout<<"so far block error rate is "<<error/(i + 1)<<endl;
    	cout<<"so far average rewriting cost "<<cost/(i + 1)<<endl;
        cout<<"so far the probability not in the range is "<<(double)count/(i+1)<<endl;
    	cout<<"-----------------------"<<endl;
  	print_info_for_noisyWEM(WForW, WForP, costM, count);	
	}
  }
  error /= trial;
  cost /= trial;
  print_info_for_noisyWEM(WForW, WForP, costM, count);
}
