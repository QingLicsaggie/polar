#include<iostream>
#include<vector>
#include<cmath>
#include<set>
#include<sys/time.h>
#include<sys/resource.h>
#include"channelPolarization.h"
#include"listdecoding.h"
#include"generator.h"
#include"crc16.h"
using namespace std;
extern int N;
extern int K;
/*-----------------------------------
 FUNCTION:
         void recursivelyUpdateB (int lambda, int varphi, vector<int> &result)
 DESCRIPTION:
          Algorithm 4.
          Note that varphi must be odd, otherwise there is no need.
 PARAMETERS:
       INPUT:
         lambda --layer
         varphi --phase information
 	 result --used to store the final result
       OUTPUT:
          result
 RETURN VALUES:
        None
   
   ------------------------------------*/
void  recursivelyUpdateB(int lambda, int varphi, vector<int> &result){
   if(varphi%2==0) //varphi must be odd
     return;
   else{
     int phi = varphi/2;
     int upper = N/(int)pow(2.0, lambda);
     for(int beta = 0; beta < upper; beta++){ 
	  result[(lambda - 1)*N + phi + 2*beta*(int)pow(2.0, lambda - 1)] = (result[lambda*N + varphi - 1 + beta *(int)pow(2.0, lambda)] + result[lambda*N + varphi + beta*(int)pow(2.0,lambda)])%2;

       result[(lambda - 1)*N + phi + (2*beta + 1)*(int)pow(2.0, lambda -1)] = result[lambda*N + varphi + beta*(int)pow(2.0,lambda)];
     }

     if(phi%2 == 1)
       recursivelyUpdateB(lambda - 1, phi, result);
   }
 }

/*-----------------------------------
 FUNCTION:
         void recursivelyUpdateC (int lambda, int varphi, vector<vector<int> *> &C0, vector<vector<int> *> &C1 )
 DESCRIPTION:
          Algorithm 7.
          Note that varphi must be odd, otherwise there is no need.
 PARAMETERS:
       INPUT:
         lambda --layer
         varphi --phase information
 	 C1, C0 --used to store the final result
       OUTPUT:
          C1, C0
 RETURN VALUES:
        None
   
   ------------------------------------*/
void recursivelyUpdateC(int lambda, int varphi, vector<vector<int> > &C0, vector<vector<int> >&C1){
   if(varphi%2==0) //varphi must be odd
     return;
   else{
     int phi = varphi/2;
     int upper = N/(int)pow(2.0, lambda);
     for(int beta = 0; beta < upper; beta++){ 
       if(phi%2 == 0){
        C0[lambda - 1][2*beta] = (C0[lambda][beta] + C1[lambda][beta])%2;
	C0[lambda - 1][2*beta + 1] = C1[lambda][beta];
       }
       else{
        C1[lambda - 1][2*beta] = (C0[lambda][beta] + C1[lambda][beta])%2;
	C1[lambda - 1][2*beta + 1] = C1[lambda][beta];
       }
     }
     if(phi%2 == 1)
       recursivelyUpdateC(lambda - 1, phi, C0, C1);
   }
 }
/*------------------------------------
 FUNCTION:
        void recursivelyCalcP(int lambda, int varphi, vector<long double>&P0,  vector<long double> &P1, vector<int> B)
 DESCRIPTION:
        Algorithm 3
 PARAMETERS:
        INPUT:
              lambda  -- layer
              varphi  -- phase
              P0      -- used to store result when the input value is 0
              P1      -- used to store result when the input value is 1
	      B       -- input bit information
        OUTPUT:
              P0--
              P1--
 RETURN VALUES
       None
       ----------------------------------*/
   void recursivelyCalcP(int lambda, int varphi, vector<long double>&P0, vector<long double>&P1, vector<int>B){
   if(lambda == 0){//stoping condition 
     return;
   }
   long double sigma = 0;
   int phi = varphi/2;
   int m = (int)log2((double)N);
   if(varphi %2 == 0){ 
     recursivelyCalcP(lambda - 1, phi, P0, P1, B);}

   for(int beta = 0; beta < (int)pow(2.0, m - lambda); beta ++){
     if(varphi%2 == 0){//apply equation 4
       P0[lambda*N + varphi + beta*(int)pow(2.0, lambda)] = 0.5*(P0[(lambda - 1)*N + phi + 2*beta*(int)pow(2.0, lambda - 1)]*P0[(lambda - 1)*N + phi + (2*beta + 1)*(int)pow(2.0, lambda - 1)] + P1[(lambda - 1)*N  + phi + 2*beta*(int)pow(2.0, lambda - 1)]*P1[(lambda - 1)*N + phi + (2*beta + 1)*(int)pow(2.0, lambda - 1)]);
       P1[lambda*N + varphi + beta*(int)pow(2.0, lambda)] = 0.5*(P1[(lambda - 1)*N + phi + 2*beta*(int)pow(2.0, lambda - 1)]*P0[(lambda - 1)*N + phi + (2*beta + 1)*(int)pow(2.0, lambda - 1)] + P0[(lambda - 1)*N  + phi + 2*beta*(int)pow(2.0, lambda - 1)]*P1[(lambda - 1)*N + phi + (2*beta + 1)*(int)pow(2.0, lambda - 1)]);
     }
     else{ //apply equation 5
       int u = B[lambda*N + varphi - 1 + beta*(int)pow(2.0, lambda)];
	
       if(u == 0){
	 P0[lambda*N + varphi + beta*(int)pow(2.0, lambda)] = 0.5*P0[(lambda - 1)*N + phi + 2*beta*(int)pow(2.0, lambda - 1)]*P0[(lambda - 1)*N + phi + (2*beta + 1)*(int)pow(2.0, lambda - 1)];
	 P1[lambda*N + varphi + beta*(int)pow(2.0, lambda)] = 0.5*P1[(lambda - 1)*N + phi + 2*beta*(int)pow(2.0, lambda - 1)]*P1[(lambda - 1)*N + phi + (2*beta + 1)*(int)pow(2.0, lambda - 1)];
       }
       else{
	 P0[lambda*N + varphi + beta*(int)pow(2.0, lambda)] = 0.5*P1[(lambda - 1)*N + phi + 2*beta*(int)pow(2.0, lambda - 1)]*P0[(lambda - 1)*N + phi + (2*beta + 1)*(int)pow(2.0, lambda - 1)];
	 P1[lambda*N + varphi + beta*(int)pow(2.0, lambda)] = 0.5*P0[(lambda - 1)*N + phi + 2*beta*(int)pow(2.0, lambda - 1)]*P1[(lambda - 1)*N + phi + (2*beta + 1)*(int)pow(2.0, lambda - 1)];
       }
     }
	if (P0[lambda*N + varphi + beta*(int)pow(2.0, lambda)]  > sigma) sigma = P0[lambda*N + varphi + beta*(int)pow(2.0, lambda)] ;
	if (P1[lambda*N + varphi + beta*(int)pow(2.0, lambda)]  > sigma) sigma = P1[lambda*N + varphi + beta*(int)pow(2.0, lambda)];
    }
   
	for(int beta = 0; beta < (int)pow(2.0, m - lambda); beta ++){
		P0[lambda*N + varphi + beta*(int)pow(2.0, lambda)] = P0[lambda*N + varphi + beta*(int)pow(2.0, lambda)]/sigma;
		P1[lambda*N + varphi + beta*(int)pow(2.0, lambda)] = P1[lambda*N + varphi + beta*(int)pow(2.0, lambda)]/sigma;	
	}
}
 /*------------------------------------
 FUNCTION:
        void recursivelyCalcPSE(int lambda, int varphi, vector<vector<long double>* >&P0,  vector<vector<long double> *> &P1, vector< vector<int> *> &C0, vector<vector<int> *> &C1)
 DESCRIPTION:
        Algorithm 6
 PARAMETERS:
        INPUT:
              lambda  -- layer
              varphi  -- phase
              P0      -- used to store result when the input value is 0
              P1      -- used to store result when the input value is 1
	      C0, C1       -- input bit information
        OUTPUT:
              P0--
              P1--
 RETURN VALUES
       None
       ----------------------------------*/
void recursivelyCalcPSE(int lambda, int varphi, vector<vector<long double> >&P0, vector<vector<long double> >&P1, vector<vector<int> >&C0, vector<vector<int> >&C1){
   if(lambda == 0){//stoping condition 
     return;
   }
   long double sigma = 0;
   int phi = varphi/2;
   int m = (int)log2((double)N);
   if(varphi %2 == 0){ 
     recursivelyCalcPSE(lambda - 1, phi, P0, P1, C0, C1);}

   for(int beta = 0; beta < (int)pow(2.0, m - lambda); beta ++){
     if(varphi%2 == 0){//apply equation 4
       P0[lambda][beta] = 0.5*(P0[lambda - 1][2*beta]*P0[lambda - 1][2*beta + 1] + P1[lambda - 1][2*beta]*P1[lambda - 1][2*beta + 1]);
       P1[lambda][beta] = 0.5*(P1[lambda - 1][2*beta]*P0[lambda - 1][2*beta + 1] + P0[lambda - 1][2*beta]*P1[lambda - 1][2*beta + 1]);
     }
     else{ //apply equation 5	
       int u = C0[lambda][beta];
       if(u == 0){
	 P0[lambda][beta] = 0.5*P0[lambda - 1][2*beta]*P0[lambda - 1][2*beta + 1];
	 P1[lambda][beta] = 0.5*P1[lambda - 1][2*beta]*P1[lambda - 1][2*beta + 1];
       }
       else{
	 P0[lambda][beta] = 0.5*P1[lambda - 1][2*beta]*P0[lambda - 1][2*beta + 1];
	 P1[lambda][beta] = 0.5*P0[lambda - 1][2*beta]*P1[lambda - 1][2*beta + 1];
       }
     }
	if (P0[lambda][beta]  > sigma) sigma = P0[lambda][beta];
	if (P1[lambda][beta]  > sigma) sigma = P1[lambda][beta];
    }
   
	for(int beta = 0; beta < (int)pow(2.0, m - lambda); beta ++){
		P0[lambda][beta] = P0[lambda][beta]/sigma;
		P1[lambda][beta] = P1[lambda][beta]/sigma;	
	}

}
/*------------------------------------------------
FUNCTION:
        void SCFirst(int * received, int *decoded, set<int>FrozenIndex, int* FrozenBit, double **W)
Description:
       Algorithm2
PARAMETERS:
       INPUT:
            received --received vector;
	    decoded  --deoded result;
	    frozenIndex
	    frozenBit
            W        --channel information
	    G        --generator
      OUTPUT:
            decoded
RETURN VALUES:
     None
 
 -------------------------------------------------*/
void SCFirst(int *received, int *decoded, set<int>FrozenIndex, int*FrozenBit, double**W){
   
   int FIndex = 0;
   int m = (int)log2((double)N);
   vector<long double>P0((1+m)*N), P1((1+m)*N);
   vector<int>B((1+m)*N);
   for(int beta = 0; beta < N; beta ++){//Initialization
     P0[beta] = W[0][received[beta]];
     P1[beta] = W[1][received[beta]];
   }
   for(int varphi = 0; varphi < N; varphi ++){//Main loop
     recursivelyCalcP(m, varphi, P0, P1, B);
     if(FrozenIndex.find(varphi) != FrozenIndex.end()){//frozen bit
       B[m*N + varphi] = FrozenBit[FIndex];
       decoded[varphi] = FrozenBit[FIndex];
       FIndex++;
     }
     else{	
       if(P0[m*N + varphi] > P1[m*N + varphi]){
	 B[m*N + varphi] = 0;
	 decoded[varphi] = 0;
	}
       else{
	 B[m*N + varphi] = 1;
         decoded[varphi] = 1;
	}
     }
     if(varphi%2 == 1){
       recursivelyUpdateB(m, varphi, B);
     }
   }
   int temp[N];
   for(int beta = 0; beta < N; beta++){
     temp[beta] = B[beta];
   }
   encode(decoded, temp);
 }

/*------------------------------------------------
FUNCTION:
        void SCFirstSE(int * received, int *decoded, set<int>FrozenIndex, int* FrozenBit, double **W)
Description:
       Algorithm5
PARAMETERS:
       INPUT:
            received --received vector;
	    decoded  --deoded result;
	    frozenIndex
	    frozenBit
            W        --channel information
      OUTPUT:
            decoded
RETURN VALUES:
     None
 
 -------------------------------------------------*/
void SCFirstSE(int *received, int *decoded, set<int>FrozenIndex, int*FrozenBit, double**W){
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
	//cout<<"varphi is "<<varphi<<endl;
     recursivelyCalcPSE(m, varphi, P0, P1, C0, C1);
     if(FrozenIndex.find(varphi) != FrozenIndex.end()){//frozen bit
	//cout<<"frozen bit"<<endl;
	decoded[varphi] = FrozenBit[FIndex];
       if (varphi%2 == 0)
          C0[m][0] = FrozenBit[FIndex++];
       else
	  C1[m][0] = FrozenBit[FIndex++];
     }
     else{
	/*cout<<"information bit"<<endl;
	 cout<<"P0 is "<<P0[m][0]<<endl;
	 cout<<"P1 is "<<P1[m][0]<<endl;*/
       if(P0[m][0] > P1[m][0]){
	  decoded[varphi] = 0;
	 if (varphi%2 == 0)
	   C0[m][0] = 0;
	 else
           C1[m][0] = 0;
       }
       else {
	 decoded[varphi] = 1;
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
 
   int temp[N];
   for(int beta = 0; beta < N; beta++){
     temp[beta] = C0[0][beta];
  }
   encode(decoded, temp);
   /*cout<<"decoding result is "<<endl;
   for(int i = 0; i < N; i++)
	cout<<decoded[i]<<" ";
   cout<<endl;*/

 }


/*------------------------------------------------
FUNCTION:
        void SCFirstSERevise(int * received1, int *received2, int *decoded, set<int>FrozenIndex, int* FrozenBit, double **W1, double **W2)
Description:
       Channel decoding when we have two codewords
PARAMETERS:
       INPUT:
            received1 --received vector 1;
            received2 --received vector 2;
	    decoded  --deoded result;
	    frozenIndex
	    frozenBit
            W1        --channel information 1
	    W2        --channel information 2
      OUTPUT:
            decoded
RETURN VALUES:
     None
 
 -------------------------------------------------*/
void SCFirstSERevised(int *received1, int *received2,  int *decoded, set<int>FrozenIndex, int*FrozenBit, double**W1, double **W2){
   int FIndex = 0;
   int m = (int)log2((double)N);
   vector<vector<long double> >P0, P1, PP0, PP1;
   vector<vector<int> >C0, C1, CC0, CC1;
   
   for(int lambda = 0; lambda < m + 1; lambda++){
	vector<int>row1;
	vector<long double>row2;
	for(int j = 0; j < (int)pow(2.0, m - lambda ); j++){
		row1.push_back(0);
		row2.push_back(0.0);
	}
     P0.push_back(row2);
     P1.push_back(row2);
     PP0.push_back(row2);
     PP1.push_back(row2);
     C0.push_back(row1);
     C1.push_back(row1);
     CC0.push_back(row1);
     CC1.push_back(row1);
   }

   for(int beta = 0; beta < N; beta ++){//Initialization
     P0[0][beta]  = W1[0][received1[beta]];
     P1[0][beta]  = W1[1][received1[beta]];
     PP0[0][beta] = W2[0][received2[beta]];
     PP1[0][beta] = W2[1][received2[beta]];
   }
   for(int varphi = 0; varphi < N; varphi ++){//Main loop
     recursivelyCalcPSE(m, varphi,  P0,  P1,  C0,  C1);
     recursivelyCalcPSE(m, varphi, PP0, PP1, CC0, CC1);
     if(FrozenIndex.find(varphi) != FrozenIndex.end()){//frozen bit
	decoded[varphi] = FrozenBit[FIndex];
	if (varphi%2 == 0)
         CC0[m][0] = C0[m][0] = FrozenBit[FIndex++];
       else
	 CC1[m][0] = C1[m][0] = FrozenBit[FIndex++];
     }
     else{
       if(P0[m][0]*PP0[m][0] > P1[m][0]*PP1[m][0]){
	  decoded[varphi] = 0;
	 if (varphi%2 == 0)
	   CC0[m][0] = C0[m][0] = 0;
	 else
           CC1[m][0] = C1[m][0] = 0;
       }
       else {
	 decoded[varphi] = 1;
	 if (varphi %2 == 0)
	   CC0[m][0] = C0[m][0] = 1;
	 else
	   CC1[m][0] = C1[m][0] = 1;
       }
     }
     if(varphi%2 == 1){
       recursivelyUpdateC(m, varphi, C0, C1);
       recursivelyUpdateC(m, varphi, CC0, CC1);
     }
   }
 
   int temp[N];
   for(int beta = 0; beta < N; beta++){
     temp[beta] = C0[0][beta];
  }
   encode(decoded, temp);
 }


//algorithm 8
void listDecoder::initializeDataStructures( ){
	int m = (int)log2((double)N);
	inactivePathIndices.reserve(L);
	activePath.resize(L);
	arrayReferenceCount.resize((m+1)*L);
	pathIndexToArrayIndex.resize((m+1)*L);
	
	for(int i = 0; i < (m+1)*L; i++){
	  vector<long double>row;
	  vector<int>B;
	  arrayPointer_P.push_back(row);
	  arrayPointer_C.push_back(B);
	}
	for(int i = 0; i < m + 1; i++){
	  vector<int>row;
	  row.reserve(L);
	  inactiveArrayIndices.push_back(row);
	}
	
	//Initialization of data structures
	for(int lambda = 0; lambda < m + 1; lambda++){
	  for(int s = 0; s < L; s++){
	    arrayPointer_P[lambda*L + s].resize((int)pow(2.0, m + 1 - lambda));
	    //cout<<"size of arrayPointer_P["<<lambda<<"]["<<s<<"] is "<<arrayPointer_P[lambda*L + s].size()<<endl;
	    arrayPointer_C[lambda*L + s].resize((int)pow(2.0, m + 1 - lambda));
	    arrayReferenceCount[lambda*L + s] = 0;
	    inactiveArrayIndices[lambda].push_back(s);
	  }
	}
	for(int l = 0; l < L; l++){
	  activePath[l] = false;
	  inactivePathIndices.push_back(l);
	}	
}
/*--------------------------------
Print out information for debuging
----------------------------------*/
void listDecoder::print(){
  cout<<"\nPrint out information:\n";
  int m = (int)log2((double)N);
  /*cout<<"Array Pointer P content: "<<endl;
  for(int lambda = 0; lambda < m + 1; lambda ++){
    for(int s = 0; s < L; s++){
      cout<<"arrayPointer_P[ "<<lambda<<" ][ "<<s<<" ]"<<endl;
      for(int i = 0; i < (int)pow(2.0, m - lambda); i++)
	cout<<arrayPointer_P[lambda*L + s][i]<<" ";
      cout<<endl;
      for(int i = 0; i < (int)pow(2.0, m - lambda); i++)
	cout<<arrayPointer_P[lambda*L + s][i + (int)pow(2.0, m - lambda)]<<" ";
      cout<<endl;
      cout<<"--------"<<endl;
    }
    }*/
  cout<<"Array Pointer C content: "<<endl;
  for(int lambda = 0; lambda < m + 1; lambda ++){
    for(int s = 0; s < L; s++){
      cout<<"arrayPointer_C[ "<<lambda<<" ][ "<<s<<" ]"<<endl;
      for(int i = 0; i < (int)pow(2.0, m - lambda); i++)
	cout<<arrayPointer_C[lambda*L + s][i]<<" ";
      cout<<endl;
      for(int i = 0; i < (int)pow(2.0, m - lambda); i++)
	cout<<arrayPointer_C[lambda*L + s][i + (int)pow(2.0, m - lambda)]<<" ";
      cout<<endl;
      cout<<"--------"<<endl;
    }
  }

  /*cout<<"arryReferenceCount content: "<<endl;
  for(int lambda = 0; lambda < m + 1; lambda ++){
    for(int s = 0; s < L; s++){
      cout<<"arrayReferenceCount["<<lambda<<"]["<<s<<"]: "<<arrayReferenceCount[lambda*L + s]<<endl;
    }
      cout<<"--------"<<endl;
  }

    cout<<"pathIndexToArrayIndex content:"<<endl;
    for(int lambda = 0; lambda < m + 1; lambda ++){
        for(int s = 0; s < L ; s++ )
            cout<<"pathIndexToArrayIndex["<<lambda<<"]["<<s<<" ] is "<<pathIndexToArrayIndex[lambda*L + s]<<endl;
    }
                   cout<<"------"<<endl;
  cout<<"inactiveArrayindices content:"<<endl;
  for(int i = 0; i < m + 1; i++){
    for(int j = 0; j < inactiveArrayIndices[i].size(); j++){
        cout<<"inactiveArrayIndices["<<i<<"]["<<j<<"] contains array "<<inactiveArrayIndices[i][j]<<" "<<endl;
    }
    cout<<"----------"<<endl;
  }
  cout<<"activePath content (stack, 0 is inactive, 1 is active):"<<endl;
  for(int i = 0; i < L; i++)
    cout<<"Path  "<<i<<": "<<activePath[i]<<endl;
  cout<<endl<<"-------------"<<endl;
		       
  cout<<"inactivePathIndices content (stack):"<<endl;
  for(int i = 0; i < inactivePathIndices.size(); i++)
    cout<<inactivePathIndices[i]<<" ";
    cout<<endl;*/
}
/*-------------------------------------------------------------------------
FUNCTION:
          int assignInitialPath()
DESCRIPTION:
          algorithm 9: the initial path of the algorithm is assigned and allocated.
          In words, we choose a path index l that is not currently in use (none of them are),
          and mark it as used. Then, for each layer lambda, we mark an index s such that
          both arrayPointer_P[lambda*L + s] and arrayPointer_C[lambda*L + s] are allocated to 
          the current path.
PARAMETERS:
          INPUT:
              None
	  OUTPUT:
	      index l of initial path
------------------------------------------------------------------------*/
int listDecoder::assignInitialPath(){
  int l;
  int m = (int)log2((double)N);
  l = inactivePathIndices[inactivePathIndices.size() - 1];
  inactivePathIndices.pop_back();
  activePath[l] = true;
  //associate arrays with path index
  for(int lambda = 0; lambda <= m; lambda++){
    int s = inactiveArrayIndices[lambda][inactiveArrayIndices[lambda].size() - 1];
    inactiveArrayIndices[lambda].pop_back();
    pathIndexToArrayIndex[lambda*L + l] = s;
    arrayReferenceCount[lambda*L + s] = 1;
  }
  return l;
}
/*------------------------------------------------------
FUNCTION:
        int clonePath(int l)
DESCRIPTION:
         algoritm 10: used to clone a path, the final step
before splitting the path in two. The logic is very similar
to that of algorithm 9, but now we make the two paths share 
bit-arrays and probability arrays.
PARAMETERS:
         INPUT:
	     l -- index of path to clone
         OUTPUT:
             l'--- index of copy
RETURN VALUES:
           l'
--------------------------------------------------------*/
int listDecoder:: clonePath(int l){
  int ll = inactivePathIndices[inactivePathIndices.size() -1];
  inactivePathIndices.pop_back();
  activePath[ll] = true;
  int m = (int)log2((double)N);
  //Make ll reference same array as l
  for(int lambda = 0; lambda <= m; lambda ++){
    int s = pathIndexToArrayIndex[lambda*L + l];
    pathIndexToArrayIndex[lambda*L + ll] = s;
    arrayReferenceCount[lambda*L + s]++;
  }
  return ll;
}
/*------------------------------------------------------
FUNCTION:
       void killPath(int l)
DESCRIPTION:
       Algorithm 11: terminate a path, which is achieved by
       marking it as inactive. After this is done, the arrays
       marked as associated with the path must be dealt with
       as follows. Since the path is inactive, we think of it
       as not having any associated arrays, and thus all the
       arrays that were previously associated with the path
       must have their reference count decreased by one.
PARAMETERS:
       INPUT: 
            l--index l of path to kill
       OUTPUT:
            all the arrays that were previously associated with
       the path must have their reference count decreased by one.
RETURN VALUES:
       None
  --------------------------------------------------------*/
void listDecoder::killPath(int l){
  //Mark the path index l as inactive
  activePath[l] = false;
  inactivePathIndices.push_back(l);
  //Disassociate arrays with path index
  int m = (int)log2((double)N);
  for(int lambda = 0; lambda <= m; lambda++){
    int s = pathIndexToArrayIndex[lambda*L + l];
    arrayReferenceCount[lambda*L + s]--;
    if(arrayReferenceCount[lambda*L + s] == 0)
      inactiveArrayIndices[lambda].push_back(s);
  }
}
/*-------------------------------------------
FUNCTION:
        vector<long double>* getArrayPointer_P(int lambda, int l);
DESCRIPTION:
         Algorithm 12:  access the probability-pair array associated
	 with a certain path l and layer lambda. There are two cases
         to consider: either the array is associated with more than 
         one path or it is not. If it is not, then nothing needs to
         be done, and we return a pointer to the array. On the other
         hand, if the array is shared, we make a private copy of the 
         path l, and return a pointer to that copy. By doing so, we 
         ensure that two paths will never write to the same array.
PARAMETERS:
         INPUT:
             lambda -- path
             l      -- layer
         OUTPUT:
             a pointer to vector<long double>
RETURN VALUES:
             a pointer to vector<long double>
  -------------------------------------------*/
vector<long double>* listDecoder::getArrayPointer_P(int lambda, int l){
  int ss;
  int s = pathIndexToArrayIndex[lambda*L + l];
  if(arrayReferenceCount[lambda*L + s] == 1){
    ss = s;
  }
  else{
    ss = inactiveArrayIndices[lambda][inactiveArrayIndices[lambda].size() - 1];
    inactiveArrayIndices[lambda].pop_back();
    //copy the contents of the array pointed to by arrayPointer_P[lambda*L + s]
    //into that pointed to by arrayPointer_P[lambda*L + ss]
    int m = (int)log2((double)N);
    for(int i = 0; i < (int)pow(2.0, m - lambda + 1); i++){
	 arrayPointer_P[lambda*L + ss][i] = arrayPointer_P[lambda*L + s][i];
	 arrayPointer_C[lambda*L + ss][i] = arrayPointer_C[lambda*L + s][i]; 
	} 
    arrayReferenceCount[lambda*L + s]--;
    arrayReferenceCount[lambda*L + ss] = 1;
    pathIndexToArrayIndex[lambda*L + l] = ss;
  }
  return &arrayPointer_P[lambda*L + ss];
}

vector<long double>* listDecoder::getArrayPointer_PCopy(int lambda, int l){
  int s = pathIndexToArrayIndex[lambda*L + l];
  return &arrayPointer_P[lambda*L + s];
}
/*----------------------------------------------
FUNCTION:
         vector<int> getArrayPointer_C(int lambda, int l)
DESCRIPTION:
           very similar to Algorithm 12
PARAMETERS:
          INPUT:
               lambda -- layer
               l      -- path index
          OUTPUT:
               a pointer to vector<int>
RETURN VALUES:
         a pointer to vector<int>
  ------------------------------------------------*/
vector<int>* listDecoder::getArrayPointer_C(int lambda, int l){
  int ss;
  int s = pathIndexToArrayIndex[lambda*L + l];
  if(arrayReferenceCount[lambda*L + s] == 1){
    ss = s;
  }
  else{
    ss = inactiveArrayIndices[lambda][inactiveArrayIndices[lambda].size() - 1];
    inactiveArrayIndices[lambda].pop_back();
    //copy the contents of the array pointed to by arrayPointer_P[lambda*L + s]
    //into that paointed to by arrayPointer_P[lambda*L + ss]
    int m = (int)log2((double)N);
    for(int i = 0; i < (int)pow(2.0, m - lambda + 1); i++){
	 arrayPointer_C[lambda*L + ss][i] = arrayPointer_C[lambda*L + s][i];
	 arrayPointer_P[lambda*L + ss][i] = arrayPointer_P[lambda*L + s][i];
	}
    arrayReferenceCount[lambda*L + s]--;
    arrayReferenceCount[lambda*L + ss] = 1;
    pathIndexToArrayIndex[lambda*L + l] = ss;
  }
  return &arrayPointer_C[lambda*L + ss];
}

vector<int>* listDecoder::getArrayPointer_CCopy(int lambda, int l){
  int s = pathIndexToArrayIndex[lambda*L + l];
  return &arrayPointer_C[lambda*L + s]; 
}
/*---------------------------------------------------
FUNCTION:
     bool pathIndexInactive(int l)
DESCRIPTION:
     Algorithm 13: simply a shorthand, meant to help readability later on
PARAMETERS:
     INPUT:
         l -- path index
     OUTPUT:
         true  -- if path l is active
         false -- otherwise
RETURN VALUES:
        See OUTPUT
  ----------------------------------------------------*/
bool listDecoder:: pathIndexInactive(int l){
  if(activePath[l])
    return false;
  else
    return true;
} 
/*----------------------------------------------------
FUNCTION:
         void recursivelyCalcP(int lambda, int varphi)
DESCRIPTION:
         Algorithm 14: list version 
PARAMETERS:
         INPUT:
              lambda -- layer
              varphi -- phase
         OUTPUT: 
  -----------------------------------------------------*/
void listDecoder::recursivelyCalcP(int lambda, int varphi){
  if(lambda == 0) return; //stopping condition
  int phi = varphi/2;
  int m = (int)log2((double)N);
  //recursive first, if needed
  if(varphi %2 == 0){ recursivelyCalcP(lambda - 1, phi); }
  //perform the calculation
  long double sigma = 0;
  for(int l = 0; l < L; l++){
    if(pathIndexInactive(l)){
      continue;
	}
    vector<long double>*Plambda   = getArrayPointer_P(lambda, l);
    vector<long double>*Plambda_1 = getArrayPointer_P(lambda - 1, l);
    vector<int>*Clambda           = getArrayPointer_C(lambda, l);
	
	
    for(int beta = 0; beta < (int)pow(2.0, m - lambda); beta++){
      if(varphi%2==0){ //apply equation 4
	//when the input is 0
	Plambda->at(beta) = 0.5*(Plambda_1->at(2*beta)*Plambda_1->at(2*beta + 1) + Plambda_1->at(2*beta + Plambda_1->size()/2)*Plambda_1->at(2*beta + 1 + Plambda_1->size()/2)); 
	
        //wen the input is 1
        Plambda->at(beta + Plambda->size()/2) =0.5*(Plambda_1->at(2*beta + Plambda_1->size()/2)*Plambda_1->at(2*beta + 1) + Plambda_1->at(2*beta)*Plambda_1->at(2*beta + 1 + Plambda_1->size()/2)); 
       
	long double max = (Plambda->at(beta + Plambda->size()/2) > Plambda->at(beta))? Plambda->at(beta+ Plambda->size()/2): Plambda->at(beta);
	sigma = (max > sigma)? max: sigma;
      }
      else{ //apply equation 5
	int u = Clambda->at(beta);
	if (u == 0){
	  //when the input is 0
	  Plambda->at(beta)                    = 0.5*(Plambda_1->at(2*beta)*Plambda_1->at(2*beta + 1));
	  //when the input is 1 
	  Plambda->at(beta + Plambda->size()/2) = 0.5*(Plambda_1->at(2*beta + Plambda_1->size()/2)*Plambda_1->at(2*beta + 1 + Plambda_1->size()/2));
	  }
	else{
	  //when the input is 0
	  Plambda->at(beta)                    = 0.5*(Plambda_1->at(2*beta + Plambda_1->size()/2)*Plambda_1->at(2*beta + 1));
	  //when the input is 1
	  Plambda->at(beta + Plambda->size()/2) = 0.5*(Plambda_1->at(2*beta)*Plambda_1->at(2*beta + 1 + Plambda_1->size()/2));
	}
	long double max = (Plambda->at(beta + Plambda->size()/2) > Plambda->at(beta))? Plambda->at(beta + Plambda->size()/2): Plambda->at(beta);
	sigma = (max > sigma)? max: sigma;
      }//end else
    }//end for (int beta...)
  }//end for (int l ...)
  
  //cout<<"after normalize probabilities..."<<endl;
  for(int l = 0; l < L; l++){
    if(pathIndexInactive(l)) continue;
    vector<long double> *Plambda = getArrayPointer_P(lambda, l);
    for(int beta = 0; beta < (int)pow(2.0, m - lambda); beta ++){
      Plambda->at(beta) = Plambda->at(beta)/sigma;
      Plambda->at(beta + Plambda->size()/2) = Plambda->at(beta + Plambda->size()/2)/sigma;
      //cout<<Plambda->at(beta)<<endl;
      //cout<<Plambda->at(beta + Plambda->size()/2)<<endl;
    }
  }
}
/*-----------------------------------------------------
FUNCTION:
      void recursivelyUpdateC(int lambda, int varphi)
DESCRIPTION:
      Algorithm 15: list version
PARAMETERS:
      INPUT:  
           lambda -- layer
           varphi -- phase
      OUTPUT:
  ------------------------------------------------------*/
void listDecoder::recursivelyUpdateC(int lambda, int varphi){
  if(varphi%2 == 0) return; //varphi must be odd 
  int m = (int)log2((double)N);
  int phi = varphi/2;
  for(int l = 0; l < L; l++){
    if(pathIndexInactive(l))continue; 
    vector<int>*Clambda   = getArrayPointer_C(lambda, l);
     vector<int>*Clambda_1 = getArrayPointer_C(lambda - 1, l);
    for(int beta = 0; beta < (int)pow(2.0, m - lambda); beta ++){
      if(phi%2==0){
	Clambda_1->at(2*beta)     = (Clambda->at(beta) + Clambda->at(beta + Clambda->size()/2))%2; //change for binary
	Clambda_1->at(2*beta + 1) = Clambda->at(beta   + Clambda->size()/2);
      }
      else{
	Clambda_1->at(2*beta     + Clambda_1->size()/2) = (Clambda->at(beta) + Clambda->at(beta + Clambda->size()/2))%2; //change for binary
	Clambda_1->at(2*beta + 1 + Clambda_1->size()/2) = Clambda->at(beta + Clambda->size()/2);
      }
    }//end for(int beta = 0; ...)
  }//end for(int l = 0; ....)

    if(phi%2 == 1)
      recursivelyUpdateC(lambda - 1, phi);
}

/*-------------------------------------------------
FUNCTION:
         void SCLDecoder(int *received, int *decoded, set<int>FrozenIndex, int *FrozenBit,double **W)
DESCRIPTION:
         Algorithm 16
PARAMETERS:
         INPUT:
              received    -- the received vector 
              decoded     -- the decoded  vector
	      FrozenIndex -- Frozen Set
	      FrozenBit   -- Frozen Bit
              W           -- channel information
	 OUTPUT:
              decoded
RETURN VALUES:
         None
  --------------------------------------------------*/
void listDecoder:: SCLDecoder(int *received, int *decoded, set<int>FrozenIndex, int *FrozenBit, double **W){
  //initialization
  int Findex = 0;
  initializeDataStructures();
  int l;
  l = assignInitialPath();
  vector<long double> *P0 = getArrayPointer_P(0, l);
  for(int beta = 0; beta < N; beta++){
    P0->at(beta)                = W[0][received[beta]];
    P0->at(beta + P0->size()/2) = W[1][received[beta]]; 
  }
     

  //main loop
  int m = (int)log2((double)N);
  for(int varphi = 0; varphi < N; varphi ++){
    //cout<<"Now, varphi is "<<varphi<<endl;
    recursivelyCalcP(m, varphi);
    if(FrozenIndex.find(varphi)!= FrozenIndex.end()){//frozen bit
      //cout<<"frozen bit "<<endl;
      continuePaths_FrozenBit(varphi, FrozenBit, Findex);
      Findex++;
    }   
    else{
      //cout<<"information bit"<<endl;
      continuePaths_UnfrozenBit(varphi);
	}
  if(varphi%2==1)
      recursivelyUpdateC(m, varphi);
  //cout<<"============="<<endl;
  }
  
  //return the best codeword in the list
  int templ = findMostProbablePathCRC(FrozenIndex);
  if(templ != -1)
  	l = templ;
  else
  	l = findMostProbablePath();
  vector<int>* C0 = getArrayPointer_C(0, l);
  int temp[N];
  for(int i = 0; i < N; i++)
	temp[i] = C0->at(i);
  encode(decoded, temp);
  /*cout<<"decoding result is "<<endl;
  for(int i = 0; i < N; i++)
	cout<<decoded[i]<<" ";
  cout<<endl;*/
}
/*----------------------------------------------------
FUNCTION:
       void continuePaths_FrozenBit(int varphi, int* FrozenBit, int Findex)
DESCRIPTION:
       Algorithm 17: it is the analog of line 6 in Algorithm 5 applied to all active paths
PARAMETERS:
       INPUT:
           varphi    -- phase
           FrozenBit -- FrozenBit information
           Findex    -- Frozen bit index
       OUTPUT:
           Update C
RETURN VALUES:
        None
  ----------------------------------------------------*/
void listDecoder:: continuePaths_FrozenBit(int varphi, int * FrozenBit, int Findex){
  int m = (int)log2((double)N);
  for(int l = 0; l < L; l++){
    if(pathIndexInactive(l)) continue;
    vector<int>*Cm = getArrayPointer_C(m, l);
    if(varphi%2 == 0)Cm->at(0) = FrozenBit[Findex];
    else  Cm->at(0 + Cm->size()/2) = FrozenBit[Findex];
  }
}
/*---------------------------------------------------
FUNCTION:
       void continuePaths_UnfrozenBit(int varphi)
DESCRIPTION:
      Algorithm 18: analog of parts of Algorithm 5. However,
      now, instead of choosing the most likely fork out of
      2 possible forks, we must typically choose the L most
      likely forks of 2L possible forks.
PARAMETERS:
      INPUT:
         varphi -- phase information
      OUTPUT:
         Complicated to describe the output...
	 Please refer to the paper
RETURN VALUES:
        None
  ----------------------------------------------------*/
void listDecoder::continuePaths_UnfrozenBit(int varphi){
  vector<long double> probForks(2*L);
  int i = 0;
  int m = (int)log2((double)N);
  //populate probForks
  for(int l = 0; l < L; l++){
    if(pathIndexInactive(l)){
      probForks[l]     = 0;
      probForks[l + L] = 0; 
    }
    else{
      vector<long double> *Pm = getArrayPointer_P(m, l);
      probForks[l]     = Pm->at(0);///(Pm->at(0) + Pm->at(Pm->size()/2))*log(channelInfo[varphi]) + pathEntrophy[l];
      probForks[l + L] = Pm->at(0 + Pm->size()/2);///(Pm->at(0) + Pm->at(Pm->size()/2))*log(channelInfo[varphi]) + pathEntrophy[l];
      i++;
    }
  }

  int rho = (2*i < L)? 2*i: L;
  vector<bool>contForks(2*L);
  //Populate contForks such that contForks[l][b] is true iff probForks[l][b] is one of the rho 
  //largest entries in probForks
  vector<long double> temp(probForks);
  QQuickSort(temp, 0, temp.size() - 1);
  long double threshold = temp[2*L - rho];
  //cout<<"threshold is "<<threshold<<endl;
 
  int index = 0;
  //first select only bigger one
  for(int l = L - 1; l >= 0; l--){
	/*cout<<"probForks["<<l<<"][0] is "<< probForks[l]<<endl;
	cout<<"probForks["<<l<<"][1] is "<< probForks[l + L]<<endl;
	cout<<"index is now "<<index<<endl;*/
    if(index < rho){
    	if(probForks[l] > threshold){
		//cout<<"set to true"<<endl;
      		contForks[l] = true;
		index ++;
	 }
    	else contForks[l] = false;
    }
    if(index < rho){
	//cout<<"contForks["<<l<<"][0] is "<< contForks[l]<<endl;
    	if(probForks[l + L] > threshold){
		//cout<<"set to true"<<endl;
		index ++;
      		contForks[l + L] = true;
	}
    	else contForks[l + L] = false;
	
	//cout<<"contForks["<<l<<"][1] is "<<contForks[l+L]<<endl;
    }
  }
   //then select equal one, until we obtain desired number
   for(int l = L - 1; l >= 0; l--){
	/*cout<<"probForks["<<l<<"][0] is "<< probForks[l]<<endl;
	cout<<"probForks["<<l<<"][1] is "<< probForks[l + L]<<endl;
	cout<<"index is now "<<index<<endl;*/
    if(index < rho){
    	if(probForks[l] == threshold){
		//cout<<"set to true"<<endl;
      		contForks[l] = true;
		index ++;
	 }
    }
    if(index < rho){
	//cout<<"contForks["<<l<<"][0] is "<< contForks[l]<<endl;
    	if(probForks[l + L] ==  threshold){
		//cout<<"set to true"<<endl;
		index ++;
      		contForks[l + L] = true;
	}
	//cout<<"contForks["<<l<<"][1] is "<<contForks[l+L]<<endl;
    }
  }

   /*for(int l = 0; l < L; l++){
	cout<<"contForks["<<l<<"][0] is "<< contForks[l]<<endl<<"contForks["<<l<<"][1] is "<<contForks[l+L]<<endl;
	cout<<"probForks["<<l<<"][0] is "<< probForks[l]<<endl<<"probForks["<<l<<"][1] is "<<probForks[l+L]<<endl;
	}*/
  //After the forks are marked, we first kill the paths for which both forks are discontinued.
  for(int l = 0; l < L; l++){
    if(pathIndexInactive(l)) {continue;}
    if(contForks[l] == false && contForks[l + L] == false){    
	killPath(l);
     }
  }
  //More details are here: continue paths for which one or both are the forks are marked.
  //In case of the latter, the path is first split, i.e., make a clone of it.
  for(int l = 0; l < L; l++){
    //cout<<"Now l is "<<l<<endl;
    //both forks are bad, or invalid
    if(contForks[l] == false && contForks[l + L] == false){
      //cout<<"both forks are bad, or invalid"<<endl;
	continue;
	}
    vector<int>*Cm = getArrayPointer_C(m,l);
    //both forks are good, we can not decide it yet...
    if(contForks[l] == true && contForks[l + L] == true){
      //cout<<"both forks are good, we can not decide it yet"<<endl;
      //pathEntrophy[l] = probForks[l];
      if(varphi%2 == 0)
	Cm->at(0) = 0;
      else
	Cm->at(0 + Cm->size()/2) = 0;
      int ll = clonePath(l);
      //pathEntrophy[ll] = probForks[ l + L];
      //cout<<"clone it to path "<<ll<<endl;
      vector<int>*TempCm = getArrayPointer_C(m, ll);
      if(varphi%2==0){
	TempCm->at(0) = 1;
	}
      else{
	TempCm->at(0 + TempCm->size()/2) = 1;
	}
    }
    
    //exactly one fork is good
    else{
      //cout<<"exactly one fork is good"<<endl;
      if(contForks[l]== true){
	//pathEntrophy[l] = probForks[l];
	if(varphi%2==0)
	  Cm->at(0) = 0;
        else
	  Cm->at(0 + Cm->size()/2) = 0;
	}
     else{
       //pathEntrophy[l] = probForks[l+L];
	if(varphi%2==0)
	  Cm->at(0) = 1;
        else
	  Cm->at(0+Cm->size()/2) = 1;
	}
    }
  }
  //print();
}
/*-------------------------------------------------
FUNCTION:
        int findMostProbablePath()
DESCRIPTION:
         Algorithm 19
PARAMETERS:
         INPUT:
	      None
	 OUTPUT:
	      the most probable path index
RETURN VALUES:
        The most probable path index
  -------------------------------------------------*/
int listDecoder::findMostProbablePath(){
  int ll = -1;
  long double pp = 0;
  int m = (int)log2((double)N);
 
  for(int l = 0; l < L ; l++){
    if(pathIndexInactive(l)){
      continue;
    }
    cout<<"path "<<l<<endl;
    vector<int> *Cm = getArrayPointer_C(m,l); 
    vector<long double> *Pm = getArrayPointer_P(m, l);
    vector<int>* C0 = getArrayPointer_C(0, l);
  	
    int temp[N], decoded[N];

    for(int i = 0; i < N; i++){
      temp[i] = C0->at(i);
    }
    encode(decoded, temp);
    cout<<"decoding result is "<<endl;
    for(int i = 0; i < N; i++)
      cout<<decoded[i]<<" ";
      cout<<endl;
    
  
      cout<<"The probability is "<<endl;
    if (Cm->at(0 + Cm->size()/2) == 1){
      cout<< Pm->at(0 + Pm->size()/2)<<endl;
      if(pp <  Pm->at(0 + Pm->size()/2)){	
	  ll = l;
	  pp = Pm->at(0 + Pm->size()/2);
	}//end of if pp..
    }//end of if Cm...
    else{
      cout<<Pm->at(0)<<endl;
      if(pp < Pm->at(0)){
	ll = l;
	pp = Pm->at(0);
      }//end of if pp...
    }//end of else

  }//end of for
  //cout<<"Most probable path is "<<ll<<endl;
    return ll;
}

/*-------------------------------------------------
FUNCTION:
        int findMostProbablePathCRC()
DESCRIPTION:
         Algorithm 19 with CRC
PARAMETERS:
         INPUT:
	      None
	 OUTPUT:
	      the most probable path index
RETURN VALUES:
        The most probable path index
  -------------------------------------------------*/
int listDecoder::findMostProbablePathCRC(const set<int>&FrozenIndex){
  int ll = -1;
  long double pp = 0;
  int m = (int)log2((double)N);

  
  for(int l = 0; l < L ; l++){
    if(pathIndexInactive(l))
      continue;
	
    vector<int>* C0 = getArrayPointer_C(0, l);
    int temp[N], decoded[N];
    for(int i = 0; i < N; i++)
	temp[i] = C0->at(i);
    encode(decoded, temp);

    int k = N - FrozenIndex.size();
    int  crc[k-16], crcBits[16];

    //obtain information bits for CRC
    for(int i = 0; i < k-16; i++)
	crc[i] = decoded[i];
    //obtain parity bits for crc16
    for(int i = 0; i < 16; i++)
	crcBits[i] = decoded[k-16+i];
	
    if(bin2int(crcBits, 16) == crc16(crc, k-16)){
	vector<int> *Cm = getArrayPointer_C(m,l); 
    	vector<long double> *Pm = getArrayPointer_P(m, l);
    	if (Cm->at(0 + Cm->size()/2) == 1){
	  if(pp < Pm->at(0 + Pm->size()/2)){	
	  	ll = l;
	  	pp = Pm->at(0 + Pm->size()/2);
		}//end of if pp..
    	}//end of if Cm...
    	else{
	  if(pp < Pm->at(0)){
		ll = l;
		pp = Pm->at(0);
      		}//end of if pp...
    	}//end of else
    }//end of if (bin...)
  
  }//end of for
	
    return ll;
}

