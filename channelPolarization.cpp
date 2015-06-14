#include<cstdio>
#include<cstdlib>
#include<ctime>
#include<algorithm>
#include<sys/time.h>
#include<sys/resource.h>
#include<math.h>
#include<iterator>
#include<iostream>
#include<fstream>
#include<sstream>
#include<vector>
#include<string>
#include<cstring>
#include<set>
#include"channelPolarization.h"
#include"generator.h"
#include"DoublyLinkedHeap.h"
using namespace std;
extern double theoretical_upper_bound;
extern int OUTPUTSIZE;
extern int Dimension;
extern int N;
extern int q;
extern int K;
/*--------------------------------------------
 FUNCTION:
             void printChannel(const vector<long double> &channel)
 DESCRIPTION:
             Print out the channel information in case we need it for debugging.
 PARAMETERS:
            INPUT:
                channel--channel information
            OUTPUT:
                 print the channel information
 RETURN VALUES:
            None
 LAST UPDATE DATES:
           05/20/2014
 ---------------------------------------------*/
void printChannel(const vector<long double> &channel){
  for(int i = 0; i < (int)channel.size()/2; i++)
        cout<<channel[i]<<" ";
    cout<<endl;
    for(int i = 0; i < (int)channel.size()/2; i++)
        cout<<channel[i+channel.size()/2]<<" ";
    cout<<endl;
}
/*----------------------------------------
FUNCTION:
	double recursiveChannelTransformations(int base, int index, int *y, int start, int end, int *u, int ubit, double **W)
DESCRIPTION:
	This function implements the Recursive Channel Transformation.
PARAMETERS:
	INPUT:
		base-2*N
		index-2*i - 1 or 2*i
		y- y[start, end]
		start-the starting point of y
		end-the ending point of y
		u-[0,...,index - 2]
		ubit - u[index-1]
		W-channel information
	OUTPUT:
		W_{2*N}^{2*i-1}(y_1^{2*N}, u_1^{2*i-2}|u_{2*i - 1})
                or
		W_{2*N}^{2*i}(y_1^{2*N}, U_1^{2*i- 2}|u_{2*i})
RETURN VALUES:
	See OUTPUT.
-----------------------------------------*/
long double recursiveChannelTransformations(int base, int index, int *y, int start, int end, int *u, int ubit, double **W){
	
	long double result = 0;
	int i, temp, *tempArrayOdd, *tempArrayEven, *tempArray, indexO, indexE;

	
	if (base == 2){//base case
		if (index %2 == 0){//even case
			for(i = 0; i < 2; i++){
				temp = (ubit + i)%2;
				result += 0.5*(W[temp][y[start]]*W[i][y[end]]); 
			}					
		}
		else{//odd case
			temp = (u[0] + ubit)%2;
			result = 0.5*(W[temp][y[start]]*W[ubit][y[end]]);
			//cout<<endl<<"----"<<endl;		
		}
	return result;
	}
	else{//recusive process begins
		if (index %2 == 0){//even case
			indexO = indexE = 0;
			tempArrayOdd  = new int[index/2]; 
			tempArrayEven = new int[index/2]; 
			tempArray     = new int[index/2]; 
			//copy tempArrayOdd and tempArrayEven
			for(i = 0; i <= index - 1; i ++){
				if(i%2 == 1){
					tempArrayOdd[indexO] = u[i];
					indexO++;
				}
				else{
					tempArrayEven[indexE] = u[i];
					indexE++;
				}			
			}
			//construct tempArray
			for(i = 0; i < index/2; i++)
				tempArray[i] = (tempArrayOdd[i] + tempArrayEven[i])%2;

			for(i = 0; i < 2; i++){
				temp = (ubit + i)%2;
				result += 0.5*(recursiveChannelTransformations(base/2, index/2, y, start, (start+end - 1)/2,   tempArrayOdd,temp, W)
					      *recursiveChannelTransformations(base/2, index/2, y, 1 + (start + end - 1)/2,  end, tempArray, i, W));
				}
		}
		else{//odd case
				indexO = indexE = 0;
				tempArrayOdd  = new int[(index - 1)/2]; 
				tempArrayEven = new int[(index - 1)/2]; 
				tempArray     = new int[(index - 1)/2]; 
				//copy tempArrayOdd and tempArrayEven
				for(i = 0; i <= index - 2; i++){
					if(i%2 == 1){
						tempArrayOdd[indexO] = u[i];
						indexO ++;
					}
					else{
						tempArrayEven[indexE] = u[i];
						indexE ++;
					}
				}
				//construct tempArray
				for(i = 0; i < (index - 1)/2 ; i++)
					tempArray[i] = (tempArrayOdd[i] + tempArrayEven[i])%2;

				result = 0.5*((recursiveChannelTransformations(base/2, (index-1)/2, y, start, (start+end -1 )/2, tempArrayOdd, (ubit+u[index-1])%2, W))
					    *recursiveChannelTransformations(base/2, (index-1)/2, y,(end+start - 1)/2 + 1, end, tempArray, ubit,  W));	
		}	
				delete[] tempArrayOdd;
				delete[] tempArrayEven;
				delete[] tempArray;				
				return result;						
	}
}

int ReverseShuffle(int original, int length){
  int tempInput[length];
  int bit[length];
  int2bin(original, bit, length);
  int2bin(original, tempInput, length);
  bit[0] = tempInput[length-1];
  for(int i = 0; i < length; i++)
    bit[i+1] = tempInput[i];
  return bin2int(bit, length);
}

void ReverseShuffleArray(int *input, int length){
  int temp[length];
  for(int i = 0; i < length; i++)
    temp[i] = input[i];
  for(int i = 0; i < length; i++){
    //cout<<i<<" -- > "<<ReverseShuffle(i, (int)log2(length))<<endl;
    input[ReverseShuffle(i, (int)log2(length))] = temp[i];
  }
  //cout<<endl;
}
/*--------------------------
bit inverse peration 
  ---------------------------*/
void bitInverse(int *input, int length){
  int tempInput[length];
  for(int i = 0; i < length; i++)
    tempInput[i] = input[i];

  int n = (int)log2(length);
  for(int i = 0; i < length; i++){
    int bit[n], flip[n];
    int2bin(i, bit, n);
    /*for(int j = 0; j < n; j++)
      cout<<bit[j]<<" ";
      cout<<endl;*/
    fliplr(bit, flip, n);
    /*for(int j = 0; j < n; j++)
      cout<<flip[j]<<" ";
      cout<<endl;*/
    int flipBit = (int)bin2int(flip, n);
    input[i] = tempInput[flipBit];
    //cout<<"output of "<<i<<" is input of "<<flipBit<<", which is "<<input[i]<<endl;
  }
}
void encodeEfficient(int *codeword, int start, int end, int *y){
	if(end - start + 1 == 2){
		y[start] = (codeword[0] + codeword[1])%2;//change this part for binary
		y[end]   = codeword[1];
		//cout<<"y["<<start<<"] is "<<y[start]<<endl;
		//cout<<"y["<<end<<"] is "<<y[end]<<endl;
		return;
	}
	else{
		int VTemp[end - start + 1];
		int half = (end - start + 1)/2;
		for(int i = 0; i < end - start + 1; i++){
		  if(i%2 == 0)
		    VTemp[i] = (codeword[i] + codeword[i+1])%2;
		  else
		    VTemp[i] = codeword[i];
		}
		/*cout<<"before reverse shuffle"<<endl;
		for(int i = 0; i < end - start + 1; i++)
		  cout<<VTemp[i]<<" ";
		  cout<<endl;*/
		ReverseShuffleArray(VTemp, end - start + 1);
		/*for(int i = 0; i < end - start + 1; i++)
		  cout<<VTemp[i]<<" ";
		  cout<<endl;*/
		int VFH[(end - start + 1)/2], VSH[(end - start + 1)/2];
		for(int i = 0; i < (end - start + 1)/2; i++){
			VFH[i] = VTemp[i];
			VSH[i]  = VTemp[i+half];
		}
		encodeEfficient(VFH, start, (start + end - 1)/2, y);
		encodeEfficient(VSH,  (start + end + 1)/2, end, y);
		
	}
}
void encode(int *codeword, int *message){
	encodeEfficient(message, 0, N - 1, codeword);
	//bitInverse(codeword);
	//cout<<"encoder result is "<<endl;
	/*for(int i = 0; i < N; i++)
	  cout<<codeword[i]<<" ";
	  cout<<endl;*/
}
/*--------------------------------------------------
 FUNCTION:
         void boxProduct(const vector<long double> &Q, vector<long double> &W)
 DESCRIPTION:
           box product of Q with itself.
 PARAMETERS:
          INPUT:
                Q--vector with probabilities
                W--store information for result
           OUTPUT:
                 W--
 RETURN VALUES:
             NONE
 COMMENTS:
              Now, just focus on binary input case
 LAST UPDATE DATE:
               05/20/2014
 ---------------------------------------------------*/
void boxProduct(const vector<long double> &Q, vector<long double> &W){
    
	int size1 = Q.size();
	int y1, y2;
	W.clear();
    
	//first W[y1,y2|0]
	for(y1 = 0; y1< size1/2; y1++)
		for(y2 = 0; y2 < size1/2; y2++)
			W.push_back(0.5*(Q[y1]*Q[y2] + Q[y1 + size1/2]*Q[y2 + size1/2]));
    
	//now W[y1,y2|1]
	for(y1 = 0; y1 < size1/2; y1++)
		for(y2 = 0; y2 < size1/2; y2++)
			W.push_back(0.5*(Q[y1 + size1/2]*Q[y2] + Q[y1]*Q[y2 + size1/2]));
}
/*---------------------------------------------------
 FUNCTION:
             void circleProduct(const vector<long double> &Q, vector<long double> &W)
 DESCRIPTION:
              circle Product Q with itself
 PARAMETERS:
            INPUT:
                 Q--vector with probabilities
                 W--store output result
             OUTPUT:
                 W--
 RETURN VALUES:
              NONE
 COMMENT:
             Now, focus on binary input channel case.
 LAST UPDATE DATE:
            05/20/2014
 ----------------------------------------------------*/
void circleProduct(const vector<long double> &Q, vector<long double> &W){
	
	int size1 = Q.size();
	int y1, y2, u1;
	W.clear();
	//first insert W[y1,y2, u1|0]
	for(y1 = 0; y1 < size1/2; y1++)
		for(y2 = 0; y2 < size1/2; y2++)
			for(u1 = 0; u1 < 2; u1++){
				if(u1 == 0){
                    W.push_back(0.5*Q[y1]*Q[y2]);
				}
				else{
                    W.push_back(0.5*Q[y1 + size1/2]*Q[y2]);
				}
			}
	
	//insert W[y1, y2, u1|1]
	for(y1 = 0; y1 < size1/2; y1++)
		for(y2 = 0; y2 < size1/2; y2++)
			for(u1 = 0; u1 < 2; u1++){
				if(u1==0){
					W.push_back(0.5*Q[y1 + size1/2]*Q[y2 + size1/2]);
				}
				else{
					W.push_back(0.5*Q[y1]*Q[y2 + size1/2]);
				}
			}
}
/*---------------------------------------------
 FUNCTION:
              int Partition(vector<long double>&W, int start, int end)
 DESCRIPTION:
              Partition the vector W[start, end] into two subvectors W[start, q -1] and W[q + 1, end] such that LR of each element of W[start, q -1] is less than to that of W[q].
 PARAMETERS:
             INPUT:
                 W--channel information
                 start--start point
                 end--end pont
              OUTPUT:
                 partition point
 RETURN VALUES:
               partition point
 COMMENT:
                this is the pre-step for quick sort for our channel.
 LAST UPDATE DATE:
             05/20/2014
 -----------------------------------------------*/
int Partition(vector<long double>&W, int start, int end){
	int i = start - 1;
	long double x = W[end]/W[W.size()/2 + end];
	for (int j = start; j <= end - 1; j++){
		if(W[j]/W[j + W.size()/2] <= x){
			i++;
			//exchange W[j] and W[i] and exchange W[j + W.size()/2] and W[i + W.size()/2]
			swap(W[j],W[i]);
			swap(W[W.size()/2 + j], W[i + W.size()/2]);
		}
	}
	//exchange W[i + 1] and W[end] and exchange W[ i + 1 + W.size()/2] and W[end + W.size()/2]
	swap(W[i + 1], W[end]);
	swap(W[i + 1 + W.size()/2], W[end + W.size()/2]);
	return i + 1;
}
/*----------------------------------------------
 FUNCTION:
             void QuickSort(vector<long double>&W, int start, int end)
 DESCRIPTION:
                  quickSort the channel information based on the LR, I do not
 think we can use current STL algorithm to solve this.
 PARAMETERS:
           INPUT:
                 W--channel information
                 start--start point
                 end--end point
           OUTPUT:
                W--increasing sorted result
 RETURN VALUES:
                None
 LAST UPDATE DATE:
                05/20/2014
 ------------------------------------------------*/
void QuickSort(vector<long double>&W, int start, int end){
	if (start < end){
		int q = Partition(W, start, end);
		QuickSort(W, start, q - 1);
		QuickSort(W, q + 1, end);
	}
}
/*-----------------------------------------------
 FUNCTION:
            void Merge(vector<long double>&W, int start, int middle, int end)
 DESCRIPTION:
               This algorithm is the final step of  merge sort the channel information based on the LR. The algorithm is added here
 to compare the time complexity of different algorithms.
 PARAMETERS:
             INPUT:
                    W     -- channel information
                     we do not elaborate the meaning of start, middle, and end, but instead, we explain the following terms
                    W[start, ...., middle] -- the first channel
                    W[middle+1, ..., end]  -- the second channel
              OUTPUT:
                    W -- increasing sorted result
 RETURN VALUES:
                  None
 LAST UPDATE DATE:
             05/20/2014
 ------------------------------------------------*/
void Merge(vector<long double>&W, int start, int middle, int end){
    int n1 = middle - start + 1;
    int n2 = end - middle;
    vector<long double> L(2*(n1 + 1));
    vector<long double> R(2*(n2 + 1));
    //copy left part of W to L
    for(int i = 0; i < n1; i++){
        L[i]               = W[start + i];
        L[i + n1 + 1]      = W[start + i  + W.size()/2];
    }
    //copy right part of W to R
    for(int j = 0; j < n2; j++){
        R[j]          = W[middle + j + 1];
        R[j + n2 + 1] = W[middle + j + 1 + W.size()/2];
    }
    int i = 0, j = 0;
    for(int k = start; k <= end; k ++){
        //decide which part to copy based on LLR information
        if( i < n1 && j < n2){
	  if ((log(L[i]) -log(L[i + n1 + 1])) <= (log(R[j])-log(R[j + n2 + 1]))){
                //copy the left
                W[k] = L[i];
                W[k + W.size()/2] = L[i + n1 + 1];
                i++;
            }
            else{
                //copy the right
                W[k] = R[j];
                W[k + W.size()/2] = R[j + n2 + 1];
                j++;
            }
        }
        else if(i == n1 && j < n2){
            //left part is done, and we copy right part
            W[k] = R[j];
            W[k + W.size()/2] = R[j + n2 + 1];
            j++;
        }
        else if(j == n2 && i < n1){
            //right part is done, and we copy left part
            W[k] = L[i];
            W[k + W.size()/2] = L[i + n1 + 1];
            i++;
        }
    }
}
/*----------------------------------------------
 FUNCTION:
             void MergeSort(vector<long double>&W, int start, int end)
 DESCRIPTION:xo
             Merge sort the channel information based on the LR
 PARAMETERS:
             INPUT:
                  W    -- channel information
                  start-- start point
                  end  -- end point
             OUTPUT:
                  W -- increasing sorted result
 RETURN VALUES:
               None
 LAST UPDATE:
             05/20/2014
 ------------------------------------------------*/
void MergeSort(vector<long double>&W, int start, int end){
    if(start < end){
        int middle = (start + end)/2;
        MergeSort(W, start, middle);
        MergeSort(W, middle + 1, end);
        Merge(W, start, middle, end);
    }
}
/*-----------------------------------------------
 FUNCTION:
          void Rearrange(vector<long double>&W)
 DESCRIPTION:
          Given the transition probability of the channel,
          we sort them mainly based on LR information,
          and when LR is 1, we further sort them based on its transition probability.
 PARAMETERS:
           INPUT:
               W--transition probability of the channel
           OUTPUT:
              sorted result of W
 RETURN VALUES:
            None
 LAST UPDATE:
          05/20/2014
 ------------------------------------------------*/
void Rearrange(vector<long double>&W){
    //cout<<"The size to sort is "<<W.size()<<endl;
    //We use merge sort here instead of quick sort as in quick sort, there is a random factor, which might effect time conplexity.
    MergeSort(W, 0, W.size()/2 - 1);
    //QuickSort(W, 0, W.size()/2 - 1);
    int start = 0, end=0;
    //find the starting point, where LR is 1
    for(int i = 0; i < (int)W.size()/2; i++){
        if(W[i] == W[i + W.size()/2]){
            start = i;
            break;
        }
    }
    //find the ending point of the section where LR is 1
    for(int i = start; i < (int)W.size()/2; i++){
        if (W[i] == W[i + W.size()/2])
            end = i;
    }
    /*time conplexity here can be further improved...*/
    if(start !=end){
        sort(W.begin() + start, W.begin() + end);
        vector<long double> copyW = W;
        for(int i = 0; i < (end - start + 1)/2; i++){
            W[start + i] = W[end - i] = copyW[start + 2*i];
            W[start + i + W.size()/2] =  W[end - i + W.size()/2 ] = copyW[start + 2*i];
        }
    }
}
/*----------------------------------------------
 FUNCTION:
             void degrading_merge(vector<long double>&W, vector<long double> &Q, int u)
 DESCRIPTION:
              This function degrades a binary-input channel to a channel with a prescribed output alphabet size(u).
 PARAMETERS:
             INPUT:
                  W-- a BMS channel W: X->Y, where |Y| = 2L,
                  u-- bound on the output alphabet size
                  Q-- the degraded channel, where its out size is bounded by u
              OUTPUT
                  Q
 RETURN VALUES:
              None
 COMMENT:
              Implement algorithm C of how to construct polar codes by Ido Tal and Alex Vardy.
 LAST UPDATE DATE:
              05/20/2014
 -----------------------------------------------*/
void degrading_merge(vector<long double>&W, vector<long double>&Q, int u){
    //if the size of the channel is already smaller than the desired size, nothing to do.
  if((int)W.size()/2 <=u) {
        Q = W;
	//cout<<"in side the function"<<endl;
        return;
    }
    else{
        //now we sort y_i based on LR value such that 1<=LR(y1)<LR(y2)<=...<=LR(yL).
        Rearrange(W);
	//cout<<"in degrading_merge"<<endl;
	//printChannel(W);
        int L = W.size()/4;
        int LL = W.size()/2;
        //not assure whether we need to change it for geneal (not binary output original channel) channel.
        DoublyLinkedHeap *dlh = new DoublyLinkedHeap();
        //insert part
        for(int i = L; i < LL - 1 ; i++){
            Node *d= new Node(W[i], W[LL + i], W[i+1], W[LL + i + 1], calcDeltaI(W[i], W[LL + i], W[i+1], W[LL + i + 1]));
            dlh->insert(d);
        }
        dlh->setIndex();
        //modify part
        int l = L;
        while( l > u/2){
            Node *d;
            Node *dLeft;
            Node *dRight;
            
            d = dlh->getMin();
            long double a1, b1;
            a1 = d->geta() + d->getaa();
            b1 = d->getb() + d->getbb();
			
            dLeft  = d->getL();
            dRight = d->getR();
            dlh->removeMin();
            l--;
            
            //cout<<"now l is "<<l<<endl;
            if (dLeft != NULL){
                dLeft->setaa(a1);
                dLeft->setbb(b1);
                dLeft->setdeltaI(calcDeltaI(dLeft->geta(), dLeft->getb(), a1,b1));
                dlh->walkDown(dLeft->geth());
            }
            if (dRight != NULL){
                dRight->seta(a1);
                dRight->setb(b1);
                dRight->setdeltaI(calcDeltaI(a1, b1, dRight->getaa(), dRight->getbb()));
                dlh->walkDown(dRight->geth());
            }
        }
        //construct channel
        dlh->constructChannel(Q);
        delete dlh;
    }
}
/*-------------------------------------
 FUNCTION:
           void scale(vector<long double> &degradecChannel)
 DESCRIPTION:
           Scale by the biggest element of degradedChannel.
 PARAMETERS:
         INPUT:
              degradedChannel-- the channel we want to scale
         OUTPUT:
              degraded Channel
 RETURN VALUES:
           None
 LAST UPDATE DATE:
          05/20/2014
 --------------------------------------*/
void scale(vector<long double> &degradedChannel){
    long double max = 0;
    long double min = 1;
    for(int i = 0; i < (int)degradedChannel.size(); i++){
        if(max < degradedChannel[i])
            max = degradedChannel[i];
	if(min > degradedChannel[i])
	   min = degradedChannel[i];
    }
    //cout<<"scale by "<<max<<endl;
    //cout<<"The smallest element is "<<min<<endl;
    for(int i = 0; i < (int)degradedChannel.size(); i++)
        degradedChannel[i] = degradedChannel[i]/max;
    for(int i = 0; i < (int)degradedChannel.size(); i++)
      if(degradedChannel[i]*pow(10.0,90)<1)
	degradedChannel[i] = (long double)1/pow(10.0,90);
}
/*---------------------------------------------
 FUNCTION:
         bitChannelDegrading(vector<long double>&W, vector<long double> &Q,  int u, int i)
 DESCRIPTION:
           degrading bit channel W_N^{(i)} with the constraint on the output alphabet size
 PARAMETERS:
           INPUT:
                W--original BMS channel
                Q--degraded BMS channel, where the result is storied
                u--bound on the output alphabet size
                i--bit channel index
           OUTPUT
                 Q--
 RETURN VALUES
              None
 LAST UPDATE DATES:
            05/20/2014
 ---------------------------------------------*/
void bitChannelDegrading(vector<long double>&W, vector<long double>&Q, int u, int i){
  vector<long double> tempQ, tempW(W);
    int n = (int)log2((float)N);
    int II[n],I[n];
    int2bin(i, II, n);
    fliplr(II, I, n);
    tempQ = tempW;
    
    for(int j = 0; j < n; j++){
        //cout<<j + 1<<" -th bit "<<endl;
        vector<long double> tempWW, tempQQ;
        if(I[j] == 0){
	  //cout<<"box product"<<endl;
            boxProduct(tempQ, tempWW);
        }
        else{
	  //cout<<"circle product"<<endl;
            circleProduct(tempQ, tempWW);
        }
	//printChannel(tempWW);
        degrading_merge(tempWW, tempQQ, u);
        scale(tempQQ);
	//printChannel(tempQQ);
        tempQ = tempQQ;
        Q = tempQQ;
    }
}
/*------------------------------------------
 FUNCTION:
             long double CalcP(vector<long double>Q)
 DESCRIPTION: 
             Given the degraded channel, we calculate the
 theoretical upper bound of error correction ability.
 That is, for each output, we choose the smallest one
 correponds to 0 or 1, and then sums that up.
 PARAMETERS:
              INPUT:
                  Q -- the channel
              OUTPUT:
                block  error rate.
 RETURN VALUES:
              block error rate.
 LAST UPDATE DATE:
 05/20/2014
 *-----------------------------------------*/
long double CalcP(vector<long double>Q){
    long double prob = 0;
    long double sum = 0;
    //since we have scaled, we calculate the sum for rescale in case we need it.
    for(int i = 0; i < (int)Q.size()/2; i++){
        sum += Q[i];
        if (Q[i] < Q[Q.size()/2 + i])
            prob += 0.5*Q[i];
        else
            prob += 0.5 *Q[Q.size()/2 + i];
    }
    return prob/sum;
}
/*--------------------------------------------
 FUNCTION:
             long double calculatePforBitChannel(double **W, int u, int i)(double **W, int u, int i)
 DESCRIPTION:
                calculate error probability of maximal likelyhood decoding for bit channel i based on degraded merge
 PARAMETERS:
            INPUT:
                  W--channel information
                  u--constraint on the output alphabet
                  i--index for bit channel
            OUTPUT
                  error probability of maximal likelyhood decoding
 RETURN VALUE
             error probability
 LAST UPDATE DATE:
          05/20/2014, update it to general binary input OUTPUTSIZE output symmetric channel
 --------------------------------------------*/
long double calculatePforBitChannel(double **W, int u, int i){
    vector<long double>WW, Q;
    for(int j = 0; j < OUTPUTSIZE; j++){//changed here for general case
        WW.push_back(W[0][j]);
	//cout<<W[0][j]<<" ";
    }
    //cout<<endl;
    for(int j = 0; j < OUTPUTSIZE; j++){//changed here for general case
        WW.push_back(W[1][j]);
	//cout<<W[1][j]<<" ";
    }
    //cout<<endl;
    bitChannelDegrading(WW, Q, u, i);
    return CalcP(Q);
}
/*----------------------------------------
 FUNCTION:
            void computePforAllBitChannels(double **W, int u)
 DESCRIPTION:
             compute error probabilities for all bit channels
 PARAMETERS:
            INTPUT:
                   W--channel information
                   u--constraint on the output size
             OUTPUT:
                  store error probabilities for all bit channels in a file
 RETURN VALUES:
              None
 LAST UPDATE DATES:
              05/20/2014
 ------------------------------------------*/
void computePforAllBitChannels(double **W, int u){
    stringstream ss;
    ss<<N;
    string nstr = ss.str();
    string name = nstr + ".txt";
    char *filename = new char[30];
    strcpy(filename, name.c_str());
    //ofstream out(filename);
    ofstream out;
    out.open(filename, ios::app);
    
    for(int i = 0; i < N; i++){
        //if (i%100==0)
        cout<<"in process of "<<i<<endl;
        long double value = calculatePforBitChannel(W, u, i);
        out<<value<<endl;
        cout<<value<<endl;
        cout<<endl<<"---"<<endl;
    }
    out.close();
    delete[] filename;
}
/*--------------------------------
 FUNCTION:
             DFSCalPforAllBitChannels(int u, double **W)
 DESCRIPTION:
             Calculate Probabilities for all bit channels
 based on DFS. This can greatly decrease t ime
 conplexity.
 PARAMETERS:
              INPUT:
                    u--constraint on the output size of channel
                    W--original channel
              OUTPUT:
                    Degraded channels
 RETURN VALUES:
               None
 LAST UPDATE:
               05/20/2014
 ---------------------------------*/
void DFSCalPforAllBitChannels(int u, double **W){
    int color[2*N - 1];
    for(int i = 0; i < 2*N -1; i++) color[i] = 0;
    //follows the convention of DFS, we create a color array.
    vector<vector<long double> > stack;
    //file where the result is stored
    stringstream ss;
    ss<<N;
    ss<<W[0][0];
    string nstr = ss.str();
    string name = nstr + ".txt";
    char *filename = new char[30];
    strcpy(filename, name.c_str());
    ofstream out(filename);
    DFSVisit(0, 2*N - 1, u, stack, color, W);
    delete[] filename;
}
/*--------------------------------
 FUNCTION:
              PreDFSVisit(int start, int total, int constaint, vector<vector<long double> >&stack, int *color, double **W)
 DESCRIPTION:
               We can start from internal node instead of the 0 node.
               This file is to do some preparaton work.
 COMMENT:
               Still needs time to fix it...
 LAST UPDATE DATES:
                05/20/2014
 ---------------------------------*/
void PreDFSVisit(int start, int total, int constaint, vector<vector<long double> >&stack, int *color, double **W){
    vector<long double>WW;
    for(int i = 0; i < OUTPUTSIZE; i++)
        WW.push_back(W[0][i]);
    for(int i = 0; i < OUTPUTSIZE; i++)
        WW.push_back(W[1][i]);
    stack.push_back(WW);
	
    int n = (int)log2((float)N);
    int II[n],I[n];
    int2bin(start, II, n);
    fliplr(II, I, n);
    for(int i = 0; i < n - 1; i++){
        vector<long double> parent = stack[stack.size() - 1];
        vector<long double> tempChannel, degradedChannel;
        if(I[i] == 0){
            boxProduct(parent, tempChannel);
        }
        else{
            circleProduct(parent, tempChannel);
        }
        degrading_merge(tempChannel, degradedChannel, constaint);
        stack.push_back(degradedChannel);
    }
    DFSVisit(start + N -1, total, constaint, stack, color, W);
}

/*-------------------------------------
 FUNCTION:
             void Rescale(vector<long double> &degradecChannel)
 DESCRIPTION:
             Scale back.
 PARAMETERS:
             INPUT:
               degradedChannel--already degraded channel
              OUTPUT:
            Rescaled channel.
 RETURN VALUES:
              None
 LAST UPDATE DATE:
               05/20/2014
 --------------------------------------*/
void Rescale(vector<long double> &degradedChannel){
    long double sum = 0;
    for(int i = 0; i < (int)degradedChannel.size()/2; i++){
        sum += degradedChannel[i];
    }
    for(int i = 0; i < (int)degradedChannel.size(); i++)
        degradedChannel[i] = degradedChannel[i]/sum;
}
/*-------------------------------------
 FUNCTION:
                 DFSVisit(int start, int total, int constaint, vector<vector<long double> >&stack, int *color, double **W)
 DESCRIPTION:
                Details of DFS implementation, and we just follow the outline provided in the book of "introduction to algorithms"
 PARAMETERS:
             INPUT:
                 start      --start index of the tree
                 total      --total nodes number
                 constraint --constraint size
                 stack      --we store some results to save time
                 color      --array indicating whether a node is done or not
                  W          --the original binary input channel, and it is supposed to symmetric.
              OUTPUT:
                  The result is stored in a file
 RETURN VALUES:
                   None
 LAST UPDATE DATE:
              05/20/2014
 --------------------------------------*/
void DFSVisit(int start, int total, int constaint, vector<vector<long double> >&stack, int *color, double **W){
    color[start] = 1; //we are working on it now...
    if (start == 0){
        vector<long double>WW;
        for(int i = 0; i < OUTPUTSIZE; i++){
            WW.push_back(W[0][i]);
	    //cout<<W[0][i]<<" ";
	}
	//cout<<endl;
        for(int i = 0; i < OUTPUTSIZE; i++){
            WW.push_back(W[1][i]);
	    //cout<<W[1][i]<<" ";
	}
	//cout<<endl;
        stack.push_back(WW);
    }
    else{
        vector<long double> parent = stack[stack.size() - 1];
        vector<long double> tempChannel, degradedChannel;
        if(start%2 == 1){
            boxProduct(parent, tempChannel);
        }
        else{
            circleProduct(parent, tempChannel);
        }
        degrading_merge(tempChannel, degradedChannel, constaint);
        scale(degradedChannel);
        //printChannel(degradedChannel);
        stack.push_back(degradedChannel);
    }
    if(start*2 + 1< total){//has child
        for(int children = 0; children < 2; children ++){
            if(color[2*start + children + 1] == 0)
                DFSVisit(2*start + children + 1, total, constaint, stack, color, W);
        }
    }
    else{ // no child, the desired result...
        cout<<"Finish node "<<start - (N -1)<<endl;
        //file
        stringstream ss;
        ss<<N;
	ss<<W[0][0];
        string nstr = ss.str();
        string name = nstr + ".txt";
        char *filename = new char[30];
        strcpy(filename, name.c_str());
        ofstream out;
        out.open(filename, ios::app);
        
        //we contruncate here
        vector<long double> channel = stack[stack.size() - 1];
        //Rescale(channel);
        long double prob = CalcP(channel);
        //printChannel(channel);
        cout<<prob<<"          "<<endl;
        /*if(prob > 0.45) {
            // printChannel(channel);
            cout<<"may not accurate, turn to traditional methods for channel "<<start- (N-1)<<endl;
            prob = calculatePforBitChannel(W,  constaint, start - (N-1));
            cout<<prob<<endl;
	    }*/
        out<<prob<<endl;
        out.close();
        delete[] filename;
    }
    stack.pop_back();
    color[start] = 2;//Finish
}
/*----------------------------------------------
 FUNCTION:
             long double testBitChannelDegrading(double **W, vector<long double> &Q, int u, int i)
 DESCRIPTION:
            test function for bitChannelDegrading()
 PARAMETERS:
           INPUT:
               W--the channel we want to test
               Q--the place we store results
               u--constraint on the size
                i--index of the channel
           OUTPUT:
               block error rate
 RETURN VALUES:
            NONE
 LAST UPDATE DATE:
          05/20/2014
 -----------------------------------------------*/
long double testBitChannelDegrading(double **W, vector<long double> &Q, int u, int i)
{
    vector<long double>WW;
    for(int j = 0; j < OUTPUTSIZE; j++)
        WW.push_back(W[0][i]);
    for(int j = 0; j < OUTPUTSIZE; j++)
        WW.push_back(W[1][i]);
    bitChannelDegrading(WW, Q, u, i);
    
    long double prob = 0;
    for(int j = 0; j < (int)Q.size()/2; j++){
        if (Q[j] < Q[Q.size()/2 + j]){
            prob += 0.5*Q[j];
        }
        else {
            prob += 0.5 *Q[Q.size()/2 + j];
        }
    }
    return prob;
}
/*-----------------------------------------------
 FUNCTION:
             void fliplr(int *original, int *flip, int length)
 DESCRIPTION:
              This function flips the given array original, and the result is stored in flip.
 PARAMETERS:
           INPUT:
               original-original array
               flip-resulting array
                length-length of original arraynd resulting array
           OUTPUT:
               flip
 RETURN VALUES:
           None
 LAST UPDATE DATE:
           05/20/2014
 ------------------------------------------------*/
void fliplr(int *original, int*flip, int length){
    int i;
    for(i=0; i < length; i++)
        flip[i]= original[length - 1 -i];
}
/*------------------------------------------------
 FUNCTION:
           void selectChannels(vector<long double> P, vector<int> result, int k)
 DESCRIPTION:
           This function selects N - K largest channels.
 PARAMETERS:
           INPUT:
               P--P parameters
               result--
               k --
           OUTPUT:
                result
 RETURN VALUES:
            None.
 LAST UPDATE DATES:
            05/20/2014
 --------------------------------------------------*/
void selectChannels(vector<long double>P, set<int>&result, int k){

    long double max = 0.0;
    int  location;
    theoretical_upper_bound = 0.0;
    vector<long double> CopyP(P);
    for(int i = 0; i < N - k; i++){
        max = -100000;
        for(int j = 0; j < N; j++ ){
	  //cout<<CopyP[j]<<endl;
            if (CopyP[j] > max && CopyP[j] !=0){
                max = CopyP[j];
                location = j;
            }
        }
        CopyP[location] = 0;
        result.insert(location);
    }
    
    for(int i = 0; i < N; i++)
        if(CopyP[i]!=0){
            theoretical_upper_bound += CopyP[i];//pow(2.0, CopyP[i]);
            //cout<<i<<" ---"<<CopyP[i]<<endl;
        }
    cout<<"Theoretical upper bound is "<<theoretical_upper_bound<<endl;
}
/*-------------------------------------
 FUNCTION:
               void selectChannellog(vector<long double>P, set<int>&result, int k)
 DESCRIPTION:
               The degraded channel is stored in log form, and we use it to select channel
 PARAMETERS:
              INPUT:
                  P--channel information
                  result--
                  k--how many channels to select
 OUTPUT
                 results is stored in the set
 RETURN VALUES:
                NONE
 LAST UPDATE DATE:
              05/20/2014
 --------------------------------------*/
void selectChannellog(vector<long double>P, set<int>&result, int k){
    long double max = 0.0;
    int  location;
    theoretical_upper_bound = 0.0;
    vector<long double> CopyP(P);
    for(int i = 0; i < N - k; i++){
        max = -100000;
        for(int j = 0; j < N; j++ ){
            if (CopyP[j] > max && CopyP[j] !=0){
                max = CopyP[j];
                location = j;
            }
        }
        //cout<<location<<" "<<max<<endl;
        CopyP[location] = 0;
        result.insert(location);
    }
    for(int i = 0; i < N; i++)
        if(CopyP[i]!=0)
            theoretical_upper_bound += pow(10.0, CopyP[i]);
    cout<<"Theoretical upper bound is "<<theoretical_upper_bound<<endl;
}
/*--------------------------------------------------
 FUNCTION:
               void readVP(vector<long double>&P, char *filename)
 DESCRIPTION:
               Read P from filename
 PARAMETERS:
               INPUT:
                     P--
                     filename--
                     OUTPUT:
                        filename
 RETURN VALUES:
              None 
 LAST UPDATE DATE
              05/20/2014
 ---------------------------------------------------*/
void readVP(vector<long double>&P, char*filename){
    ifstream file(filename);
    string line;
    while(getline(file, line)){
        P.push_back((long double)atof(line.c_str()));	
    }
}


void readFrozenset(set<int>&result){	
    ifstream file("frozenIndex.txt");
    string line;
    while(getline(file, line)){
        result.insert((long double)atof(line.c_str()));	
    }
}

int bitReverse(int num){
  int n = (int)log2(N);
  int bit[n];
  int2bin(num, bit, n);
  int flip[n];
  fliplr(bit, flip, n);
  return (int)bin2int(flip, n);
}
