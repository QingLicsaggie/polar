#ifndef LISTDECODING_H
#define LISTDECODING_H
#include<set>
#include<vector>
#include<cmath>
#include<iostream>
#include"channelPolarization.h"
#include"sorting.h"
#include"crc16.h"
extern int N;
using namespace std;
void recursivelyUpdateB(int lambda, int varphi, vector<int>&result);

void recursivelyUpdateC(int lambda, int varphi, vector<vector<int> > &C0, vector<vector<int> >&C1);

void recursivelyCalcP(int lambda, int varphi, vector<long double>&P0, vector<long double>&P1, vector<int>B);

void recursivelyCalcPSE(int lambda, int varphi, vector<vector<long double> >&P0, vector<vector<long double> >&P1, vector<vector<int> >&C0, vector<vector<int> >&C1);

void SCFirst(int *received, int *decoded, set<int>FrozenIndex, int *FrozenBit, double **W);

void SCFirstSE(int *received, int *decoded, set<int>FrozenIndex, int*FrozenBit, double**W);

void SCFirstSERevised(int *received1, int *received2, int *decode, set<int>FrozenIndex, int *FrozenBit, double **W1, double **W2);

class listDecoder{
  int L;
  vector<int>inactivePathIndices; 
  vector<bool>activePath;
  vector<vector<long double> > arrayPointer_P;
  vector<vector<int> > arrayPointer_C;
  vector<int>pathIndexToArrayIndex;
  vector<vector<int> > inactiveArrayIndices;
  vector<int> arrayReferenceCount;
 public:
 listDecoder(int l):L(l){};
 void initializeDataStructures();
 int assignInitialPath(); 
 int clonePath(int l);
 void killPath(int l);
 vector<long double>* getArrayPointer_P(int lambda, int l);
 vector<long double>* getArrayPointer_PCopy(int lambda, int l);
 vector<int>*         getArrayPointer_C(int lambda, int l);
 vector<int>*         getArrayPointer_CCopy(int lambda, int l);
 bool pathIndexInactive(int l);
 void recursivelyCalcP(int lambda, int varphi);
 void recursivelyUpdateC(int lambda, int varphi);
 void SCLDecoder(int *received, int *decoded, set<int>FrozenIndex, int *FrozenBit, double **W);
 void continuePaths_FrozenBit(int varphi, int *FrozenBit, int Findex);
 void continuePaths_UnfrozenBit(int varphi);
 int findMostProbablePath();
 int findMostProbablePathCRC(const set<int>&FrozenIndex);
 void print();
};
#endif//LISTDECODING_H
