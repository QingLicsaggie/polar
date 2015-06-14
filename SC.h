#ifndef SCENCODING_H
#define SCENCODING_H
#include<set>
using namespace std;
//function list
int compareBit(int c1, int c2);

long double measureDistortion(int *message, int *decodedcodeword, double **M);

long double sourceCodingVerify(double **W, double **M, set<int> frozenIndex);

void SCDecoding(int *u, int *decodedresult);

void SCEncoding(int *received, int *decoded, set<int>FrozenIndex, int*FrozenBit, double**W);

double compute_D(double **W, double **M);
#endif//SCENCODING_H
