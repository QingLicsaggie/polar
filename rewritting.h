#ifndef REWRITTING_H_
#define REWRITTING_H_
//function list
#include<set>
using namespace std;

void rewritter(int *current, int *dither, int *updated, int *data, set<int>FrozenIndex, double **W);

int decoder(int *received, int *dither, set<int>frozenIndex);

void findFrozenBits(int *codeword, int *result, set<int>FrozenIndex);

void experiment_for_rewriting(double **W, double **M, set<int>frozenIndex);

double  measureDistribution(int *codeword);

void flipToBalanced(int *codeword);
#endif//REWRITTING_H_
