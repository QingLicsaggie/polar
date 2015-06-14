#ifndef CHANNELPOLARIZATION_H
#define CHANNELPOLARIZATION_H
#include<iostream>
#include<vector>
using namespace std;
//Function list
void encode(int *codeword, int *message);

long double recursiveChannelTransformations(int base, int index, int *y, int start, int end, int *u, int ubit, double **W);

void EquivalentChannel(const vector<long double> &original, vector<long double> &equival, vector<int> &map);

int Partition(vector<long double>&W, int start, int end);

void QuickSort(vector<long double>&W, int start, int end);

void bitChannelDegrading(const vector<long double>&W, vector<long double>&Q, int u, int i);

void testEquivalentChannel(int i, double **W);

void testdegrading_merge(int i, double **W, int u);

long double testBitChannelDegrading(double **W, vector<long double> &Q, int u, int i);

long double calculatePforBitChannel(double **W, int u, int i);

void computePforAllBitChannels(double **W, int u);

void DFSCalPforAllBitChannels(int u, double **W);

void DFSVisit(int start, int total, int constaint, vector<vector<long double> >&stack, int *color, double **W);

void fliplr(int *original, int*flip, int length);

int bitReverse(int num);

void selectChannels(vector<long double>P, set<int>&result, int k);

long double selectChannelP(vector<long double>P, set<int>&result, int k);

void selectChannellog(vector<long double>P, set<int>&result, int k);

void readVP(vector<long double>&P, char *filename);

void readFrozenset(set<int>&result);

double compute_elapsed( struct timeval* start, struct timeval* end);
#endif //CHANNELPOLARIZATION_H
