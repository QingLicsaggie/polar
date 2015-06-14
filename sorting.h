#include<iostream>
#include<cstdio>
#include<cstdlib>
#include<time.h>
#include<vector>
#include<sys/time.h>
#include<sys/resource.h>
using namespace std;
int PPartition(vector<long double>&W, int start, int end);

void QQuickSort(vector<long double>&W, int start, int end);

long double selection(int m, vector<long double>sequence);
