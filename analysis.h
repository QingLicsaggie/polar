#ifndef ANALYSIS_H
#define ANALYSIS_H
#include<set>
#include<vector>
using namespace std;
double compute_I(double **W);

double errorBitRate(int *originalMessage, int *decodedMessage);

double blockErrorRate(int *originalMessage, int *decodedMessage);

void compute_average_performance();

void print_info(set<int> frozenset, double prob, double prob1);

void  concatePrint_info(set<int> frozenset, double prob, double prob1);
long double computeTheoreticalBER(set<int>frozenIndex, vector<long double>Probability);
#endif //ANALYSIS_H
