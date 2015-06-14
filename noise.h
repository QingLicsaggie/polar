#ifndef NOISE_H
#define NOISE_H

//Function list
void BSC_corrupted(int *codeword,int*corruptedcodeword, double **W);

void mgrns(double mean, double sigma, int n, double * a);

void AWGN(double *original, double *noisy, double sigma, int length);

void AWGNModel(double sigma, int num);

void AWGNMatrix(double **matrix, double sigma);


int countFlip(int *codeword, int *noisycodeword);

void flipBits(int *codeword, int *corruptedcodeword, int FlipNum);

int quantize(double receive, double sigma, int length);

double quantizeProbBeing1(int num, int length);

void quantizeMatrix(int **matrix, int length);
#endif//NOISE_H
