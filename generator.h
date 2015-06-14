#ifndef GENERATOR_H
#define GENERATOR_H
//Function list
void int2bin(int intstat, int *tempstat, int length);

long double bin2int(int *binseq, int length);

void generator(int **generator, int n);

void polarFMatrix(int **generator, int n);

void generatorC(int *generator, int j);

void generatorR(int *generator, int i);
#endif //GENERATOR_H



