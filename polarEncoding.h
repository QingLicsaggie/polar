#ifndef POLARENCODING_H
#define POLARENCODING_H
#include<set>
using namespace std;
//Function list
void polarEncode(int *codeword, int *information, int *frozenBit, set<int>frozenIndex);

void setInput(int *information, set<int>frozenIndex, int *frozenBit, int *Message);

#endif//POLARENCODING_H
