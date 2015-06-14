#ifndef NOISYWEM_H_
#define NOISYWEM_H_
#include<set>
using namespace std;
//function list
class noisyWEM{
  int kForP;
  int kForW;
  int L;
  long double channelP;
  long double wemChannelP;
  long double error;
  long double cost;
  set<int>frozenIndex;
  set<int>frozenIndexForP;
  set<int>frozenIndexForW;
 public:
 noisyWEM(int p, int w, int l, long double cp, long double wp): kForP(p), kForW(w), L(l), channelP(cp), wemChannelP(wp), error(0), cost(0) {}; 
  int getL(){ return L;}
  void setFrozenIndex(const set<int>&frozen){frozenIndex = frozen;}
  void setFrozenIndexForP(const set<int>&frozenP){frozenIndexForP = frozenP;}
  void setFrozenIndexForW(const set<int>&frozenW){frozenIndexForW = frozenW;}
  void determineFrozenSet(char *fileforP, char *fileforW );
  void rewritterForNoisyWEM(int *current, int *dither, int *updated, int *data, double **WForWEM);
  void SCEncodingForNoisyWEM(int *received, int *result, int *frozenBit, double **WForWEM);
  void decoderForNoisyWEM(int *received, int *dither, int *result, double **WForP);
  int computeError(int *data, int *decoded);
  void experiment_of_noisyWEM(char *fileforP, char *fileforW, int trial, double **WForP, double **WForW, double **costM);
  void print_info_for_noisyWEM(double **WForW, double **WForP, double **costM, int count);
};
#endif//noisywem.h

