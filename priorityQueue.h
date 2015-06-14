#ifndef PRIORITYQUEUE.H
#define PRIORITYQUEUE.H
#include "binaryHeap.h"
#include<iostream>
using namespace std;
class PriorityQueue{
 private:
  BinaryHeap  T;
 public:
  void insertItem(const long double k,  const int &e){
   Item  it(k,e);
   T.insert(it);
 }

 Item& getElem(int i){return T.getElem(i);}

 void print(){
   int size = T.getSize();
   for(int i = 0; i < size; i++){
     cout<<"( "<< T.getElem(i).getElem()<<", " <<T.getElem(i).getKey()<<") ";
   }
   cout<<endl;
 }

 int search(const int  &e){ return T.search(e);}

 const int & minElement(){
   return T.findMin().getElem();
 }

 const int &maxElement(){
   return T.findMax().getElem();
 }

 const long double &maxKey(){
  return T.findMax().getKey();
}
 const int minKey(){
   return T.findMin().getKey();
 }
 
 int getSize(){return T.getSize();}

 void removeMin(){
   T.deleteMin();
 }

 void removeMax(){
   T.deleteMax();
 }

 void erase(int i){
   T.erase(i);
 }
};
#endif //priorityQueue.h
