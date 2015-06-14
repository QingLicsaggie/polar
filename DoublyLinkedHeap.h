#ifndef DOUBLYLINKEDHEAP_H_
#define DOUBLYLINKEDHEAP_H_
#include<iostream>
#include<exception>
#include<math.h>
#include<vector>
#include<cstdlib>
#include<cstdio>
int u = 256;
using namespace std;

long double calC(long double a, long double b) { return (-(a+b)*log2l((a+b)/2) + a*log2l(a) + b*log2l(b)); }

long double calcDeltaI(long double a, long double b, long double a1, long double b1){return (calC(a,b) + calC(a1, b1) - calC(a + a1, b+ b1));}

class DoublyLinkedHeap; //class declaration

class EmptyDoublyLinkedHeapException{};

class Node{

private:
	long double a;       //W(y_i|0)
	long double b;      //W(\bar{y}_i|0)
	long double aa;     //W(y_{i+1}|0)
	long double bb;     //W(\bar{y}_{i+1}|0)
	long double deltaI; //contains the difference in capacity that would result from applying Lemma 7(how to construct polar codes) to y_i and y_{i+1}.
	Node *dLeft, *dRight; //dLeft is a pointer to the data element corresponding to the pair y_{i-1} and y_i or null if i = 1; likewise, dRight is a pointerto the element corresponding to the pair y_{i+1} and y_{i+2}.
	int h; //index in the heap
	friend class DoublyLinkedHeap;
public:

	Node(long double a_ = 0.0, long double b_ = 0.0, long double aa_= 0.0, long double bb_ = 0.0, long double delta = 0.0, Node *dl = NULL, Node *dr = NULL, int hh = 0): a(a_), b(b_), aa(aa_), bb(bb_), deltaI(delta), dLeft(dl), dRight(dr), h(hh){}
	~Node(){}
	long double geta() {return a; }
	long double getb() {return b; }
	long double getaa() {return aa;}	
	long double getbb() {return bb;}
	long double getDeltaI() {return deltaI;}
	int geth(){return h; }
	void seta(long double a1) { a = a1; }	
	void setb(long double b1) { b = b1; }
	void setaa(long double aa1) { aa = aa1; }
	void setbb(long double bb1) { bb = bb1;}
	void setdeltaI(long double delta) {deltaI = delta; }
	Node *getL() { return dLeft;}
	Node *getR() { return dRight;} 
};

class DoublyLinkedHeap{
protected:
	Node header, trailer;
private:
	//heap element
	int curSize; //number of elements in list
	bool orderOK; //true if array has heap-order property
	Node **array; //(dynamic) heap array
	int length; //the length of the array
	static const int DEF_SIZE = 1024*1024;
	void getNewArray(int newSize){
		array = new Node*[newSize];
		length = newSize;
	}
public: 
	void checkSize(); //double the heap array if it is full
	void toss(Node *x); //insert into array without maintaining heap order
	void walkDown(int hold);
	DoublyLinkedHeap(int size = DEF_SIZE): header(0), trailer(0){
	//default constructor	
		//heap part
		curSize = 0;
		orderOK = true;
		getNewArray(size);
		//doubly-link part
		header.dLeft  = &trailer;
		trailer.dRight = &header;
	}
	~DoublyLinkedHeap();
	Node *getFirst() const { return header.dRight; }
	Node *getLast() const { return trailer.dLeft; }
	Node *getNode(int i) const {return array[i];}
	int getcurSize() const {return curSize;}
	bool isEmpty() const {return curSize == 0;}
	void insertDLL(Node *newNode);	
	void removeDLL(Node *node);
	void buildHeap();
	void insertHeap(Node *newNode);
	void insert(Node *newNode){ insertDLL(newNode);insertHeap(newNode); }
	void setIndex();
	void removeMin() throw(EmptyDoublyLinkedHeapException);
	Node *getMin(){ return array[0]; }//return data with minimal deltaI 
	void print();
	void constructChannel(vector<long double> &channel);

};

void DoublyLinkedHeap::checkSize(){
	if(curSize == length){
		Node **oldArray = array;
		getNewArray(2*curSize);
		for(int i = 0; i < curSize; i++){
			array[i] = oldArray[i];
			delete oldArray[i];
		}
		delete [] oldArray;
	}
}

void DoublyLinkedHeap::toss(Node *x){
	//checkSize();
	array[curSize - 1] = x;
	/*if(array[(curSize-1)/2]->getDeltaI() >  x->getDeltaI())
		orderOK = false;*/
}

void DoublyLinkedHeap:: walkDown(int hole){
	int child;
	Node *key;
        key = array[hole];
	for(; 2*hole + 1 < curSize; hole = child){
		child = 2*hole + 1;
		if(child != curSize - 1 && array[child]->getDeltaI() > array[child + 1]->getDeltaI())
			child++; //right child = 2*hole + 2
		if(key->getDeltaI() > array[child]->getDeltaI()) {
			array[hole] = array[child];
			array[hole]->h = hole;
		}
		else break;
	}	
	array[hole] = key;
	array[hole]->h = hole;
	//record the element index in the heap
	//for(int i = location; i < curSize; i++)
	//	array[i]->h = i;
}

DoublyLinkedHeap::~DoublyLinkedHeap(){
	Node *prev_node, *node = header.dRight;
	while(node != &trailer){
		prev_node = node;
		node = node->dRight;
		delete prev_node;
	}
	delete [] array;
}


//insert into DLL based on the order of LR, LR = a/aa
void DoublyLinkedHeap::insertDLL(Node *newNode){
	if(isEmpty()){
		trailer.dLeft = newNode;
		header.dRight = newNode;
		newNode->dLeft = &header;
		newNode->dRight = &trailer;
		curSize++;
		return;
	}
  	//insert into the trail
	else {
		Node *current = getLast();
		current->dRight = newNode;
		newNode->dLeft  =  current;
		newNode->dRight = &trailer;
		trailer.dLeft   = newNode;	
		curSize++;
	}
}

void DoublyLinkedHeap::removeDLL(Node *node){
	if(isEmpty())
		throw EmptyDoublyLinkedHeapException();	
	node->dLeft->dRight = node->dRight;
	node->dRight->dLeft = node->dLeft;
}

void DoublyLinkedHeap:: buildHeap(){
	for(int i = (curSize - 2)/2; i >= 0; i--){
		walkDown(i);
	}
	orderOK = true;
}

void DoublyLinkedHeap::insertHeap(Node *node){
	int hole = curSize - 1;
	for(; hole > 0 && array[(hole - 1)/2]->getDeltaI() > node->getDeltaI(); hole = (hole -1)/2){
		array[hole] = array[(hole-1)/2];
	}
	array[hole] = node; 
}

void DoublyLinkedHeap::setIndex(){
	for(int i = 0; i < curSize; i++)
		array[i]->h = i;
}
void DoublyLinkedHeap::removeMin()throw(EmptyDoublyLinkedHeapException){
	removeDLL(array[0]);
	delete array[0];
	array[0] = array[curSize - 1]; 	
	curSize -- ;
	walkDown(0);		
}

void DoublyLinkedHeap::print(){
	Node *node;
	node = header.dRight;
	cout<<"linked list information is as follows:"<<endl;
	while(node != &trailer){
		cout<<node->geta()<<" "<<node->getb();
		cout<<"; "<<node->getaa()<<" "<<node->getbb();
		cout<<" "<<node->getDeltaI();
		cout<<endl;	
		node = node->dRight;	
	}
	/*cout<<"Heap information is as follows: "<<endl;
	for(int i = 0; i < curSize; i++){
		cout<<"current Node "<<i<<":"<<endl;
		cout<<array[i]->geta()<<" "<<array[i]->getb()<<" "<<array[i]->getaa()<<" "<<array[i]->getbb()<<" ";
		cout<<array[i]->getDeltaI()<<endl;
		cout<<"left node in the linked list "<<endl;
		cout<<array[i]->dLeft->geta()<<" "<<array[i]->dLeft->getb()<<" "<<array[i]->dLeft->getaa()<<" "<<array[i]->dLeft->getbb()<<" ";
		cout<<array[i]->dLeft->getDeltaI()<<endl;
		cout<<"Right node in the linked list "<<endl;
		cout<<array[i]->dRight->geta()<<" "<<array[i]->dRight->getb()<<" "<<array[i]->dRight->getaa()<<" "<<array[i]->dRight->getbb()<<" ";
		cout<<array[i]->dRight->getDeltaI()<<endl;
	}*/
}

void DoublyLinkedHeap:: constructChannel(vector<long double> &channel){
	vector<long double> p0, p1;
	Node *node;
	node = header.dRight;
	while(node != &trailer){
		p0.push_back(node->geta());
		p0.push_back(node->getb());
		p1.push_back(node->getb());
		p1.push_back(node->geta());
		node = node->dRight;	
	}
	
	node = trailer.dLeft;
	p0.push_back(node->getaa());
	p0.push_back(node->getbb());
	p1.push_back(node->getbb());
	p1.push_back(node->getaa());
	long double prob = 0;
	for(int i = 0; i < (int)p0.size(); i++){
		channel.push_back(p0[i]);
		//cout<<p0[i]<<" ";
		prob += p0[i];
	}
	//cout<<endl;
	
	for(int i = 0; i < (int)p1.size(); i++){
		channel.push_back(p1[i]);
		//cout<<p1[i]<<" ";
	}
	//cout<<endl;
	//cout<<"sum prob "<< prob<<endl;
}
#endif//DoublyLinkedHeap.h
