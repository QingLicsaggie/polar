CC 	= g++ -std=c++0x        #The c++ compiler to use
CFLAGS  = -c -Wall -g -lm	        #The optimize code

all: polar

polar: generator.o channelPolarization.o polarEncoding.o noise.o polarDecoding.o sorting.o listdecoding.o analysis.o SC.o rewritting.o noisywem.o  crc16.o  main.o
	$(CC)  generator.o channelPolarization.o polarEncoding.o noise.o polarDecoding.o sorting.o  listdecoding.o analysis.o  SC.o rewritting.o noisywem.o  crc16.o  main.o  -o polar

generator.o: generator.cpp generator.h
	$(CC) $(CFLAGS) generator.cpp
 
channelPolarization.o: channelPolarization.cpp channelPolarization.h DoublyLinkedHeap.h
	$(CC) $(CFLAGS) channelPolarization.cpp

polarEncoding.o: generator.h generator.cpp channelPolarization.h channelPolarization.cpp polarEncoding.h polarEncoding.cpp
	$(CC) $(CFLAGS) polarEncoding.cpp

noise.o: noise.cpp noise.h
	$(CC) $(CFLAGS) noise.cpp

sorting.o: sorting.cpp sorting.h
	$(CC) $(CFLAGS) sorting.cpp

listdecoding.o: listdecoding.cpp listdecoding.h channelPolarization.cpp channelPolarization.h sorting.cpp sorting.h crc16.cpp crc16.h
	$(CC) $(CFLAGS) listdecoding.cpp

polarDecoding.o: generator.h generator.cpp channelPolarization.h channelPolarization.cpp polarEncoding.h polarEncoding.cpp  polarDecoding.h polarDecoding.cpp 
	$(CC) $(CFLAGS) polarDecoding.cpp

analysis.o: analysis.h analysis.cpp
	$(CC) $(CFLAGS) analysis.cpp

SC.o: SC.cpp SC.h analysis.cpp analysis.h channelPolarization.cpp generator.cpp generator.h channelPolarization.h polarEncoding.cpp polarEncoding.h polarDecoding.cpp polarDecoding.h listdecoding.cpp listdecoding.h
	$(CC) $(CFLAGS) SC.cpp

rewritting.o: rewritting.cpp rewritting.h analysis.cpp analysis.h channelPolarization.cpp channelPolarization.h SC.cpp SC.h generator.cpp generator.h listdecoding.cpp listdecoding.h
	$(CC) $(CFLAGS) rewritting.cpp

noisywem.o: noisywem.cpp noisywem.h channelPolarization.cpp channelPolarization.h
	$(CC) $(CFLAGS) noisywem.cpp


crc16.o: crc16.cpp crc16.h
	$(CC) $(CFLAGS) crc16.cpp

main.o: main.cpp generator.cpp generator.h channelPolarization.cpp channelPolarization.h DoublyLinkedHeap.h polarEncoding.h polarEncoding.cpp noise.cpp noise.h polarDecoding.cpp polarDecoding.h sorting.cpp sorting.h listdecoding.cpp listdecoding.h analysis.h analysis.cpp SC.h SC.cpp rewritting.cpp rewritting.h noisywem.h noisywem.cpp crc16.cpp crc16.h
	$(CC) $(CFLAGS) main.cpp

#
#clean up
#

clean:
	rm *.o 
	rm polar >/dev/null 2>&1
	

