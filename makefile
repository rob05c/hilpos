CUDA_CC=nvcc
CUDA_FLAGS=-O3 -I /usr/local/cuda/include -I ../../cuda/cub

LINK_CC=g++
LINK_FLAGS=-O3 -pthread -ldl -L/usr/local/cuda/lib64 -lcuda -lcudart #-ltbb

CC_CPP=g++
CPP_FLAGS=-O3 -std=c++03 -Wall -Werror -Wfatal-errors

CC_FLAGS=-O3 -std=c99 -Wall -Werror -Wfatal-errors
CC=gcc

all: hilpos
hilpos: main.o hilpos.o hilposcuda.o sqlite.o
	$(LINK_CC) main.o sqlite.o hilpos.o hilposcuda.o -o hilpos -lm $(LINK_FLAGS)
main.o: 
	$(CC_CPP) $(CPP_FLAGS) -c main.cpp -o main.o
hilpos.o: 
	$(CC_CPP) $(CPP_FLAGS) -c hilpos.cpp -o hilpos.o
sqlite.o: 
	$(CC) $(CC_FLAGS) -c sqlite3.c -o sqlite.o
hilposcuda.o: 
	$(CUDA_CC) $(CUDA_FLAGS) -c hilpos.cu -o hilposcuda.o
clean:
	rm -f *o hilpos
