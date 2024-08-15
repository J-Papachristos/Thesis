CPP=g++
CPP_FLAGS=-Wall -Wpedantic -O2

CUDA=nvcc
CUDA_FLAGS=-O2

all: sensingArray
main: sensingArrayCuda
sensingArray: sensingArray.cpp
	$(CPP) $(CPP_FLAGS) sensingArray.cpp -o sensingArray.exe
sensingArrayCuda: sensingArray.cu
	$(CUDA) $(CUDA_FLAGS) sensingArray.cu -o sensingArray.exe

clean:
	rm -f *.exe
	rm -f *.exp
	rm -f *.lib
	rm -f *.txt