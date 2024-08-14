CPP=g++
CPP_FLAGS=-Wall -Wpedantic

CUDA=nvcc
CUDA_FLAGS=-O2

all: sensingArray
main: sensingArray
sensingArray: sensingArray.cu
	$(CUDA) $(CUDA_FLAGS) sensingArray.cu -o sensingArray.exe

clean:
	rm -f *.exe
	rm -f *.exp
	rm -f *.lib
	rm -f *.txt