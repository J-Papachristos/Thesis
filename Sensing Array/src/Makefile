CPP=g++
CPP_FLAGS=-Wall -Wpedantic -O2

all: sensingArray
main: sensingArray
sensingArray: sensingArray.cpp
	$(CPP) $(CPP_FLAGS) sensingArray.cpp -o sensingArray.exe

clean:
	rm -f *.exe
	rm -f *.exp
	rm -f *.lib
	rm -f *.txt