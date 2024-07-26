CPP=g++
CPP_FLAGS=-Wall -Wpedantic

all: FiniteDifferences testerModule
main: FiniteDifferences
FiniteDifferences: FiniteDifferences.cpp
	$(CPP) $(CPP_FLAGS) FiniteDifferences.cpp -O2 -o FiniteDifferences.exe

testerModule:
	$(CPP) $(CPP_FLAGS) .\libs\FunctionDiscretizeTestModule.cpp -O2 -o .\libs\FunctionDiscretizeTestModule.exe

clean:
	rm -f *.exe
	rm -f libs\*.exe
	rm -f A_sparse.txt b.txt x.txt