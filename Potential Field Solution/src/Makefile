CPP=g++
CPP_FLAGS=-Wall -Wpedantic

all: FiniteDifferences testerModule
main: FiniteDifferences
FiniteDifferences: FiniteDifferences.cpp
	$(CPP) $(CPP_FLAGS) FiniteDifferences.cpp -O2 -o FiniteDifferences.exe

testerModule:
	$(CPP) $(CPP_FLAGS) .\libs\FunctionDiscretizeTestModule.cpp -O2 -o .\libs\FunctionDiscretizeTestModule.exe

clean:
	erase .\*.exe
	erase .\libs\*.exe