CXX = mpic++
CXXFLAGS = -Wall

main: main.o BICGstab.o fonction.o operation.o param.o
	$(CXX) main.o BICGstab.o fonction.o operation.o param.o -o main

main.o: main.cpp BICGstab.h fonction.h operation.h param.h
	$(CXX) $(CXXFLAGS) -c main.cpp

BICGstab.o: BICGstab.cpp BICGstab.h fonction.h param.h operation.h
	$(CXX) $(CXXFLAGS) -c BICGstab.cpp

operation.o: operation.cpp operation.h fonction.h
	$(CXX) $(CXXFLAGS) -c operation.cpp

fonction.o: fonction.cpp fonction.h
	$(CXX) $(CXXFLAGS) -c fonction.cpp

param.o: param.cpp param.h
	$(CXX) $(CXXFLAGS) -c param.cpp

clean:
	rm -f *.o main