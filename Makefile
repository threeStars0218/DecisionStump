GPP = g++ -std=c++14
DSTUMPS = decisionstump.cpp decisionstump_1d.cpp

exec: test.cpp $(DSTUMPS)
	$(GPP) -o test test.cpp $(DSTUMPS)
