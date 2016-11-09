CXXFLAGS= -std=c++11 -O3 $(shell pkg-config --cflags jellyfish-2.0)
LDFLAGS=-static -lsdsl -ldivsufsort -ldivsufsort64 $(shell pkg-config --libs jellyfish-2.0) -lpthread -lz
#-DNDEBUG (add to CXXFLAGS to cancel assert lines)
# add -pg to CXXFLAGS / after $(LDFLAGS) below to restore gprof
# add -Wall to CXXFLAGS if you want to have warnings
bt: main.o Build.o Query.o Kmers.o BloomTree.o BF.o util.o Count.o SplitBloomTree.o 
#SplitQuery.o
	$(CXX) -o $@ $^ $(LDFLAGS)

clean:
	rm -f *.o bt
# DO NOT DELETE
