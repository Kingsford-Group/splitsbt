#INCLUDES = -I ~/include
#LIBS = -L ~/lib 

CXXFLAGS= -std=c++11 -O3 $(shell pkg-config --cflags jellyfish-2.0) $(INCLUDES) $(LIBS)
LDFLAGS=-lsdsl -ldivsufsort -ldivsufsort64 $(shell pkg-config --libs jellyfish-2.0) -lpthread -lz $(shell pkg-config --libs-only-L jellyfish-2.0 | sed -e 's/-L/-Wl,-rpath,/')
#-DNDEBUG (add to CXXFLAGS to cancel assert lines)
# add -pg to CXXFLAGS / after $(LDFLAGS) below to restore gprof
# add -Wall to CXXFLAGS if you want to have warnings
bt: main.o Build.o Query.o Kmers.o BloomTree.o BF.o util.o Count.o SplitBloomTree.o minhash.o 
#SplitQuery.o
	$(CXX) -o $@ $^ $(LDFLAGS) 

clean:
	rm -f *.o bt
# DO NOT DELETE
