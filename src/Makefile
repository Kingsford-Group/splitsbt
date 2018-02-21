INCLUDES = -I /usr/include/boost 
#LIBS = -L ~/lib 
#CXXFLAGS= -std=c++11 -O0 -g $(shell pkg-config --cflags htslib && echo "-DHAVE_HTSLIB") $(shell pkg-config --cflags jellyfish-2.0) $(shell pkg-config --cflags sdsl-lite) $(INCLUDES) $(LIBS)
CXXFLAGS= -std=c++11 -O3 $(shell pkg-config --cflags htslib && echo "-DHAVE_HTSLIB") $(shell pkg-config --cflags jellyfish-2.0) $(shell pkg-config --cflags sdsl-lite) $(INCLUDES) $(LIBS)
LDFLAGS=$(shell pkg-config --libs htslib) $(shell pkg-config --libs sdsl-lite) $(shell pkg-config --libs jellyfish-2.0) -lpthread -lz $(shell pkg-config --libs-only-L jellyfish-2.0 | sed -e 's/-L/-Wl,-rpath,/') $(shell pkg-config --libs-only-L htslib | sed -e 's/-L/-Wl,-rpath,/')
#-DNDEBUG (add to CXXFLAGS to cancel assert lines)
# add -pg to CXXFLAGS / after $(LDFLAGS) below to restore gprof
# add -Wall to CXXFLAGS if you want to have warnings
bt: main.o Build.o Query.o Kmers.o BloomTree.o BF.o util.o Count.o SplitBloomTree.o minhash.o 
#SplitQuery.o
	$(CXX) -o $@ $^ $(LDFLAGS) 

clean:
	rm -f *.o bt
# DO NOT DELETE
