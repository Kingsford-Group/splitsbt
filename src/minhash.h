#ifndef MINHASH_H
#define MINHASH_H

#include <string>
#include <sys/mman.h>
#include <sdsl/bit_vectors.hpp>
#include <jellyfish/mer_dna_bloom_counter.hpp>
#include "Kmers.h"
#include "BloomTree.h"

using HashPair = jellyfish::hash_pair<jellyfish::mer_dna>;

class mh_node {
public:
    mh_node(std::string file, int numhash);
    //mh_node(mh_node left, mh_node right);
    //void write_instruct();

    std::string fname;
    std::vector<uint64_t> minhash;
    //mh_node* left, right;
};
//This MHcluster class is a reimplementation of the allsome SBT build
class MHcluster{
public:
    MHcluster(std::string minlist, int numhash, std::string outfile);
    ~MHcluster();
    void gcluster(const std::string &outfile, const std::string &hash_file);
    bool calcDMatrix();
    void unionAB2A_copyC2B(mh_node* a, mh_node* b, mh_node* c, std::string n);
    void unionAB2A(mh_node* a, mh_node* b, std::string n);
protected:
    std::string outfile;
    int nh;
    int num_elements;
    double** distance_matrix;
    std::vector<mh_node*> mhvector;
};


std::vector<uint64_t> minhash(sdsl::bit_vector & bv, uint64_t seed, int hash);
std::vector<uint64_t> minhash_fast(sdsl::bit_vector & bv, uint64_t seed, int hash);
void sbthash(sdsl::bit_vector & bv, uint64_t* outarray);
uint64_t bitshift_hash(uint64_t bit_index, uint64_t seed, int hash);
void write_minhash(std::vector<uint64_t> minhash, std::string bvfile1);
std::vector<uint64_t> read_minhash(std::string file, int numhash);
//int minhash_sim(std::vector<uint64_t> & mh1, std::vector<uint64_t> & mh2, int nh);
double minhash_sim(std::vector<uint64_t> & mh1, std::vector<uint64_t> & mh2, int nh);
void allsome_gclust(int nelements, double** distmatrix, std::string outfile, int nh);
double find_closest_pair(int n, double** distmatrix, int* ip, int* jp);
#endif
