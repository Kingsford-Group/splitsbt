#include "Kmers.h"
#include "util.h"
#include "BF.h"
/*
// return the number for each DNA base
int acgt(char c) {
    switch (c) {
    case 'A': case 'a': return 0;
    case 'C': case 'c': return 1;
    case 'G': case 'g': return 2;
    case 'T': case 't': return 3;
    default:
        DIE("Unknown nucleotide");
        return 0;
    }
}

// convert a kmer to a uint64_t
Kmer kmer_to_bits(const std::string & str) {
    Kmer b;
    for (int i = 0; i < str.size(); i++) {
        b = (b<<2) | acgt(str[i]);
    }
    return b;
}
*/

std::set<jellyfish::mer_dna> kmers_in_string(const std::string & str) {
    auto k = jellyfish::mer_dna::k();
    std::set<jellyfish::mer_dna> s;
    for (size_t i = 0; i <= str.size() - k; i++) {
        s.insert(jellyfish::mer_dna(str.substr(i, k)));
    }
    return s;
}
/*
std::list<jellyfish::mer_dna> list_kmers_in_string(const std::string & str){
    auto k = jellyfish::mer_dna::k();
    std::list<jellyfish::mer_dna> l;
    for (size_t i = 0; i <= str.size() -k; i++) {
        l.insert(jellyfish::mer_dna(str.substr(i, k)));
    }
    return l;
}
*/


// to exactly equal the SBT implementation as a set of kmers
std::vector<jellyfish::mer_dna> vector_kmers_in_string(const std::string & str){
    auto k = jellyfish::mer_dna::k();
    std::set<jellyfish::mer_dna> s;
    for (size_t i = 0; i <= str.size() - k; i++) {
        s.insert(jellyfish::mer_dna(str.substr(i, k)));
    }
    std::vector<jellyfish::mer_dna> v(s.size());
    std::copy(s.begin(), s.end(), v.begin());
    return v;
}

/*
std::vector<jellyfish::mer_dna> vector_kmers_in_string(const std::string & str){
    auto k = jellyfish::mer_dna::k();
    std::vector<jellyfish::mer_dna> v(str.size()-k+1);
    for (size_t i = 0; i <= str.size() -k; i++) {
        v[i]=jellyfish::mer_dna(str.substr(i,k));
    }
    return v;
}
*/
// Make sure npos is always right of opos
// Also make sure npos is inside array
void swap_kmer_position(std::vector<jellyfish::mer_dna> & v, int opos, int npos){
    if (npos > opos) {
        assert(npos < v.size());
        auto temp = v[npos];
        v[npos]=v[opos];
        v[opos]=temp;
    }
}

void swap_kmer_position(std::vector<size_t> & v, int opos, int npos){
    if (npos > opos) {
        assert(npos < v.size());
        auto temp = v[npos];
        v[npos]=v[opos];
        v[opos]=temp;
    }
}

// This assumes num_hashes is one and stores the direct size_t position.
// **More accurate to say it scales with num_hashes but is inefficient if num_hashes > 1
// The alternative would be to store two size_t values (base and inc) [more efficient than many size_t]
/*std::vector<size_t> vector_hash_conversion(BF* root, std::vector<jellyfish::mer_dna> v){
    HashPair hashes = root->get_hashes();
    unsigned long num_hash = root->get_num_hash();
    uint64_t size = root->size();

    std::vector<size_t> out(v.size()*num_hash);

    uint64_t index=0;
    for (size_t i=0; i<v.size(); i++){
        jellyfish::mer_dna m = v[i];
        m.canonicalize();
        uint64_t h0 = hashes.m1.times(m);
        uint64_t h1 = hashes.m2.times(m);
        size_t base = h0 % size;
        size_t inc = h1 % size;
        for (unsigned long j = 0; j < num_hash; ++j){
            const size_t pos = (base + j * inc) % size;
            out[index]=pos;
            index++;
        }        
    }
    if (index != out.size()){
        DIE("Index mismatch with output vector");
    }

    return out;
}*/
