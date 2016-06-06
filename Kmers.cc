#include "Kmers.h"
#include "util.h"

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
std::vector<jellyfish::mer_dna> vector_kmers_in_string(const std::string & str){
    auto k = jellyfish::mer_dna::k();
    std::vector<jellyfish::mer_dna> v(str.size()-k+1);
    for (size_t i = 0; i <= str.size() -k; i++) {
        v[i]=jellyfish::mer_dna(str.substr(i,k));
    }
    return v;
}
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
