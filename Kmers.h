#ifndef KMERS_H
#define KMERS_H
#include <set>
#include <string>

#include <jellyfish/mer_dna.hpp>
#include <jellyfish/rectangular_binary_matrix.hpp>

//int acgt(char c);
//Kmer kmer_to_bits(const std::string & str);
// set_kmers_in_string just sounds like a different operation
std::set<jellyfish::mer_dna> kmers_in_string(const std::string & str);
//std::list<jellyfish::mer_dna> list_kmers_in_string(const std::string & str);
std::vector<jellyfish::mer_dna> vector_kmers_in_string(const std::string & str);

void swap_kmer_position(std::vector<jellyfish::mer_dna> & v, int opos, int npos);
#endif
