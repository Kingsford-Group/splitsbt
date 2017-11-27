#ifndef MINHASH_H
#define MINHASH_H

#include <string>
#include <sys/mman.h>
#include <sdsl/bit_vectors.hpp>
#include <jellyfish/mer_dna_bloom_counter.hpp>
#include "Kmers.h"

using HashPair = jellyfish::hash_pair<jellyfish::mer_dna>;

int minhash(sdsl::bit_vector & bv, uint64_t seed, int hash);
int minhash_fast(sdsl::bit_vector & bv, uint64_t seed, int hash);
void sbthash(sdsl::bit_vector & bv, uint64_t* outarray);
uint64_t bitshift_hash(uint64_t bit_index, uint64_t seed, int hash);
#endif
