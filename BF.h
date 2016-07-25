#ifndef BF_H
#define BF_H

#include <string>
#include <sys/mman.h>
#include <sdsl/bit_vectors.hpp>
#include <jellyfish/mer_dna_bloom_counter.hpp>
#include "Kmers.h"

using HashPair = jellyfish::hash_pair<jellyfish::mer_dna>;

// a kmer bloom filter
class BF {
public:
    BF(const std::string & filename, HashPair hp, int nh);
    BF(const std::string & filename, BF* copy);
    //Added overloaded method to store second filename
    //BF(const std::string & f1, const std::string & f2, HashPair hp, int nh);
    virtual ~BF();

    virtual void load();
    virtual void save();

    virtual int operator[](uint64_t pos) const;
    virtual void set_bit(uint64_t p);
    virtual uint64_t size() const;

    virtual std::string get_name() const;
    virtual HashPair get_hashes() const;
    virtual int get_num_hash() const;

    virtual bool contains(const jellyfish::mer_dna & m) const;
    virtual bool contains(const size_t pos) const;
    bool contains(const std::string & str) const;

    virtual bool contains(const jellyfish::mer_dna & m, int type) const;
    virtual bool contains(const size_t pos, int type) const;
    bool contains(const std::string & str, int type) const;

    void add(const jellyfish::mer_dna & m);

    virtual uint64_t similarity(const BF* other, int type) const;
    virtual std::tuple<uint64_t, uint64_t> b_similarity(const BF* other) const;
    virtual BF* union_with(const std::string & new_name, const BF* f2) const;
    virtual void union_into(const BF* f2);
    virtual void union_into(const BF* f2, int type);
    virtual uint64_t count_ones() const;
    virtual uint64_t count_ones(int type) const;
    virtual void compress();
    virtual void compress(BF* rm);

    virtual BF* sim_with(const std::string & new_name, const BF* f2) const;
    virtual void sim_into(const BF* f2);

    virtual BF* dif_with(const std::string & new_name, const BF* f2) const;
    virtual void dif_into(const BF* f2);
    
    virtual void update_mask(const BF* update);
    virtual void update_mask(const BF* u1, const BF* u2);
    // Things which are not used in all cases but were added to save time in interface
protected:
    std::string filename;
    sdsl::rrr_vector<255>* bits;

    HashPair hashes;
    unsigned long num_hash;
};

class compressedSBF : public BF {
public:
    compressedSBF(const std::string & filename, HashPair hp, int nh, uint64_t size =0);
    virtual ~compressedSBF();

    virtual void load();
    virtual void save();

    virtual std::string get_sim_name();
    virtual std::string get_dif_name();

    virtual uint64_t size() const;

    virtual int operator[](uint64_t pos) const;

    virtual bool contains(const jellyfish::mer_dna & m, int type) const;
    virtual bool contains(const size_t pos, int type) const;
    bool contains(const std::string & str, int type) const;

protected:
    sdsl::rrr_vector<255>* sim_bits;
    sdsl::rrr_vector<255>* dif_bits;
};

class UncompressedBF : public BF {
friend class SBF;
public:
    UncompressedBF(const std::string & filename, HashPair hp, int nh, uint64_t size = 0);
    UncompressedBF(const std::string & filename, BF* copy);
    virtual ~UncompressedBF();

    virtual void load();
    virtual void save();

    virtual int operator[](uint64_t pos) const;
    virtual void set_bit(uint64_t p);
    virtual uint64_t size() const;
    virtual uint64_t similarity(const BF* other, int type) const;
    virtual std::tuple<uint64_t, uint64_t> b_similarity(const BF* other) const;
    virtual BF* union_with(const std::string & new_name, const BF* f2) const;
    virtual void union_into(const BF* f2);
    virtual void union_into(const BF* f2, int type);
    virtual uint64_t count_ones() const;
    virtual uint64_t count_ones(int type) const;
    virtual void compress();

    virtual BF* sim_with(const std::string & new_name, const BF* f2) const;
    virtual void sim_into(const BF* f2);

    virtual BF* dif_with(const std::string & new_name, const BF* f2) const;
    virtual void dif_into(const BF* f2);

protected:
    sdsl::bit_vector* bv;
};

class SBF : public BF {
friend class UncompressedBF;
public:
    SBF(const std::string & filename, HashPair hp, int nh, uint64_t size = 0);
    SBF(const std::string & filename, BF* copy);
    SBF(const std::string & filename, BF* mask, BF* copy);
    virtual ~SBF();

    virtual void load();
    virtual void save();

    virtual bool contains(const jellyfish::mer_dna & m, int type) const;
    virtual bool contains(const size_t pos, int type) const;
    bool contains(const std::string & str, int type) const;
    //bool contains(const size_t & m, int type) const;

    virtual std::string get_sim_name();
    virtual std::string get_dif_name();

    virtual int operator[](uint64_t pos) const;
    virtual void set_bit(uint64_t p);
    virtual void set_difbit(uint64_t p);
    virtual void unset_bit(uint64_t p);
    virtual void unset_difbit(uint64_t p);
    virtual uint64_t size() const;
    virtual uint64_t size(int type) const;
    virtual uint64_t similarity(const BF* other, int type) const;
    virtual std::tuple<uint64_t, uint64_t> b_similarity(const BF* other) const;
    virtual BF* union_with(const std::string & new_name, const BF* f2) const;
    virtual void union_into(const BF* f2);
    virtual void union_into(const BF* f2, int type);
    virtual uint64_t count_ones() const;
    virtual uint64_t count_ones(int type) const;
    virtual void compress();
    virtual void compress(BF* rm);
    // Consider using bit_vectors for both.
    // Would like for both of these to be constant but having problems
    virtual void remove_duplicate(BF* f2);
    virtual void add_different(const sdsl::bit_vector & new_dif);
    // Finds the elements which were in this->sim and not in f2->sim
    virtual sdsl::bit_vector* calc_new_dif_bv(const BF* f2);
    // Accessor functions to perform simple and or xor operations on sim or dif filters
    virtual sdsl::bit_vector* calc_sim_bv(const BF* f2, int type);
    virtual sdsl::bit_vector* calc_dif_bv(const BF* f2, int type);
    virtual sdsl::bit_vector* calc_union_bv(const BF* f2, int type);

    virtual BF* sim_with(const std::string & new_name, const BF* f2) const;
    virtual void sim_into(const BF* f2);

    virtual BF* dif_with(const std::string & new_name, const BF* f2) const;
    virtual void dif_into(const BF* f2);

    virtual void update_mask(const BF* update);
    virtual void update_mask(const BF* u1, const BF* u2);
    //virtual bool contains(const jellyfish::mer_dna & m, int type) const;
    //bool contains(std::string str, int type) const;
protected:
    sdsl::bit_vector* sim;
    sdsl::bit_vector* dif;
};


sdsl::bit_vector* union_bv_fast(const sdsl::bit_vector & b1, const sdsl::bit_vector& b2);
BF* load_bf_from_file(const std::string & fn, HashPair hp, int nh);

// sim == and
sdsl::bit_vector* sim_bv_fast(const sdsl::bit_vector & b1, const sdsl::bit_vector& b2);
// dif == xor
sdsl::bit_vector* dif_bv_fast(const sdsl::bit_vector & b1, const sdsl::bit_vector& b2);
sdsl::bit_vector* copy_bv_fast(const sdsl::bit_vector & b1);

std::string split_filename(const std::string & new_name);
#endif
