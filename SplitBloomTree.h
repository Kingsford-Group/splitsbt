#ifndef SPLITBLOOMTREE_H
#define SPLITBLOOMTREE_H

#include <string>
#include <queue>
#include "Heap.h"
#include "BF.h"

// this is the max number of BF allowed in memory at once.
extern int BF_INMEM_LIMIT;

class SplitBloomTree : public BloomTree {
public:
    SplitBloomTree(const std::string & f, HashPair hp, int nh);
    virtual ~SplitBloomTree();

    virtual SplitBloomTree* union_bloom_filters(const std::string & new_name, SplitBloomTree* f2);
    virtual void union_into(const SplitBloomTree* other);

    //virtual int usage() const;
    //virtual void increment_usage() const;
    //virtual void set_usage(int val) const;
    //virtual int cache_size();
    //virtual void set_heap_ref(Heap<const SplitBloomTree>::heap_reference* hr);
    //virtual Heap<const SplitBloomTree>::heap_reference* get_heap_ref();
    //static void protected_cache(bool b);

protected:
    bool load() const;
    void unload() const;
};


//Keeps one filename but uses "sim" and "diff" suffixes? 
/*class splitBloomTree : public BloomTree {
public:
    splitBloomTree(const std::string & f, HashPair hp, int nh);
    ~splitBloomTree();
    std::string name() const;

    uint64_t similarity(BloomTree* other, int type) const;
    std::tuple<uint64_t, uint64_t> b_similarity(BloomTree* other) const;
    BF* bf() const;

    BloomTree* union_bloom_filters(const std::string & new_name, BloomTree* f2);
    void union_into(const BloomTree* other);

    int usage() const;
    void increment_usage() const;
    static void protected_cache(bool b);

private:
    bool load() const;
    void unload() const;

    static Heap<const BloomTree> bf_cache;
    static void drain_cache();

    std::string filename;
    HashPair hashes;
    int num_hash;
    mutable BF* bloom_filter;
    mutable Heap<const BloomTree>::heap_reference* heap_ref;

    BloomTree* children[2];
    BloomTree* parent;
    mutable int usage_count;
    mutable bool dirty;
}
*/

SplitBloomTree* read_split_bloom_tree(const std::string & filename, bool read_hashes=true);
/*
void write_bloom_tree(const std::string & outfile, SplitBloomTree* root, const std::string & matrix_file);
void write_compressed_bloom_tree(const std::string & outfile, SplitBloomTree* root, const std::string & matrix_file);
*/
#endif
