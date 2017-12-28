#ifndef BLOOMTREE_H
#define BLOOMTREE_H

#include <string>
#include <queue>
#include "Heap.h"
#include "BF.h"

// this is the max number of BF allowed in memory at once.
extern int BF_INMEM_LIMIT;

class BloomTree {
public:
    BloomTree(std::string & f);
    BloomTree(const std::string & f, HashPair hp, int nh);
    virtual ~BloomTree();
    virtual std::string name() const;

    virtual BloomTree* child(int which) const;
    virtual void set_child(int which, BloomTree* c);
    virtual int num_children() const;
    virtual void set_parent(const BloomTree* p);
    virtual const BloomTree* get_parent() const;
    virtual uint64_t similarity(BloomTree* other, int type) const;
    virtual std::tuple<uint64_t, uint64_t> b_similarity(BloomTree* other) const;
    virtual BF* bf() const;

    virtual BloomTree* union_bloom_filters(const std::string & new_name, BloomTree* f2);
    virtual void union_into(const BloomTree* other);

    virtual int usage() const;
    virtual void increment_usage() const;
    virtual void set_usage(int val) const;
    virtual int cache_size() const;
    virtual void set_dirty(bool b) const;
    static void protected_cache(bool b);

    virtual HashPair get_hashes() const;
    virtual int get_num_hash() const;

protected:
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
HashPair* get_hash_function(const std::string & matrix_file, int & nh);
BloomTree* read_bloom_tree(const std::string & filename, bool read_hashes=true);
void write_bloom_tree(const std::string & outfile, BloomTree* root, const std::string & matrix_file);
void write_compressed_bloom_tree(const std::string & outfile, BloomTree* root, const std::string & matrix_file);
void convert_bloom_to_build(BloomTree* root, const std::string & outfile);
#endif
