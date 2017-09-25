#include "BloomTree.h"
#include "SplitBloomTree.h"
#include "util.h"
#include "BF.h"
#include "gzstream.h"

#include <fstream>
#include <list>
#include <cassert>
#include <jellyfish/file_header.hpp>

//Heap<const SplitBloomTree> SplitBloomTree::bf_cache;

// construct a bloom filter with the given filter backing.
SplitBloomTree::SplitBloomTree(
    const std::string & f, 
    HashPair hp,
    int nh
) :
    BloomTree(f, hp, nh)
    //filename(f),
    //hashes(hp),
    //num_hash(nh),
    //bloom_filter(0),
    //heap_ref(nullptr),
    //parent(0),
    //usage_count(0),
    //dirty(false)
{
    //children[0] = nullptr;
    //children[1] = nullptr;
}

// free the memory for this node.
SplitBloomTree::~SplitBloomTree() {
    unload();
}

// Frees the memory associated with the bloom filter
/*
void SplitBloomTree::unload() const { 
    // you can't unload something until you remove it from the cache
    // DEBUG 
    //std::cerr << "Unloading " << name() << std::endl;
    
    // free the memory
    if (bloom_filter != nullptr) {
        if (dirty) {
            bloom_filter->save();
        }
        delete bloom_filter; 
        bloom_filter = nullptr; 
    }
    dirty = false;
}
*/
// Loads the bloom filtering into memory
/* Remove this to try to eliminate overloaded methods
bool SplitBloomTree::load() const {
    if (bloom_filter == nullptr) {
        //std::cerr << "Loading BF: " << filename << std::endl;

        // if the cache isn't protected from deleting elements, remove enough
        // elements so that there is 1 cache spot free (if the cache is
        // protected, we're allowed to go over the cache limit)
        if(!bf_cache.is_protected()) SplitBloomTree::drain_cache();

        // read the BF file and set bloom_filter
        bloom_filter = load_bf_from_file(filename, hashes, num_hash);
        const SBF* b = dynamic_cast<const SBF*>(bloom_filter);
        if (b == nullptr) {
            std::cerr << filename << std::endl;
            DIE("Split Bloom Filter can only load SSBF.");
        }

        bloom_filter->load();
        heap_ref = bf_cache.insert(this, usage());
        dirty = false;

        // since we had to load, we bump up the usage to make it less likely we
        // load again in the near future.
        increment_usage();
    }
    return true;
}
*/
// Create a new node that is the union of the bloom filters
// in two other nodes;
// These nodes were not previously linked and thus we dont have any 'new differences' to add back in
SplitBloomTree* SplitBloomTree::union_bloom_filters(const std::string & new_name, SplitBloomTree* f2) {

    //std::cerr << "Cache size (start): " << cache_size() << std::endl;
    protected_cache(true);

    SBF* other_bf = dynamic_cast<SBF*>(f2->bf());
    if (other_bf == nullptr) {
        DIE("Split Bloom Filter bf() should be SBF.");
    }
    SBF* my_bf = dynamic_cast<SBF*>(bf());
    if (my_bf == nullptr) {
        DIE("Split Bloom Filter can only load SSBF.");
    }

    SplitBloomTree* bt = new SplitBloomTree(new_name, hashes, num_hash);

    // We first union the two key nodes and save it as a new node.
    //protected_cache(true);
    bt->bloom_filter = my_bf->union_with(new_name, other_bf); 
    //std::cerr << "Num 1s (New Filter Sim): " << bt->bf()->count_ones(0) << std::endl;
    //std::cerr << "Num 1s (New Filter Dif): " << bt->bf()->count_ones(1) << std::endl;
    //std::cerr << "Num 1s (Left Node Sim): " << my_bf->count_ones(0) << std::endl;
    //std::cerr << "Num 1s (Left Node Dif): " << my_bf->count_ones(1) << std::endl;
    //std::cerr << "Num 1s (Right Node Sim): " << other_bf->count_ones(0) << std::endl;
    //std::cerr << "Num 1s (Right Node Dif): " << other_bf->count_ones(1) << std::endl;
    bt->set_child(0, this);
    bt->set_child(1, f2);

    //Here we remove elements from the union's sim from each child's sim
    my_bf->remove_duplicate(bt->bf());
    other_bf->remove_duplicate(bt->bf());
    //std::cerr << "AFTER REMOVING DUPLICATES " << std::endl;    
    //std::cerr << "Num 1s (Left Node Sim): " << my_bf->count_ones(0) << std::endl;
    //std::cerr << "Num 1s (Left Node Dif): " << my_bf->count_ones(1) << std::endl;
    //std::cerr << "Num 1s (Right Node Sim): " << other_bf->count_ones(0) << std::endl;
    //std::cerr << "Num 1s (Right Node Dif): " << other_bf->count_ones(1) << std::endl;

    this->dirty = true;
    f2->dirty = true;

    //bf_cache.insert(bt, bt->usage());
    bt->dirty = true;
    bt->unload();

    //std::cerr << "Cache size (end): " << cache_size() << std::endl; 
    protected_cache(false);
    return bt; 
}

// Even when we are unioning into a node, we must still make a new filter
// This filter is temporary
// Could have zero or two children. (IT SHOULD BE IMPOSSIBLE TO HAVE ONE CHILD)
void SplitBloomTree::union_into(const SplitBloomTree* other) {
    //SplitBloomTree* temp = new SplitBloomTree("temporary.sim.bf.bv", hashes, num_hash);

    //std::cerr << "Cache size (start): " << cache_size() << std::endl;
    protected_cache(true);

    SBF* my_bf = dynamic_cast<SBF*>(bf());
    if (my_bf == nullptr) {
        DIE("Split Bloom Filter can only load SSBF.");
    }
    SBF* other_bf = dynamic_cast<SBF*>(other->bf());
    if (other_bf == nullptr) {
        DIE("Split Bloom Filter bf() should be SBF.");
    }
    // new differences are things which were similar at this node until other was added to the tree.
    // BITS WHICH ARE NEW TO THE TREE SOURCED FROM OTHER_BF ARE ADDED TO THE APPROPRIATE CHILD LATER
    sdsl::bit_vector* new_dif = my_bf->calc_new_dif_bv(other_bf); 

    // After determing what is different we can then perform the original union.
    //protected_cache(true);
    my_bf->union_into(other_bf);
    other_bf->remove_duplicate(my_bf); 
    dirty = true;
    other->dirty = true;
    //protected_cache(false);

    // We now add back elements which are no longer universally similar in the children
    // This is inefficient - if we know the insert path we could just place this at the lowest level of divergence.
    if (child(0) != nullptr){
        SBF* cbf0 = dynamic_cast<SBF*>(child(0)->bf());
        if (cbf0 == nullptr) {
            DIE("child bf() should be SBF.");
        }
        cbf0->add_different(*new_dif);
        child(0)->set_dirty(true);
    }
    if (child(1) != nullptr){
        SBF* cbf1 = dynamic_cast<SBF*>(child(1)->bf());
        if (cbf1 == nullptr) {
            DIE("child bf() should be SBF.");
        }
        cbf1->add_different(*new_dif);
        child(1)->set_dirty(true);
    }
    //std::cerr << "Cache size (end): " << cache_size() << std::endl; 
    protected_cache(false);

    delete new_dif;

}

/* read a file that defines the bloom tree structure. The
   file has lines of the form:
    Root
    *Child1
    ***Child3
    ****Child4
    *Child2
   where the '*' indicates the level and where "Child1" etc
   are comma-separated lists of attributes, the first one 
   of which must be the BF filename.

   This function will return a pointer to the root of the
   constructed bloom tree.
*/
SplitBloomTree* read_split_bloom_tree(const std::string & filename, bool read_hashes) {
    std::ifstream in(filename.c_str());

    std::list<SplitBloomTree*> path;
    SplitBloomTree* tree_root = 0;
    int n = 0;
    // if read_hashes is false, you must promise never to access the bloom filters
    HashPair* hashes = new HashPair; // useless hashpair used if read_hashes is false
    int num_hashes = 0;

    std::string node_info;
    while (getline(in, node_info)) {
        node_info = Trim(node_info);
        if (node_info.size() == 0) continue;
        size_t level = node_info.find_first_not_of("*");
        node_info.erase(0, level);

        // each node info is a comma separated list
        std::vector<std::string> fields;
        SplitString(node_info, ',', fields);
        std::string bf_filename = fields[0];
        //std::cerr << "Reading BN info: " << bf_filename << " level = " << level << std::endl;

        n++;

        SplitBloomTree* bn = nullptr;

        // if we're at the root
        if (path.size() == 0) {
            DIE_IF(level != 0, "Root must start in column 0");
            DIE_IF(tree_root != 0, "Can't set root twice!");

            // set the hash function up
            if (read_hashes) {
                DIE_IF(fields.size() < 2, "Must specify hash file for root.");
                hashes = get_hash_function(fields[1], num_hashes);
            }

            // create the root node
            bn = new SplitBloomTree(bf_filename, *hashes, num_hashes); 
            tree_root = bn;
            
        // if we're adding a child
        } else {
            bn = new SplitBloomTree(bf_filename, *hashes, num_hashes); 

            while (path.size() > level) {
                path.pop_back();
            }
            DIE_IF(level != path.size(), 
                "Must increase level by <= 1");

            if (path.back()->child(0) == nullptr) {
                path.back()->set_child(0, bn);
            } else if (path.back()->child(1) == nullptr) {
                path.back()->set_child(1, bn);
            } else {
                DIE("Tried to add >= 2 children to a node.");
            }
        }
        path.push_back(bn);
    }
    delete hashes;

    std::cerr << "Read " << n << " nodes in Bloom Tree" << std::endl;
    
    return tree_root;
}

void convert_sbt_filters(BloomTree* T, BF* cumul, std::string out_loc){
    SplitBloomTree* ST = dynamic_cast<SplitBloomTree*>(T);
    if (ST == nullptr){
        DIE("Could not convert BloomTree to SplitBT");
    }
    // handle the case of inserting into an empty tree
    if (T == nullptr) {
        DIE("Empty tree cannot be rebuilt.");
    }

    std::string base_name = test_basename(T->name(),".sim.bf.bv");
    std::string new_name = out_loc + "/" + base_name + ".bf.bv";
    if (ST->num_children() == 0) {
        BF* new_bf = new UncompressedBF(new_name, cumul);
        new_bf->union_into(ST->bf(),2);
        new_bf->save();
        delete new_bf;
    } else {
        BF* new_bf = new UncompressedBF(new_name, cumul);
        new_bf->union_into(ST->bf(),0); // Just union similarity
        convert_sbt_filters(ST->child(0), new_bf, out_loc);
        convert_sbt_filters(ST->child(1), new_bf, out_loc);
        delete new_bf;
    }

}


// Validates a compressed SSBT
void validate_SSBT(BloomTree* T){
    SplitBloomTree* ST = dynamic_cast<SplitBloomTree*>(T);
    if (ST==nullptr){
        DIE("Must use SplitBloomTree");
    }

    BF* root_bf = ST->bf();

    compressedSBF* rbf = dynamic_cast<compressedSBF*>(root_bf);
    if (rbf==nullptr){
        DIE("Failed to convert root to cbf");
    }   
    uint64_t sim_size = rbf->size(0);
    uint64_t dif_size = rbf->size(1);

    sdsl::rank_support_rrr<1,255> rbv_sim(rbf->sim_bits);
    uint64_t sim_ones = rbv_sim(sim_size);


    //Each filter's dif is sized by the number of ones in the sim
    if (sim_size!=sim_ones+dif_size){
        std::cerr << "Bad internal filters at: " << rbf->get_name() << std::endl;
        std::cerr << "sim_size: " << sim_size << std::endl;
        std::cerr << "dif_size: " << dif_size << std::endl;
        std::cerr << "sim_ones: " << sim_ones << std::endl;
        DIE("Dif filter not sized to sim_size - sim_ones");
    }
    
    // Validate children
    if (ST->child(0)){
        //Make sure children load as compressedSBF
        BF* child0_bf = ST->child(0)->bf();
        compressedSBF* c0bf = dynamic_cast<compressedSBF*>(child0_bf);
        if (c0bf==nullptr){
            DIE("Failed to convert child to cbf");
        } 

        BF* child1_bf = ST->child(1)->bf();
        compressedSBF* c1bf = dynamic_cast<compressedSBF*>(child1_bf);
        if (c1bf==nullptr){
            DIE("Failed to convert child to cbf");
        }

        if(rbf->get_name() != c0bf->get_name()){
            DIE("Somehow parent is left child!");
        }
        if(rbf->get_name() != c1bf->get_name()){
            DIE("Somehow parent is right child!");
        }

        // Make sure both children of my node have the same sim size
        sdsl::rank_support_rrr<0,255> rbv_dif(rbf->dif_bits);
        uint64_t dif_ones = rbv_dif(dif_size);
        uint64_t child0_size = c0bf->size(0);
        uint64_t child1_size = c1bf->size(0);
        assert(dif_ones >= 0);

        if(child0_size != child1_size){
            std::cerr <<"Bad external filters at: " << rbf->get_name() << std::endl;
            std::cerr<<"Parent dif filter: " << dif_size <<std::endl;
            std::cerr<<"Child 0 sim filter: " << child0_size <<std::endl;
            std::cerr<<"Child 1 sim filter: " << child1_size <<std::endl;
            DIE("Child sim size mismatch");
        }

        // Make sure sim filters are exactly equal to dif-dif ones
        if(dif_size!=child0_size+dif_ones){
            std::cerr <<"Child not sized properly at: " << rbf->get_name() <<std::endl;
            std::cerr<<"Parent dif filter: " << dif_size <<std::endl;
            std::cerr<<"Parent dif num ones: " << dif_ones <<std::endl;
            std::cerr<<"Child 0 sim filter: " << child0_size <<std::endl;
            std::cerr<<"Child 1 sim filter: " << child1_size <<std::endl;
            DIE("Sim + dif ones should equal dif filter");
        }

        validate_SSBT(ST->child(0));
        validate_SSBT(ST->child(1));
    }
    std::cerr << "Validated " << rbf->get_name() << " (Size: " << rbf->size(0) <<")"  <<std::endl;
}

/*
void write_bloom_tree_helper(std::ostream & out, SplitBloomTree* root, int level=1) {
    std::string lstr(level, '*');

    for (int i = 0; i < 2; i++) {
        if (root->child(i) != nullptr) {
            out << lstr << root->child(i)->name() << std::endl;
            write_bloom_tree_helper(out, root->child(i), level+1);
        }
    }
}

// write the bloom tree file format in a way that can be read by
// read_bloom_tree()
void write_bloom_tree(
    const std::string & outfile, 
    SplitBloomTree* root, 
    const std::string & matrix_file
) {
    std::cerr << "Writing to " << outfile << std::endl;
    std::ofstream out(outfile.c_str());
    out << root->name() << "," << matrix_file << std::endl;
    write_bloom_tree_helper(out, root);
    std::cerr << "Done." << std::endl;
}

void write_compressed_bloom_tree_helper(std::ostream & out, SplitBloomTree* root, int level=1) {
    std::string lstr(level, '*');

    for (int i = 0; i < 2; i++) {
        if (root->child(i) != nullptr) {
            out << lstr << root->child(i)->name() << ".rrr" << std::endl;
            write_compressed_bloom_tree_helper(out, root->child(i), level+1);
        }
    }
}

// write the bloom tree file format in a way that can be read by
// read_bloom_tree()
void write_compressed_bloom_tree(
    const std::string & outfile,
    SplitBloomTree* root,
    const std::string & matrix_file
) {
    std::cerr << "Writing to " << outfile << std::endl;
    std::ofstream out(outfile.c_str());
    out << root->name() << ".rrr," << matrix_file << std::endl;
    write_compressed_bloom_tree_helper(out, root);
    std::cerr << "Done." << std::endl;
}
*/
