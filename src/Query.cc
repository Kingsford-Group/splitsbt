#include "Query.h"
#include "Kmers.h"
#include "util.h"
#include <cassert>

float QUERY_THRESHOLD = 0.9;

// ** THIS IS NOW PARTIALLY DEPRICATED. ONLY WORKS WITH HARDCODED SIMILARITY TYPE
void assert_is_union(BloomTree* u) {
    BloomTree::protected_cache(true);
    BF* ubf = u->child(0)->bf()->union_with("union", u->child(1)->bf());
    BloomTree::protected_cache(false);

    // if the similaritity isn't 100%, we stop
    std::ostringstream oss;
    uint64_t sim = ubf->similarity(u->bf(),0);
    if (sim != ubf->size()) {
        std::cerr << "Filter at " << u->name() << " is not the union of its two children!" << std::endl;
        std::cerr << "Sim= " << sim << "Size= " << ubf->size() << std::endl;
        std::cerr << "Children:" << u->child(0)->name() << " " << u->child(1)->name() << std::endl;
        std::cerr << std::endl;
    } else {
        std::cerr << "Filter at " << u->name() << " looks good." << std::endl;
        std::cerr << "Children:" << u->child(0)->name() << " " << u->child(1)->name() << std::endl;
        std::cerr << std::endl;
    }
    delete ubf;
}

void check_bt(BloomTree* root) {
    if (root == nullptr) return;

    if (root->child(0) && root->child(1)) {
        assert_is_union(root);
    }

    check_bt(root->child(0));
    check_bt(root->child(1));
}


void draw_bt_recur(BloomTree* root, std::ostream& out) {
    if (root == nullptr) return;

    std::string current = quote(test_basename(root->name(), ".bf.bv"));
    
    if (root->child(0)) {
        out << current << " -> " << quote(test_basename(root->child(0)->name(), ".bf.bv")) << " ; " << std::endl;
        draw_bt_recur(root->child(0), out);
    } 
    if (root->child(1)) {
        out << current << " -> " << quote(test_basename(root->child(1)->name(), ".bf.bv")) << " ; " << std::endl;
        draw_bt_recur(root->child(1), out);
    } 
}

void draw_bt(BloomTree* root, std::string outfile) {
    std::ofstream out(outfile.c_str());
    DIE_IF(!out, "Couldn't open output file");
    out << "digraph BloomTree {" << std::endl;
    draw_bt_recur(root, out);
    out << "}" << std::endl;
}

void popcount_bt(BloomTree* root) {
	if (root == nullptr) return;

	std::cerr << root->bf()->count_ones(0) << " " << root->bf()->size(0) << ", ";

	if (root->child(0)) {
        std::cerr << root->bf()->count_ones(1) << " " << root->bf()->size(1) << std::endl;
		popcount_bt(root->child(0));
	} else{
        std::cerr << std::endl;
    }
	if (root->child(1)){
		popcount_bt(root->child(1));
	}
}

void compress_bt(BloomTree* root) {
	if (root == nullptr) return;

	root->bf()->compress();

	if (root->child(0)) {
		compress_bt(root->child(0));
	}
	if (root->child(1)){
		compress_bt(root->child(1));
	}
}

void compress_splitbt(BloomTree* root, sdsl::bit_vector* noninfo){
    root->bf()->compress(noninfo);

    if (root->child(0)){
        //update noninformative bits
        root->bf()->get_noninfo(noninfo);
        compress_splitbt(root->child(0),noninfo);
        root->bf()->get_noninfo(noninfo);
        compress_splitbt(root->child(1),noninfo);
    }
}

void compress_splitbt(BloomTree* root, BF* rbf){
    std::cerr << "Cache Size: " << root->cache_size() <<std::endl;

    if (root == nullptr) {
        DIE("We should no longer pass through leaves.");
        return;
    }
   
    //SplitBloomTree* - We might not need to dynamic cast this
 
    SBF* remove_bf = dynamic_cast<SBF*>(rbf);
    if (remove_bf == nullptr){
        DIE("remove_bf doesn't really need to be SBF but it is for now.");
    }
    
    //compress
    // Test whether compression is memory leak
    root->bf()->compress(remove_bf);

    // Adjust rbf for next level
    // If we are a leaf, no need to adjust filter
    if(root->child(0)){
        //debug
        if (root->name()==root->child(0)->name()){
            DIE("Node cannot be its own parent!");
        }
        if (root->name()==root->child(1)->name()){
            DIE("Node cannot be its own parent!");
        }
        //root->protected_cache(true);
        //add current, compress children
        remove_bf->update_mask(root->bf());
    
        if(root->child(0)){
            compress_splitbt(root->child(0), remove_bf);
        }
        if(root->child(1)){
            compress_splitbt(root->child(1), remove_bf);
        }
        
        root->protected_cache(true);
        //remove current
        if (root->get_parent() != nullptr){
            BloomTree* mySiblingBT;
            //std::cerr << "My name is: '" << root->name() << "'\n";
            //std::cerr << "My bloom filter's name is: '" << root->bf()->get_name() << "'\n"; 
            if (root->name() == root->get_parent()->child(0)->name()){
                mySiblingBT = root->get_parent()->child(1);
                //std::cerr << "My name is: '" << root->get_parent()->child(0)->name() << "'\n";
                //std::cerr << "My sibling's name is: '" << root->get_parent()->child(1)->name() << "'\n";
                BF* myStupidBF = root->bf();
                BF* myStupidSibling = mySiblingBT->bf();
                //std::cerr << "My stupid name is: " << myStupidBF->get_name() <<std::endl;
                //std::cerr << "My stupid sibling's name is: " << myStupidSibling->get_name() <<std::endl;
                remove_bf->update_mask(myStupidBF, myStupidSibling);
                //remove_bf->update_mask(root->bf(), root->get_parent()->child(1)->bf());
            } else {
                mySiblingBT = root->get_parent()->child(0);
                //std::cerr << "My name is: '" << root->get_parent()->child(1)->name() << "'\n";
                //std::cerr << "My siblin's name is: '" << root->get_parent()->child(0)->name() << "'\n";
                assert(root->get_parent()->child(1)->name()==root->name());
                BF* myStupidBF = root->bf();
                BF* myStupidSibling = mySiblingBT->bf();
                //std::cerr << "My stupid name is: " << myStupidBF->get_name() <<std::endl;
                //std::cerr << "My stupid sibling's name is: " << myStupidSibling->get_name() <<std::endl;
                remove_bf->update_mask(myStupidBF, myStupidSibling);
                
                //remove_bf->update_mask(root->bf(), root->get_parent()->child(0)->bf());
            }
        } else{
            std::cerr<< "GLOBAL ROOT NAME: " << root->name() <<std::endl;
        }
        root->set_usage(0);
        root->protected_cache(false);
    }
} 

bool query_passes(BloomTree* root, const std::set<jellyfish::mer_dna> & q) {
    assert(q.size() > 0);
    auto bf = root->bf();
    unsigned c = 0;
    for (const auto & m : q) {
        //DEBUG: std::cout << "checking: " << m.to_str();
        if (bf->contains(m)) c++;
        //DEBUG: std::cout << c << std::endl;
    }
    return (c >= QUERY_THRESHOLD * q.size());
}

bool query_passes(BloomTree* root, batchInfo* bi){
    auto bf = root->bf();
    bool has_children = root->child(0) || root->child(1);
    float c, min_pass;    
   
    // for each splitQueryInfo 
    for (const auto sqi : bi->sqi_list){
        // for each queryInfo 
        min_pass = 0.0;
        for (const auto & q : sqi->partial_queries){
            c = 0.0;

            min_pass += q->total_kmers * q->q_thresh;
            
            // determine matching in query
            for (const auto & m : q->query_kmers){
                if(bf->contains(m)) c++;
            }
        }
    
        // if 'not contains' is true or 'contains' is false, return false
        if (c >= min_pass){
            if(!has_children && sqi->type == 1) return false;
        } else{
            if(sqi->type == 0) return false;
        }
    
    } 
    return true;
}
// return true if the filter at this node contains > QUERY_THRESHOLD kmers
// XXX: Fix query to work for BF, SBF, compressedSBF
bool query_passes(BloomTree* root, QueryInfo*  q) {//const std::set<jellyfish::mer_dna> & q) {
    float weight = 1.0;
    //assert(q.size() > 0);
    auto bf = root->bf();
    float c = 0;
    unsigned n = 0;
    bool weighted = 0;

    SplitBloomTree* sroot = dynamic_cast<SplitBloomTree*>(root);
    if (sroot == nullptr) { // start of normal query
        if (q->weight.empty()){
	        weighted=0;
        } else { weighted = 1; }
        for (const auto & m : q->query_kmers) {
        //DEBUG: std::cout << "checking: " << m.to_str();
        	if (weighted){
	    	    if(q->weight.size() > n){ 
		        	weight=q->weight[n]; 
        		} else {
    			std::cerr << "Number of weights (" << q->weight.size() <<") less than query kmers (" << q->query_kmers.size() << ")."  << std::endl;
			    exit(3);
		        }
	        }
            if (bf->contains(m)) c+=weight;
	        n++;
            //DEBUG: std::cout << c << std::endl;
        }

    } else {// end of normal
        // tail index 
        c = q->matched_kmers;
        std::cerr << "Prematch: " << c << ", Tail Index: " << q->tail_index <<std::endl;
        //std::cerr << bf->size(0) << " " << bf->size(1) << std::endl;
        //std::set<jellyfish::mer_dna> passedKmers;
        //for (const auto & m : q->query_kmers) {
        // We dont add or delete elements. Re-arranging their order based on tail index
        // should allow storage of a single int rather then the full kmer set
        for (int i = 0; i <= q->tail_index; i++) {
            auto & m = q->query_kmers[i];
            // We can break whenever the total matched is above threshold.
            //if (q->matched_kmers >= QUERY_THRESHOLD * q->total_kmers){
            //    return true;
            //}
            //std::cerr << i << " " << q->tail_index << std::endl;
            compressedSBF* cbf = dynamic_cast<compressedSBF*>(bf);
            if (cbf == nullptr){
                DIE("Could not convert to compressedSBF");
            }
            if (m >= cbf->size(0)){
                std::cerr << "Tail Index: " << i << " " << q->tail_index << std::endl;
                std::cerr << m << " " << cbf->size(0) << " " << cbf->size(1) << std::endl;
                DIE("Hash value out of bounds");
            }
        
            if (bf->contains(m,0)) { //kmer found in similarity filter
                q->matched_kmers+=weight; //record the hit (weighted for future weighted functionality)
                swap_kmer_position(q->query_kmers,i, q->tail_index);
                //something we havent seen now occupies the same position as us
                //something we HAVE seen now occupies the tail position. 
                //decrease both so we stay on same index but don't look at tail again.
                q->tail_index--;
                i--;
                //q->query_kmers.erase(m); //we no longer need to store kmer.
                c+=weight;
            } else if (bf->contains(m,1)) { //kmer found in difference filter
                c+=weight; //treat it like regular SBT query.
                //passedKmers.insert(m); //Insert kmers which exist somewhere (and weren't in sim)
            } else{ //kmer doesnt exist at this depth of the tree
                //we remove it like the sim filter but dont add to matched_kmers
                swap_kmer_position(q->query_kmers,i, q->tail_index);
                q->tail_index--;
                i--;
            }
            
        }
        //q->query_kmers=passedKmers;
    }
    //std::cerr << root->name() << std::endl;
    std::cerr << "Total weight: " << c << ", thresh: " << QUERY_THRESHOLD << ", maxkmer: " << q->total_kmers << std::endl;
    return (c >= q->q_thresh * q->total_kmers);
}

// recursively walk down the tree, proceeding to children only
// if their parent passes the query threshold; 
void query(
    BloomTree* root, 
    const std::set<jellyfish::mer_dna> & q, 
    std::vector<BloomTree*> & out
) {
    root->increment_usage();
    if (query_passes(root, q)) {
        //DEBUG: std::cout << "passed at " << root->name() << std::endl;
        int children = 0;
        if (root->child(0)) {
            query(root->child(0), q, out);
            children++;
        }
        if (root->child(1)) {
            query(root->child(1), q, out);
            children++;
        }
        if (children == 0) {
            out.push_back(root);
        }
    }
}

// same as query() but the string is first converted into a set of kmers.
void query_string(
    BloomTree* root, 
    const std::string & q,
    std::vector<BloomTree*> & out
) {
    query(root, kmers_in_string(q), out);
}

// read 1 query per line, execute it, and print to the output stream o the
// results in the format:
//      *QUERY number_results
//      BF names
//
void query_from_file(
    BloomTree* root, 
    const std::string & fn,
    std::ostream & o
) { 
    std::vector<BloomTree*> out;
    std::string line;

    std::ifstream in(fn);
    DIE_IF(!in.good(), "Couldn't open query file.");
    while (getline(in, line)) {
        line = Trim(line);
        if (line.size() < jellyfish::mer_dna::k()) continue;

        o << "*" << line;
        query_string(root, line, out);
        o << " " << out.size() << std::endl;

        for (const auto& n : out) {
            o << n->name() << std::endl;
        }

        out.clear();
    }
}

/******
 * Batch querying
 ******/
void print_query_results(const batchSet & bs, std::ostream & out) {
    int count = 0;
    for (auto& bi : bs) {
        //out << "split query: " << count << std::endl;
        out << "*";
        for(auto & sqi : bi->sqi_list){
            if (sqi->type==0) out << "+";
            else if (sqi->type==1) out << "-";
            for (auto & q : sqi->partial_queries){
                out << "[" << q->query << "]";
            }
        }
        out << " " << bi->matching.size() << std::endl;
        for (const auto& n : bi->matching) {
            out << n->name() << std::endl;
        }
        count++;
    }
}

void print_query_results(const QuerySet & qs, std::ostream & out) {
    for (auto& q : qs) {
        out << "*" << q->query << " " << q->matching.size() << std::endl;
        for (const auto& n : q->matching) {
            out << n->name() << std::endl;
        }
    }
}


void split_query_batch(BloomTree* root, batchSet & bs){
    bool has_children = root->child(0) || root->child(1);

    batchSet pass;
    unsigned n = 0;
    for (auto & bi : bs) {
        if (query_passes(root, bi)) { //q->query_kmers)) {
            if (has_children) {
                pass.emplace_back(bi);
            } else {
                bi->matching.emplace_back(root);
                n++;
            }
        }
    }
     if (has_children) { //Changing format
            std::cout << root->name() << " internal " << pass.size() << std::endl;
        } else {
            std::cout << root->name() << " leaf " << n << std::endl;
        }

        if (pass.size() > 0) {
            if (root->child(0)) {
                split_query_batch(root->child(0), pass);
            }

            if (root->child(1)) {
                split_query_batch(root->child(1), pass);
            }
        }
}

void query_batch(BloomTree* root, QuerySet & qs) {
    //std::cerr << "Batch Query!" << std::endl;

    SplitBloomTree* sroot = dynamic_cast<SplitBloomTree*>(root);
    if (sroot == nullptr) { //the standard query_batch.
        // how many children do we have?
        std::cerr << "Old Bloom Tree" << std::endl;
        bool has_children = root->child(0) || root->child(1);

        // construct the set of queries that pass this node
        QuerySet pass;
        unsigned n = 0;
        for (auto & q : qs) {
            if (query_passes(root, q)) { //q->query_kmers)) {
                if (has_children) {
                    pass.emplace_back(q);
                } else {
                    q->matching.emplace_back(root);
                    n++;
                }
            } 
        }

        // $(node name) $(internal / leaf) $(number of matches)
        if (has_children) { //Changing format
            std::cout << root->name() << " internal " << pass.size() << std::endl;
        } else {
            std::cout << root->name() << " leaf " << n << std::endl;
        }

        if (pass.size() > 0) {
            // if present, recurse into left child
            if (root->child(0)) {
                query_batch(root->child(0), pass);
            }
            
            // if present, recurse into right child
            if (root->child(1)) {
                query_batch(root->child(1), pass);
            }
        }
    } else { //Right now we have two general query types.
        QuerySet pass;
        //QuerySet copy;
        unsigned n = 0;
        bool search_flag = false;
        int num_search = 0;
        std::cerr << "Split Bloom Tree " << root->bf()->get_name() << std::endl;
        bool has_children = root->child(0) || root->child(1);
        for (auto & q : qs) {
            //If query has enough matched kmers (from sim), don't waste time searching
            //std::cerr << "TAIL INDEX: " << q->tail_index << std::endl;
            //std::cerr << q->matched_kmers << " " << q->total_kmers * QUERY_THRESHOLD << std::endl;
            if (q->matched_kmers >= QUERY_THRESHOLD * q->total_kmers){
                //std::cerr << "Skipped search (enough similar kmers already found)" << std::endl;
                if (has_children) {
                    pass.emplace_back(q);
                } else {
                    q->matching.emplace_back(root);
                    n++;
                }
            } else if (query_passes(root, q)) { //q->query_kmers)) {
                //std::cerr << "Query passes \n";
                num_search++;
                if (has_children) {
                    pass.emplace_back(q);
                } else {
                    q->matching.emplace_back(root);
                    n++;
                }
            } else{ //Things which did not pass were still searched in query_passes
                num_search++;
            }
        }

        //Copy the queries that pass locally
        //We can then restore certain values after a depth traversal.
        //Each node is storing the query_kmers and # matching at that position
        std::vector<int> mk_vector; //matched_kmers
        std::vector<int> ti_vector; //tail_index
        //std::vector<std::vector<size_t>> query_vector;
        //std::vector<std::vector<size_t>> adj_query_vector;
        //std::vector<std::set<jellyfish::mer_dna>> qk_vector;
        //std::vector<std::string> qs_vector;

        // Check to see if we have to load files or not 
        for (auto qc : pass){
            if (qc->matched_kmers < QUERY_THRESHOLD * qc->total_kmers){
                search_flag=true;
            }
        }

        // *** Report number of queries which passed and how many searches were needed ***
        // $(node name) $(internal / leaf) $(number of passes) $(number of searches)
        if (has_children) { 
            std::cout << root->name() << " internal " << pass.size() << " " << num_search << std::endl;
        } else {
            std::cout << root->name() << " leaf " << n << " " << num_search << std::endl;
        }      
 
        // If we need ANY query to process, we have to load this current filter
        // and we have to build rank_support vectors
        compressedSBF* cbf;
        if (search_flag){
            cbf = dynamic_cast<compressedSBF*>(sroot->bf());
            if (cbf == nullptr){
                DIE("Add cases later to clean this up");
            }

            sdsl::rank_support_rrr<1,255> rbv_sim(cbf->sim_bits);
            sdsl::rank_support_rrr<0,255> rbv_dif(cbf->dif_bits);

            // Everything in pass was a hit match in the dif. 
            // A hit match in the sim is recorded by tail_index and matched_kmer

            // *** Store matched_kmers and tail_index *** 
            for (auto qc : pass) {
            //Because we can't know if mk is needed or not we'll just keep it always
                mk_vector.emplace_back(qc->matched_kmers);
                ti_vector.emplace_back(qc->tail_index);
                //query_vector.emplace_back(qc->query_kmers);
                //std::cerr << "Internal loop TAIL INDEX: " << qc->tail_index << std::endl;
            
                if (qc->matched_kmers >= QUERY_THRESHOLD * qc->total_kmers){
                    //pass
                } else{ //Only adjust queries which still need to be checked
                    for (int i = 0; i <= qc->tail_index; i++) {
                        auto m = qc->query_kmers[i];
                        assert(m <= cbf->size(0)); 
                        size_t sim_ones = rbv_sim(m);
                        assert(m-sim_ones <= cbf->size(1));
                        size_t dif_ones = rbv_dif(m-sim_ones);
                        qc->query_kmers[i]=m-sim_ones-dif_ones;
                        //std::cerr << i << " " << sim_ones << " " << dif_ones << " " << qc->query_kmers[i] << std::endl;  
                    }
                }
            } //for (auto qc : pas) bracket
        } // Search flag    


        int copy_it = 0;

        if (pass.size() > 0) {
            // if present, recurse into left child
            //std::cerr << "Non-zero pass size: " << pass.size() << std::endl;
            if (root->child(0)) {
                //std::cerr << "Querying left child \n";
                query_batch(root->child(0), pass);
            }

            // If search_flag is false, we made no edits and don't need to save any values
            if (search_flag){
                if (mk_vector.size() != pass.size()){
                    DIE("DELETED QUERY NEEDS TO BE RESTORED");
                }
                copy_it = 0;
                for (auto & q : pass){
                    if (mk_vector[copy_it] >= QUERY_THRESHOLD * q->total_kmers){
                        //This is debug and can be removed later
                        if (q->matched_kmers!=mk_vector[copy_it]){
                            DIE("Touched a passable query's matched_kmers!");
                        }
                        if (q->tail_index!=ti_vector[copy_it]){
                            DIE("Touched a passable query's tail index!");
                        }
                        //pass
                    } else{ 
                        q->matched_kmers=mk_vector[copy_it];//temp->matched_kmers;
                        q->tail_index=ti_vector[copy_it];
                    }
                    copy_it++; 
                }
            }

            // if present, recurse into right child
            if (root->child(1)) {
                //std::cerr << "Querying right child \n";
                query_batch(root->child(1), pass);
            }
            //const BloomTree* r_parent = root->get_parent();
            /*
            if (r_parent){
                std::cerr << "Tree position: \n";
                std::cerr << root->name() << " " << r_parent->name() << std::endl;
                std::cerr << root->child(0)->name() << " " << root->child(1)->name() << std::endl;
            } 
            */
            /*if (r_parent){
                copy_it=0;
                for (auto & q : pass){
                    std::cerr << "Restoring query_kmers" << std::endl;
                    q->query_kmers=query_vector[copy_it];
                    copy_it++;
                }
            }
            */

            // search_flag records if we needed to load cbf. If we did above, we need to now
            if (search_flag){
                copy_it = 0;
                cbf = dynamic_cast<compressedSBF*>(sroot->bf());
                sdsl::select_support_rrr<0,255> sbv_sim(cbf->sim_bits);
                sdsl::select_support_rrr<1,255> sbv_dif(cbf->dif_bits);
                for (auto & q : pass) {
                    if (mk_vector[copy_it] >= QUERY_THRESHOLD * q->total_kmers){
                        //pass
                    } else {
                        for (int i = 0; i <= ti_vector[copy_it]; i++) {
                            auto & m = q->query_kmers[i];
                            if (m >= cbf->dif_bits->size()){
                                std::cerr << "m: " << m << std::endl;
                                std::cerr << "cbf->dif_bits->size(): " << cbf->dif_bits->size() << std::endl;
                                assert(m<cbf->dif_bits->size());
                            }
                            size_t dif_pos = sbv_dif(m+1);
                            assert(dif_pos<cbf->sim_bits->size());
                            q->query_kmers[i]=sbv_sim(dif_pos+1);
                        }
                    }
                copy_it++;
               }
            } 
            /*
                    std::set<size_t> adjs(query_vector[copy_it].begin(), query_vector[copy_it].end());
                    std::set<size_t> qks(q->query_kmers.begin(), q->query_kmers.end());
                    std::set<size_t>::iterator j = qks.begin();
                    int cnt =0;
                    for(std::set<size_t>::iterator i=adjs.begin(); i!=adjs.end(); ++i){
                        if (*i != *j){
                            std::cerr << "Error at rollback kmer "<<  cnt << std::endl;
                            std::cerr << *i << " " << *j << std::endl;
                            std::cerr << cbf->sim_bits->size() << std::endl;
                            std::cerr << cbf->dif_bits->size() << std::endl;
                            DIE("Should match up!");
                        }
                        ++j;
                        ++cnt;
                    }
                    copy_it++;
                }
            */
            //}
        
        } // if (pass.size > 0) 
    
        std::cerr << "Completed batch_query instance " << root->bf()->get_name() << std::endl;
        root->set_usage(0);
    } // else case
} // function


void query_leaves (BloomTree* root, QuerySet & qs) {
    // how many children do we have?
    bool has_children = root->child(0) || root->child(1);

    // construct the set of queries that pass this node
	// But only for leaf nodes
	unsigned n=0;
	if (!has_children) {
    		for (auto & q : qs) {
		        if (query_passes(root, q)) {
		                q->matching.emplace_back(root);
		                n++;
       	 		}
    		}
	}

    // $(node name) $(internal / leaf) $(number of matches)
    if (!has_children) { //Changing format	
        std::cout << root->name() << " leaf " << n << std::endl;
    }

        // if present, recurse into left child
        if (root->child(0)) {
            query_leaves(root->child(0), qs);
        }

        // if present, recurse into right child
        if (root->child(1)) {
            query_leaves(root->child(1), qs);
        }
    
}

void batch_query_from_file(
    BloomTree* root, 
    const std::string & fn,
    std::ostream & o
) { 
    // read in the query lines from the file.
    std::string line;
    QuerySet qs;
    std::ifstream in(fn);
    DIE_IF(!in.good(), "Couldn't open query file.");
    std::size_t n = 0;
    while (getline(in, line)) {
        line = Trim(line);
        if (line.size() < jellyfish::mer_dna::k()) continue;
        qs.emplace_back(new QueryInfo(root->bf(), line));
        n++;
    }
    in.close();
    std::cerr << "Read " << n << " queries." << std::endl;

    // batch process the queries
    query_batch(root, qs);
    print_query_results(qs, o);

    // free the query info objects
    for (auto & p : qs) {
        delete p;
    }
}


// split on line
// split on '+'/'-'
void batch_splitquery_from_file(
    BloomTree* root, 
    const std::string & fn,
    std::ostream & o
) { 
    // read in the query lines from the file.
    std::string line;
    batchSet bs;
    std::ifstream in(fn);
    DIE_IF(!in.good(), "Couldn't open query file.");
    std::size_t n = 0;
    while (getline(in, line)) {
        //splitQuerySet* sqs = new splitQuerySet;
        line = Trim(line);
        // This boundary should be adjusted but really matters later
        if (line.size() < jellyfish::mer_dna::k()) continue;
        
        batchInfo* bi = new batchInfo(root->bf(),line);
/* 
        std::size_t left_end = 0;
        std::size_t right_end = 0;
        std::size_t index = 0;
        int next_type = 0;
        std::stringstream ss(line);
        char c;
        while (ss >> c){
            if (c == '+'){
                left_end = right_end;
                right_end = index;
                sqs->emplace_back(new splitQueryInfo(root->bf(), line.substr(left_end, right_end-left_end+1),next_type));
                next_type = 0;
                n++;
            } else if (c == '-'){
                left_end = right_end;
                right_end = index;
                sqs->emplace_back(new splitQueryInfo(root->bf(), line.substr(left_end, right_end-left_end+1),next_type));
                next_type = 1;
                n++;
            }
            index++;
        }

        left_end = right_end;
        right_end = index;
        //std::cerr << left_end << " " << right_end << std::endl;
        //std::cerr << "Substr: " << line.substr(left_end, right_end-left_end+1) <<  std::endl;
        sqs->emplace_back(new splitQueryInfo(root->bf(), line.substr(left_end, right_end-left_end+1),next_type));
*/
        bs.emplace_back(bi);
        n++;
    }
    in.close();
    std::cerr << "Read " << n << " queries." << std::endl;

    // batch process the queries
    split_query_batch(root, bs);
    print_query_results(bs, o);

    // free the query info objects
    int count = 0;
    for (auto & bi : bs) {
        std::cerr << "Query " << count << ": " << std::endl;
        for (auto & sqi : bi->sqi_list){
            for (auto & q : sqi->partial_queries){
                std::cerr << "Total kmers " << q->total_kmers << std::endl;
                std::cerr << "query: " << q->query << std::endl;
                delete q;
            }
            delete sqi;
        }
        count++;
        delete bi;
    }
}
void batch_weightedquery_from_file(
    BloomTree* root,
    const std::string & fn,
    const std::string & wf,
    std::ostream & o
) {
    // read in the query lines from the file.
    std::string line;
    std::string wfline; 
    QuerySet qs;
    std::ifstream in(fn);
    std::ifstream wfin(wf);
    DIE_IF(!in.good(), "Couldn't open query file.");
    DIE_IF(!wfin.good(), "Couldn't open weight file.");
    std::size_t n = 0;
    while (getline(in, line)) {
	getline(wfin, wfline);
        line = Trim(line);
	wfline = Trim(wfline);
	
        if (line.size() < jellyfish::mer_dna::k()) continue;
        qs.emplace_back(new QueryInfo(root->bf(), line, wfline));
        n++;
    }
    in.close();
    std::cerr << "Read " << n << " queries." << std::endl;

    // batch process the queries
    query_batch(root, qs);
    print_query_results(qs, o);

    // free the query info objects
    for (auto & p : qs) {
        delete p;
    }
}


void leaf_query_from_file(
	BloomTree* root,
	const std::string & fn,
	std::ostream & o
) {
	std::string line;
	QuerySet qs;
	std::ifstream in(fn);
	DIE_IF(!in.good(), "Couldn't open query file.");
	std::size_t n=0;
	while (getline(in, line)) {
		line = Trim(line);
		if (line.size() < jellyfish::mer_dna::k()) continue;
		qs.emplace_back(new QueryInfo(root->bf(), line));
		n++;
	}
	in.close();
	std::cerr << "Read " << n << " queries." << std::endl;

	// batch process the queries on ONLY the leaves
	query_leaves(root, qs);
	print_query_results(qs, o);
	
	for (auto & p : qs) {
		delete p;
	}
}

// This assumes num_hashes is one and stores the direct size_t position.
// **More accurate to say it scales with num_hashes but is inefficient if num_hashes > 1
// The alternative would be to store two size_t values (base and inc) [more efficient than many size_t]
std::vector<size_t> vector_hash_conversion(BF* root, std::vector<jellyfish::mer_dna> v){
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
}


// ********************* New code for split query parsing *********************
//
//
template<typename QueryStruct>

void split(const std::string &s, char delim, QueryStruct result) {
    std::stringstream ss(s);
    std::string item;
    while (std::getline(ss, item, delim)) {
        *(result++) = item;
    }
}

std::vector<std::string> split(const std::string &s, char delim) {
    std::vector<std::string> elems;
    split(s, delim, std::back_inserter(elems));
    return elems;
}


std::vector<std::vector<std::string> > splitQuery(const std::string &s){
    auto k = jellyfish::mer_dna::k();
    
    std::stringstream ss(s);
    std::stringstream subs;
    char c;

    std::vector<std::vector<std::string> > outVect;
    std::back_insert_iterator<std::vector<std::vector<std::string> > > temp = std::back_inserter(outVect);
    bool subs_b = false;
    char delim = ',';

    while (ss >> c){
        if (c == ']'){
            subs_b=false;
            *(temp++) = split(subs.str(), delim);
            subs.str("");
        }

        if (subs_b){ subs << c; }
        if (c == '['){ subs_b=true; }

    }

    // ** Add k-1 nucleotides to the end of each block except for last block **
    for(int i = 0; i < outVect.size() - 1; i++){
        std::string substring = outVect[i+1][0].substr(0,19);
        outVect[i][0].append(substring);
    }

    /*
    for (std::vector<std::vector<std::string> >::iterator it1 = outVect.begin(); it1 != outVect.end(); it1++){
        for (std::vector<std::string>::iterator it2 = (*it1).begin(); it2 != (*it1).end(); it2++){
            std::cerr << *it2 << ", " ;
        }
        std::cerr << std::endl;
    }
    */
    return outVect;
}

