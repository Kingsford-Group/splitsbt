#ifndef QUERY_H
#define QUERY_H

#include <set>
#include <string>
#include <vector>
#include <list>
#include <iostream>

#include "BloomTree.h"
#include "SplitBloomTree.h"

extern float QUERY_THRESHOLD;

std::vector<size_t> vector_hash_conversion(BF* root, std::vector<jellyfish::mer_dna> v);

struct QueryInfo {
    QueryInfo(BF* root, const std::string & q) : 
    query(q), 
    query_kmers(vector_hash_conversion(root, vector_kmers_in_string(q))),
    total_kmers(query_kmers.size()),
    tail_index(total_kmers-1),
    matched_kmers(0)  {}

    QueryInfo(BF* root, const std::string & q, const std::string & w){
	    query = q;
    	query_kmers = vector_hash_conversion(root, vector_kmers_in_string(q));
        total_kmers = query_kmers.size();
        tail_index = total_kmers-1;
        matched_kmers = 0;
	    std::vector<std::string> fields;
    	SplitString(w, ' ', fields);
    	unsigned n = 0;
	    for (const auto & w : fields){
		    if (w!="") {
			    // If string w has invalid letters after numbers, composite_string will return pointer (else null)
    			// std::string::size_type not working here?
	    		std::size_t composite_string;
		    	//std::cerr << typeid(w).name() << std::endl;
			    try{
				    float value = std::stof(w, &composite_string);
	                if (composite_string != w.size()) {
                        std::cerr << "Invalid weight \'" << w << "\' at position " << n << std::endl;
        	            exit(3);
                    }
	                weight.emplace_back(value); //Currently zero error handling here!
		    	}
			    catch(...){
				    std::cerr << "Invalid weight \'" << w << "\' at position " << n << std::endl;
    				exit(3);
	    		}
	    	}
	    	n++;
    	}
    }

/*    
    QueryInfo(QueryInfo & copy) :
    query(copy.query),
    query_kmers(copy.query_kmers),
    total_kmers(copy.total_kmers),
    matched_kmers(copy.matched_kmers){}
*/
    ~QueryInfo() {}
   
    std::string query;
    std::vector<size_t> query_kmers;
    //std::vector<jellyfish::mer_dna> query_kmers;
    std::vector<const BloomTree*> matching;
    std::vector<float> weight;
    int total_kmers;
    int tail_index;
    int matched_kmers;
};

using QuerySet = std::list<QueryInfo*>;

void query_from_file(BloomTree* root, const std::string & fn, std::ostream & o);
void batch_query_from_file(BloomTree* root, const std::string & fn, std::ostream & o);
void batch_weightedquery_from_file(BloomTree* root, const std::string & fn, const std::string & wf, std::ostream & o); 
void query_string(BloomTree* root, const std::string & q, std::vector<BloomTree*> & out);
void query(BloomTree* root, const std::set<jellyfish::mer_dna> & q, std::vector<BloomTree*> & out);
void check_bt(BloomTree* root);
void draw_bt(BloomTree* root, std::string outfile);
void compress_bt(BloomTree* root);
void compress_splitbt(BloomTree* root, BF* rbf); 
void leaf_query_from_file(BloomTree* root, const std::string & fn, std::ostream & o);

#endif