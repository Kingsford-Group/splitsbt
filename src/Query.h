#ifndef QUERY_H
#define QUERY_H

#include <set>
#include <string>
#include <vector>
#include <list>
#include <iostream>
#include <iterator>

#include "BloomTree.h"
#include "SplitBloomTree.h"

extern float QUERY_THRESHOLD;

std::vector<size_t> vector_hash_conversion(BF* root, std::vector<jellyfish::mer_dna> v);

struct QueryInfo {
    QueryInfo(){}

    QueryInfo(BF* root, const std::string & q, const float qt){
        query = q;
        query_kmers = vector_hash_conversion(root, vector_kmers_in_string(q));
        total_kmers = query_kmers.size();
        tail_index=total_kmers-1;
        matched_kmers=0;
        q_thresh=qt;
    }

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
    float q_thresh = QUERY_THRESHOLD;
};


using QuerySet = std::list<QueryInfo*>;

void query_from_file(BloomTree* root, const std::string & fn, std::ostream & o);
void batch_query_from_file(BloomTree* root, const std::string & fn, std::ostream & o);
void batch_weightedquery_from_file(BloomTree* root, const std::string & fn, const std::string & wf, std::ostream & o); 
void query_string(BloomTree* root, const std::string & q, std::vector<BloomTree*> & out);
void query(BloomTree* root, const std::set<jellyfish::mer_dna> & q, std::vector<BloomTree*> & out);
void check_bt(BloomTree* root);
void draw_bt(BloomTree* root, std::string outfile);
void popcount_bt(BloomTree* root);
void compress_bt(BloomTree* root);
void compress_splitbt(BloomTree* root, sdsl::bit_vector* noninfo);
void compress_splitbt(BloomTree* root, BF* rbf); 
void leaf_query_from_file(BloomTree* root, const std::string & fn, std::ostream & o);

void batch_splitquery_from_file(BloomTree* root, const std::string & fn, std::ostream & o);

std::vector<std::vector<std::string> > splitQuery(const std::string &s);

// type == 0 ; contains
// type == 1 ; does not contain
struct splitQueryInfo {
    splitQueryInfo(BF* root, const std::string & sq, int t){
        type = t;

        std::vector<std::vector<std::string> > parsed_sq = splitQuery(sq);
        int count = 0;
        for (auto & psq : parsed_sq){
            const std::string temps = psq[0];
            const float tempf = std::stof(psq[1]);
            QueryInfo* tempQI = new QueryInfo(root, temps, tempf);
            partial_queries.emplace_back(tempQI);
            // std::copy(tempQI->query_kmers.begin(), tempQI->query_kmers.end(), std::inserter(total_kmers, total_kmers.end()));
            count++;
        } 
    }
    //std::set<size_t> total_kmers; // disabled until used
    QuerySet partial_queries;
    //std::vector<const BloomTree*> matching;
    //std::vector<QueryInfo> partial_queries;
    int type;
};


using splitQuerySet = std::list<splitQueryInfo*>;

struct batchInfo{
    batchInfo(BF* root, const std::string & line){
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
                sqi_list.emplace_back(new splitQueryInfo(root, line.substr(left_end, right_end-left_end+1),next_type));
                next_type = 0;
            } else if (c == '-'){
                left_end = right_end;
                right_end = index;
                sqi_list.emplace_back(new splitQueryInfo(root, line.substr(left_end, right_end-left_end+1),next_type));
                next_type = 1;
            }
            index++;
        }

        left_end = right_end;
        right_end = index;
        sqi_list.emplace_back(new splitQueryInfo(root, line.substr(left_end, right_end-left_end+1),next_type));
    }

    splitQuerySet sqi_list;
    std::vector<const BloomTree*> matching;
};

using batchSet = std::list<batchInfo*>;

#endif
