#ifndef QUERY_H
#define QUERY_H

#include "boost/spirit/include/qi.hpp"
#include "boost/spirit/include/phoenix.hpp"
#include "boost/spirit/include/phoenix_operator.hpp"
#include "boost/variant/recursive_wrapper.hpp"


#include <set>
#include <string>
#include <vector>
#include <list>
#include <iostream>
#include <iterator>

#include "BloomTree.h"
#include "SplitBloomTree.h"
#include "BF.h"

extern float QUERY_THRESHOLD;

std::vector<size_t> vector_hash_conversion(BF* root, std::vector<jellyfish::mer_dna> v);

/***** PART OF NEW QUERY STUFF***/
namespace qi    = boost::spirit::qi;
namespace phx   = boost::phoenix;

struct op_or  {};
struct op_and {};
struct op_xor {};
struct op_not {};

typedef std::string var;
template <typename tag> struct binop;
template <typename tag> struct unop;

typedef boost::variant<var,
        boost::recursive_wrapper<unop <op_not> >,
        boost::recursive_wrapper<binop<op_and> >,
        boost::recursive_wrapper<binop<op_xor> >,
        boost::recursive_wrapper<binop<op_or> >
        > expr;

template <typename tag> struct binop
{
    explicit binop(const expr& l, const expr& r) : oper1(l), oper2(r) { }
    expr oper1, oper2;
};

template <typename tag> struct unop
{
    explicit unop(const expr& o) : oper1(o) { }
    expr oper1;
};

template <typename It, typename Skipper = qi::space_type>
    struct parser : qi::grammar<It, expr(), Skipper>
{
    parser() : parser::base_type(expr_)
    {
        using namespace qi;

        expr_  = or_.alias();

        or_  = (xor_ >> "|"  >> or_ ) [ _val = phx::construct<binop<op_or >>(_1, _2) ] | xor_   [ _val = _1 ];
        xor_ = (and_ >> "^" >> xor_) [ _val = phx::construct<binop<op_xor>>(_1, _2) ] | and_   [ _val = _1 ];
        and_ = (not_ >> "&" >> and_) [ _val = phx::construct<binop<op_and>>(_1, _2) ] | not_   [ _val = _1 ];
        not_ = ("!" > simple       ) [ _val = phx::construct<unop <op_not>>(_1)     ] | simple [ _val = _1 ];

        simple = (('(' > expr_ > ')') | var_);
        var_ = qi::lexeme[ +alpha | +digit ];

        BOOST_SPIRIT_DEBUG_NODE(expr_);
        BOOST_SPIRIT_DEBUG_NODE(or_);
        BOOST_SPIRIT_DEBUG_NODE(xor_);
        BOOST_SPIRIT_DEBUG_NODE(and_);
        BOOST_SPIRIT_DEBUG_NODE(not_);
        BOOST_SPIRIT_DEBUG_NODE(simple);
        BOOST_SPIRIT_DEBUG_NODE(var_);
    }

  private:
    qi::rule<It, var() , Skipper> var_;
    qi::rule<It, expr(), Skipper> not_, and_, xor_, or_, simple, expr_;
};
/**** END OF NEW QUERY STUFF**/

struct QueryInfo {
    QueryInfo(){}

    QueryInfo(BF* root, const std::string & q, const float qt){
        query = q;
        query_kmers = vector_hash_conversion(root, vector_kmers_in_string(q));
        total_kmers = query_kmers.size();
        tail_index=total_kmers-1;
        matched_kmers=0;
        local_kmers=0;
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
    int local_kmers;
    float q_thresh = QUERY_THRESHOLD;
};


using QuerySet = std::vector<QueryInfo*>;

void query_from_file(BloomTree* root, const std::string & fn, std::ostream & o);
void batch_query_from_file(BloomTree* root, const std::string & fn, std::ostream & o);
void batch_weightedquery_from_file(BloomTree* root, const std::string & fn, const std::string & wf, std::ostream & o); 
void query_string(BloomTree* root, const std::string & q, std::vector<BloomTree*> & out);
void query(BloomTree* root, const std::set<jellyfish::mer_dna> & q, std::vector<BloomTree*> & out);
void check_bt(BloomTree* root);
void draw_bt(BloomTree* root, std::string outfile);
void popcount_bt(BloomTree* root);
void compress_bt(BloomTree* root);
void compress_splitbt_helper(BloomTree* root, sdsl::bit_vector* noninfo);
void compress_splitbt(BloomTree* root, bool new_only);
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

using kmer_type = size_t;
// Stores pointers to querySet objects for each kmer
// If tail_index is -1, no queries are using this kmer
struct queryVec{
    queryVec(){
        tail_index=-1;
    }
    queryVec(QueryInfo* matchQ){
        mySet.push_back(matchQ);
        tail_index++;// = 0;//mySet.size();
    }
    int tail_index=-1;
    std::vector<QueryInfo*> mySet;
};

//stores the original value of the 'kmer'
//and the current value to map back to queries
struct splitKmer{
    splitKmer(){}

    splitKmer(kmer_type val){
        curr = val;
        orig = val;
    }
    kmer_type curr;
    kmer_type orig;
};

// SUPPORT FUNCTIONS FOR NEW QUERY (USING BOOST)
bool local_match(QuerySet& qs, const var& v);
bool global_match(QuerySet& qs, const var& v);
bool local_query(QuerySet& qs, const expr& e);
bool global_query(QuerySet& qs, const expr& e);

struct boolInfo{
    boolInfo(BF* root, const std::string & line){
        std::stringstream ss(line);
        std::stringstream subs;
        std::stringstream editted;
        char c;

        bool subs_b = false;
        char delim = ',';
        int count = 0;

        while (ss >> c){
            if (c == '['){ 
                subs_b=true; 
            }
            else if (c == ']'){
                subs_b=false;

                std::vector<std::string> parse;
                int args = SplitString(subs.str(), delim, parse);
                if(args != 2){ DIE("Query should be of type [TCGA,0.9]"); }
                part_q.emplace_back(new QueryInfo(root, parse[0], std::stof(parse[1])));
                editted << "(" << count << ")";
                count++;
                subs.str("");
            } else if (subs_b){ 
                subs << c; 
            }
            else{ editted << c; }

        }
       
        std::string store_string = editted.str();
        std::cerr << "Read " << count << " queries in bool string." << std::endl;
        std::cerr << "Processed boolean query: " << store_string << std::endl;
        auto f(std::begin(store_string)), l(std::end(store_string));
        parser<decltype(f)> p;

        try
        {
            bool ok = qi::phrase_parse(f,l,p > ';',qi::space,bool_expr);

            if (!ok){ DIE("boolean expression invalid!"); }
        } catch(const qi::expectation_failure<decltype(f)>& e)
        {
            std::cerr << "expectation_failure at '" << std::string(e.first, e.last) << "'\n";
            DIE("Canceling query.");
        }
            

    }
    QuerySet part_q;
    expr bool_expr;
    std::vector<const BloomTree*> matching;

    bool localPasses()
    {
        return local_query(part_q, bool_expr);
    } 
    bool globalPasses()
    {
        return global_query(part_q, bool_expr);
    }
};

using boolQuerySet = std::list<boolInfo*>;

struct boolQuery{
    boolQuery(BloomTree* root, const std::string & queryFile){
        std::string line;
        std::ifstream in(queryFile);
        DIE_IF(!in.good(), "Couldn't open query file.");

        std::size_t n = 0;
        while (getline(in, line)) {
            line = Trim(line);
            if (line.size() < jellyfish::mer_dna::k()) continue;
            bs.emplace_back(new boolInfo(root->bf(), line));
            n++;
        }
        in.close();
        std::cerr << "Read " << n << " queries." << std::endl;
    
        std::set<std::size_t> temp;
        for (auto bi : bs){
            for (auto q : bi->part_q){
                qs.emplace_back(q);
                for (auto km : q->query_kmers){
                    //add_query_to_splitkmer(kmer_set[km],q);
                    query_vec[km].mySet.push_back(q);
                    query_vec[km].tail_index++;
                    temp.insert(km);
                }
            }
        }
        std::cerr << "Built query map with " << query_vec.size() << " elements" << std::endl;

        //Convert set to vector
        std::vector<size_t> temp_v(temp.size());
        std::copy(temp.begin(), temp.end(), temp_v.begin());

        //Build splitKmer vector
        std::vector<splitKmer> temp_sk(temp.size());
        batch_kmers = temp_sk;
        for(int i =0; i< temp_v.size(); i++){
            batch_kmers[i]=splitKmer(temp_v[i]);
        }

        std::cerr << "Built kmer vector with " << batch_kmers.size() << " elements" << std::endl;
        
        assert(query_vec.size() == batch_kmers.size());
        tail_index = query_vec.size()-1;
        std::vector<bool> myVec(query_vec.size());
        hit_vector = myVec;
    }
    QuerySet qs;
    std::unordered_map<kmer_type,queryVec> query_vec;
    std::vector<splitKmer> batch_kmers; // the equivalent of query_kmers
    int tail_index; 
    std::vector<bool> hit_vector;
    boolQuerySet bs;
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


struct batchQuery{
    batchQuery(BloomTree* root, const std::string & fn){
        std::string line;
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

        std::set<std::size_t> temp;
        for (auto q : qs){
            for (auto km : q->query_kmers){
                //add_query_to_splitkmer(kmer_set[km],q);
                query_vec[km].mySet.push_back(q);
                query_vec[km].tail_index++;
                temp.insert(km);
            }
        }
        std::cerr << "Built query map with " << query_vec.size() << " elements" << std::endl;

        std::vector<size_t> temp_v(temp.size());
        std::copy(temp.begin(), temp.end(), temp_v.begin());

        std::vector<splitKmer> temp_sk(temp.size());
        batch_kmers = temp_sk;
        for(int i =0; i< temp_v.size(); i++){
            batch_kmers[i]=splitKmer(temp_v[i]);
        }

        std::cerr << "Built kmer vector with " << batch_kmers.size() << " elements" << std::endl;
        
        assert(query_vec.size() == batch_kmers.size());
        tail_index = query_vec.size()-1;
        std::vector<bool> myVec(query_vec.size());
        hit_vector = myVec;
    }
    QuerySet qs;
    std::unordered_map<kmer_type,queryVec> query_vec;
    std::vector<splitKmer> batch_kmers; // the equivalent of query_kmers
    int tail_index; 
    std::vector<bool> hit_vector;
    
};

void split_query_batch(BloomTree* root, batchQuery & bq);
void split_noquery_batch(BloomTree* root, batchQuery & bq);
void print_batch_query(const batchQuery & bq, std::ostream & out);
void bool_query_batch(BloomTree* root, boolQuery & bq);
void bool_noquery_batch(BloomTree* root, boolQuery & bq);
void print_bool_query(const boolQuery & bq, std::ostream & out);
void batchkmer_swap(std::vector<splitKmer>& v, int opos, int npos);

/******************** Code for binary tree boolean logic parser *************************/



struct localMatch : boost::static_visitor<bool>{
    localMatch(QuerySet& qs) : _qs(qs) { }
    QuerySet& _qs;

    bool operator()(const var& v) const { return local_match(_qs, v); }
    bool operator()(const binop<op_and>& b) const { return process(0, b.oper1, b.oper2); }
    bool operator()(const binop<op_or >& b) const { return process(1, b.oper1, b.oper2); }
    bool operator()(const binop<op_xor>& b) const { return process(2, b.oper1, b.oper2); }

    bool process(const int& type, const expr& l, const expr& r) const{
        bool lproc = boost::apply_visitor(*this, l);
        bool rproc = boost::apply_visitor(*this, r);
        if(type==0){ return (lproc && rproc); }
        if(type==1){ return (lproc || rproc); }
        if(type==2){ return (lproc != rproc); }

        std::cerr <<"Only 3 supported binops (and, or, xor)" << std::endl;
        return false;
    }

    bool operator()(const unop<op_not>& u) const{
        bool v = boost::apply_visitor(*this, u.oper1);
        return !v;
    }
};

struct globalMatch : boost::static_visitor<bool>{
    globalMatch(QuerySet& qs) : _qs(qs) { }
    QuerySet& _qs;

    bool operator()(const var& v) const { return global_match(_qs, v); }
    bool operator()(const binop<op_and>& b) const { return process(0, b.oper1, b.oper2); }
    bool operator()(const binop<op_or >& b) const { return process(1, b.oper1, b.oper2); }
    bool operator()(const binop<op_xor>& b) const { return process(2, b.oper1, b.oper2); }

    bool process(const int& type, const expr& l, const expr& r) const{
        bool lproc = boost::apply_visitor(*this, l);
        bool rproc = boost::apply_visitor(*this, r);
        if(type==0){ return (lproc && rproc); }
        if(type==1){ return (lproc || rproc); }
        if(type==2){ return (lproc != rproc); }

        std::cerr <<"Only 3 supported binops (and, or, xor)" << std::endl;
        return false;
    }

    bool operator()(const unop<op_not>& u) const{
        bool v = boost::apply_visitor(*this, u.oper1);
        return !v;
    }
};



#endif
