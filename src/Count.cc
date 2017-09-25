#include "Count.h"
#include "Kmers.h"
#include "BF.h"
#include <jellyfish/mer_dna_bloom_counter.hpp>
#include <jellyfish/file_header.hpp>
#include <sdsl/bit_vectors.hpp>

#include <jellyfish/mer_dna.hpp>
#include <jellyfish/thread_exec.hpp>
#include <jellyfish/hash_counter.hpp>
#include <jellyfish/stream_manager.hpp>
#include <jellyfish/mer_overlap_sequence_parser.hpp>
#include <jellyfish/mer_iterator.hpp>

#include <vector>

/*==== COPIED FROM THE JF count_dump example ====*/

typedef jellyfish::cooperative::hash_counter<jellyfish::mer_dna> mer_hash_type;
typedef jellyfish::stream_manager<std::vector<std::string>::iterator > stream_manager_type;
typedef jellyfish::mer_overlap_sequence_parser<stream_manager_type> sequence_parser_type;
typedef jellyfish::mer_iterator<sequence_parser_type, jellyfish::mer_dna> mer_iterator_type;


class mer_counter : public jellyfish::thread_exec {
  mer_hash_type&       mer_hash_;
  sequence_parser_type parser_;
  const bool           canonical_;
  std::mutex           lock_;

public:
  mer_counter(int nb_threads, mer_hash_type& mer_hash,
              stream_manager_type  streams,
              bool canonical)
    : mer_hash_(mer_hash)
    , parser_(jellyfish::mer_dna::k(), streams.nb_streams(), 3 * nb_threads, 4096, streams)
    , canonical_(canonical)
  { }

  virtual void start(int thid) {
    mer_iterator_type mers(parser_, canonical_);
    size_t count = 0;
    for( ; mers; ++mers) {
      mer_hash_.add(*mers, 1);
        ++count;
    }
    mer_hash_.done();
    {
    std::lock_guard<std::mutex> lck(lock_);
    std::cerr << "thread " << thid << " added " << count << std::endl;
    }
  }
};

/*=== END COPY ===*/


// run JF count, and build a BF and save it to disk

enum OPERATION { COUNT, PRIME, UPDATE };

bool count(
    std::string infilen,
    std::string outfilen,
    HashPair hp,
    int nh,
    uint64_t bf_size,
    int num_threads,
    unsigned cutoff_count
    ) {
   
#ifdef HAVE_HTSLIB
    std::cerr << "HAVE_HTSLIB FLAG: TRUE" << std::endl;
#endif
 
    // jellyfish default counting values
    const uint64_t hash_size = 10000000;
    const uint32_t num_reprobes = 126;
    const uint32_t counter_len = 7;
    const bool canonical = true;

    // create the hash
    mer_hash_type mer_hash(hash_size, jellyfish::mer_dna::k()*2, counter_len, num_threads, num_reprobes);

    // create a mock up of the array of file names
    std::vector<std::string> files;
    files.push_back(infilen);

    // count the kmers
    int count_type=0;
    stream_manager_type streams;
    if (infilen.substr(infilen.size()-4) == ".bam"){
        std::cerr << "Processing bam file." << std::endl;
        streams.sams(files.begin(), files.end());
    } else {
      std::cerr << "Processing fast[aq] file." << std::endl;
      streams.paths(files.begin(), files.end());
    }
    mer_counter counter(num_threads, mer_hash, streams, canonical);
    counter.exec_join(num_threads);

    // build the BF
    UncompressedBF bf(outfilen, hp, nh, bf_size);

    // add each kmer to the BF
    const auto jf_ary = mer_hash.ary();
    const auto end = jf_ary->end();
    std::cerr << "Right before cutoff count: " << cutoff_count << std::endl;
    size_t in_hash = 0, in_bloom = 0;
    for(auto kmer = jf_ary->begin(); kmer != end; ++kmer) {
        ++in_hash;
        auto& key_val = *kmer;
        if (key_val.second >= cutoff_count) {
            bf.add(key_val.first);
            ++in_bloom;
        }
    }
    std::cerr << "in hash: " << in_hash << " in blom: " << in_bloom << std::endl;
    bf.save();
    return true;
}

bool split_count(
    std::string infilen,
    std::string outfilen,
    HashPair hp,
    int nh,
    uint64_t bf_size,
    int num_threads,
    unsigned cutoff_count
    ) {

    // jellyfish default counting values
    const uint64_t hash_size = 10000000;
    const uint32_t num_reprobes = 126;
    const uint32_t counter_len = 7;
    const bool canonical = true;

    // create the hash
    mer_hash_type mer_hash(hash_size, jellyfish::mer_dna::k()*2, counter_len, num_threads, num_reprobes);

    // create a mock up of the array of file names
    std::vector<std::string> files;
    files.push_back(infilen);

    // count the kmers
    stream_manager_type streams(files.begin(), files.end());
    mer_counter counter(num_threads, mer_hash, streams, canonical);
    counter.exec_join(num_threads);

    // build the BF
    SBF bf(outfilen, hp, nh, bf_size);

    // add each kmer to the BF
    const auto jf_ary = mer_hash.ary();
    const auto end = jf_ary->end();
    std::cerr << "Right before cutoff count: " << cutoff_count << std::endl;
    for(auto kmer = jf_ary->begin(); kmer != end; ++kmer) {
        auto& key_val = *kmer;
        if (key_val.second >= cutoff_count) {
            bf.add(key_val.first);
        }
    }
    bf.save();
    return true;
}
// bt count k FASTA OUTFILE

 /*
    // create the hash
    jellyfish::mer_hash ary(hash_size, k*2, counter_len, num_threads, num_reprobes);

    // prime the kmer values
    jellyfish::mer_counter primer(num_threads, ary,
            files.begin(), files.end(),
            files.end(), files.end(), 1, PRIME); 
    primer.exec_join(num_threads);

    // count the kmers
    jellyfish::mer_counter counter(num_threads, ary,
            files.begin(), files.end(),
            files.end(), files.end(),
            1, UPDATE);
    counter.exec_join(num_threads);
    */
