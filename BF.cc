#include "BF.h"
#include "Kmers.h"
#include "util.h"

#include <jellyfish/file_header.hpp>

BF::BF(const std::string & f, HashPair hp, int nh) :
    filename(f),
    bits(nullptr),
    hashes(hp),
    num_hash(nh)
{ 
}

BF::BF(const std::string & f, const std::string & f_2, HashPair hp, int nh) :
    filename(f),
    f2(f_2),
    bits(nullptr),
    hashes(hp),
    num_hash(nh)
{
}

BF::~BF() {
    if (bits != nullptr) {
        delete bits;
    }
}

void BF::add(const jellyfish::mer_dna & m) {
    jellyfish::mer_dna can(m);
    can.canonicalize();
    uint64_t h0 = hashes.m1.times(can);
    uint64_t h1 = hashes.m2.times(can);
    
    const size_t base = h0 % size();
    const size_t inc = h1 % size();

    for (unsigned long i = 0; i < num_hash; ++i) {
        const size_t pos = (base + i * inc) % size();
        this->set_bit(pos);
    }
}

// set the bit at position 1 (can't use operator[] b/c we'd need
// to return a proxy, which is more trouble than its worth)
void BF::set_bit(uint64_t p) {
    DIE("Compressed BF are not mutable!");
}

// returns true iff the bloom filter contains the given kmer
bool BF::contains(const jellyfish::mer_dna & m) const {
    //std::cout << "TESTING STUFF: " << m.to_str() << std::endl;
    std::string temp = m.to_str();
    jellyfish::mer_dna n = jellyfish::mer_dna(temp);
    n.canonicalize();
    //std::cout << "Canonical version! " << n.to_str() << std::endl;
    uint64_t h0 = hashes.m1.times(n);
    uint64_t h1 = hashes.m2.times(n);

    //DEBUG: std::cout << "size = " << bits->size() << std::endl;
    
    const size_t base = h0 % size();
    const size_t inc = h1 % size();

    for (unsigned long i = 0; i < num_hash; ++i) {
        const size_t pos = (base + i * inc) % size();
        //DEBUG: std::cout << "pos=" << pos << std::endl;
        if ((*this)[pos] == 0) return false;
    }
    return true;
}

// convience function
bool BF::contains(const std::string & str) const {
    //jellyfish::mer_dna temp = jellyfish::mer_dna(str);
    //temp.canonicalize();
    //std::cout << "checking: " << str << "\n";
    //std::cout << "canonical: " << temp.to_str() << "\n";
    return contains(jellyfish::mer_dna(str));
}

uint64_t BF::size() const {
    return bits->size();
}

int BF::operator[](uint64_t pos) const {
    return (*bits)[pos];
}


// read the bit vector and the matrices for the hash functions.
void BF::load() {
    // read the actual bits
    bits = new sdsl::rrr_vector<255>();
    sdsl::load_from_file(*bits, filename);
}

void BF::save() {
    std::cerr << "Saving BF to " << filename << std::endl;
    sdsl::store_to_file(*bits, filename);
}


// create a new RRR bloom filter that is the union of this BF and the given BF.
// Will re-use the hashes from this and both BFs must use exactly the same hash
// (not checked).
BF* BF::union_with(const std::string & new_name, const BF* f2) const {
    assert(bits->size() == f2->size());

    // create an uncompressed version of this bf
    sdsl::bit_vector b(bits->size(), 0);

    std::cerr << "Performing OR... (size " << b.size() << ")" << std::endl;

    // union it with f2
    for (unsigned long i = 0; i < b.size(); i++) {
        b[i] = (*bits)[i] | (*f2->bits)[i];
        if (i % 1000000 == 0) std::cerr << "i=" << i << std::endl;
    }

    // create a new BF wraper for the new BF
    std::cerr << "Building BF object..." << std::endl;
    BF* out = new BF(new_name, hashes, num_hash);

    // "load" the new BF by converting it to a RRR
    std::cerr << "Building RRR vector..." << std::endl;
    out->bits = new sdsl::rrr_vector<255>(b);
    return out;
}

uint64_t BF::similarity(const BF* other, int type) const {
    DIE("not yet implemented");
    return 0;
}

std::tuple<uint64_t, uint64_t> BF::b_similarity(const BF* other) const {
    DIE("not yet implemented");
    return std::make_tuple(0,0);
}


void BF::union_into(const BF* other) {
    DIE("not yet implemented");
}

uint64_t BF::count_ones() const {
    DIE("not yet implemented");
    return 0;
}

void BF::compress() {
	DIE("Cant compress rrr further with existing code base");
}

BF* BF::sim_with(const std::string & new_name, const BF* f2) const{
    assert(bits->size() == f2->size());

    // create an uncompressed version of this bf
    sdsl::bit_vector b(bits->size(), 0);

    std::cerr << "Performing AND... (size " << b.size() << ")" << std::endl;

    // union it with f2
    for (unsigned long i = 0; i < b.size(); i++) {
        b[i] = (*bits)[i] & (*f2->bits)[i];
        if (i % 1000000 == 0) std::cerr << "i=" << i << std::endl;
    }

    // create a new BF wraper for the new BF
    std::cerr << "Building BF object..." << std::endl;
    BF* out = new BF(new_name, hashes, num_hash);

    // "load" the new BF by converting it to a RRR
    std::cerr << "Building RRR vector..." << std::endl;
    out->bits = new sdsl::rrr_vector<255>(b);
    return out;
}
void BF::sim_into(const BF* f2){
    DIE("not yet implemented");
}

BF* BF::dif_with(const std::string & new_name, const BF* f2) const{
    assert(bits->size() == f2->size());

    // create an uncompressed version of this bf
    sdsl::bit_vector b(bits->size(), 0);

    std::cerr << "Performing XOR... (size " << b.size() << ")" << std::endl;

    // union it with f2
    for (unsigned long i = 0; i < b.size(); i++) {
        b[i] = (*bits)[i] ^ (*f2->bits)[i];
        if (i % 1000000 == 0) std::cerr << "i=" << i << std::endl;
    }

    // create a new BF wraper for the new BF
    std::cerr << "Building BF object..." << std::endl;
    BF* out = new BF(new_name, hashes, num_hash);

    // "load" the new BF by converting it to a RRR
    std::cerr << "Building RRR vector..." << std::endl;
    out->bits = new sdsl::rrr_vector<255>(b);
    return out;
}
void BF::dif_into(const BF* f2){
    DIE("not yet implemented");
}

/*
Uncompressed Bloom Filter (SBT 1.0 Implementation)
*/
/*============================================*/

UncompressedBF::UncompressedBF(const std::string & f, HashPair hp, int nh, uint64_t size) :
    BF(f, hp, nh),
    bv(nullptr)
{
    if (size > 0) {
        bv = new sdsl::bit_vector(size);
    }
}


UncompressedBF::~UncompressedBF() {
    // call to base destructor happens automatically
    if (bv != nullptr) {
        delete bv;
    }
}

void UncompressedBF::load() {
    // read the actual bits
    assert(bv == nullptr);
    bv = new sdsl::bit_vector();
    sdsl::load_from_file(*bv, filename);
}

void UncompressedBF::save() {
    std::cerr << "Saving BF to " << filename << std::endl;
    sdsl::store_to_file(*bv, filename);
}


uint64_t UncompressedBF::size() const {
    return bv->size();
}

int UncompressedBF::operator[](uint64_t pos) const {
    return (*bv)[pos];
}

void UncompressedBF::set_bit(uint64_t p) {
    (*bv)[p] = 1;
}

BF* UncompressedBF::union_with(const std::string & new_name, const BF* f2) const {
    assert(size() == f2->size());
    const UncompressedBF* b = dynamic_cast<const UncompressedBF*>(f2);
    if (b == nullptr) {
        DIE("Can only union two uncompressed BF");
    }
    UncompressedBF* out = new UncompressedBF(new_name, hashes, num_hash);
    out->bv = union_bv_fast(*this->bv, *b->bv);
    return out;
}

void UncompressedBF::union_into(const BF* f2) {
    assert(size() == f2->size());

    const UncompressedBF* b = dynamic_cast<const UncompressedBF*>(f2);
    if (b == nullptr) {
        DIE("Can only union two uncompressed BF");
    }

    uint64_t* b1_data = this->bv->data();
    const uint64_t* b2_data = b->bv->data();

    sdsl::bit_vector::size_type len = size()>>6;
    for (sdsl::bit_vector::size_type p = 0; p < len; ++p) {
        *b1_data = (*b1_data) | (*b2_data++);
        b1_data++;
    }
}

uint64_t UncompressedBF::similarity(const BF* other, int type) const {
    assert(other->size() == size());

    const uint64_t* b1_data = this->bv->data();
    const UncompressedBF* o = dynamic_cast<const UncompressedBF*>(other);
    if (o == nullptr) {
        DIE("Can only compute similarity on same type of BF.");
    }

    const uint64_t* b2_data = o->bv->data();
    	if (type == 1) {
		uint64_t xor_count = 0;
	   	uint64_t or_count = 0;
    
		sdsl::bit_vector::size_type len = size()>>6;
		for (sdsl::bit_vector::size_type p = 0; p < len; ++p) {
        		xor_count += __builtin_popcountl((*b1_data) ^ (*b2_data));
	        	or_count += __builtin_popcountl((*b1_data++) | (*b2_data++));
			//count += __builtin_popcountl((*b1_data++) ^ (*b2_data++));
  	  	}
   
    		return uint64_t(float(or_count - xor_count) / float(or_count) * size() );
	}
	else if (type == 0) {
		uint64_t count = 0;
		sdsl::bit_vector::size_type len = size()>>6;
                for (sdsl::bit_vector::size_type p = 0; p < len; ++p) {
                        count += __builtin_popcountl((*b1_data++) ^ (*b2_data++));
                }
		return size() - count;
	}
	
	DIE("ERROR: ONLY TWO TYPES IMPLEMENTED");
	return 0;
}


// 1-15-16 This was corrected from 'and_count' to 'xor_count'. The bit operators were unchanged
std::tuple<uint64_t, uint64_t> UncompressedBF::b_similarity(const BF* other) const {
    assert(other->size() == size());

    const uint64_t* b1_data = this->bv->data();
    const UncompressedBF* o = dynamic_cast<const UncompressedBF*>(other);
    if (o == nullptr) {
        DIE("Can only compute similarity on same type of BF.");
    }

    const uint64_t* b2_data = o->bv->data();

    uint64_t xor_count = 0;
    uint64_t or_count = 0;
    sdsl::bit_vector::size_type len = size()>>6;
    for (sdsl::bit_vector::size_type p = 0; p < len; ++p) {
        xor_count += __builtin_popcountl((*b1_data) ^ (*b2_data));
        or_count += __builtin_popcountl((*b1_data++) | (*b2_data++));
        //count += __builtin_popcountl((*b1_data++) ^ (*b2_data++));
    }
    return std::make_tuple(xor_count, or_count);
//size() - count;
}


// 00 = 0
// 01 = 1
// 10 = 1
// 11 = 0

int popcount(uint64_t b) {
    int count = 0;
    for (int i = 0; i < 64; i++) {
        if (b % 2) count++;
        b >>= 1;
    }
    return count;
}

// return the # of 1s in the bitvector
uint64_t UncompressedBF::count_ones() const {
    uint64_t* data = bv->data();
    sdsl::bit_vector::size_type len = bv->size()>>6;
    uint64_t count = 0;
    for (sdsl::bit_vector::size_type p = 0; p < len; ++p) {
        int pc = popcount(*data);
        DIE_IF(pc != __builtin_popcountl(*data), "popcountl and us disagree about popcount");
        data++;
        count += pc;
    }

    sdsl::rank_support_v<> rbv(bv);
    DIE_IF(rbv.rank(bv->size()) != count, "SDSL and us disagree about number of 1s");
    return count;
}

void UncompressedBF::compress() {
	sdsl::rrr_vector<255> rrr(*bv);
	std::cerr << "Compressed RRR vector is " << sdsl::size_in_mega_bytes(rrr) << std::endl;
	sdsl::store_to_file(rrr,filename+".rrr");
}


BF* UncompressedBF::sim_with(const std::string & new_name, const BF* f2) const{
    assert(size() == f2->size());
    const UncompressedBF* b = dynamic_cast<const UncompressedBF*>(f2);
    if (b == nullptr) {
        DIE("Can only union two uncompressed BF");
    }
    UncompressedBF* out = new UncompressedBF(new_name, hashes, num_hash);
    out->bv = sim_bv_fast(*this->bv, *b->bv);
    return out;
}


void UncompressedBF::sim_into(const BF* f2){
   assert(size() == f2->size());

    const UncompressedBF* b = dynamic_cast<const UncompressedBF*>(f2);
    if (b == nullptr) {
        DIE("Can only union two uncompressed BF");
    }

    uint64_t* b1_data = this->bv->data();
    const uint64_t* b2_data = b->bv->data();

    sdsl::bit_vector::size_type len = size()>>6;
    for (sdsl::bit_vector::size_type p = 0; p < len; ++p) {
        *b1_data = (*b1_data) & (*b2_data++);
        b1_data++;
    }
}

BF* UncompressedBF::dif_with(const std::string & new_name, const BF* f2) const{
    assert(size() == f2->size());
    const UncompressedBF* b = dynamic_cast<const UncompressedBF*>(f2);
    if (b == nullptr) {
        DIE("Can only union two uncompressed BF");
    }
    UncompressedBF* out = new UncompressedBF(new_name, hashes, num_hash);
    out->bv = dif_bv_fast(*this->bv, *b->bv);
    return out;
}

void UncompressedBF::dif_into(const BF* f2){
   assert(size() == f2->size());

    const UncompressedBF* b = dynamic_cast<const UncompressedBF*>(f2);
    if (b == nullptr) {
        DIE("Can only union two uncompressed BF");
    }

    uint64_t* b1_data = this->bv->data();
    const uint64_t* b2_data = b->bv->data();

    sdsl::bit_vector::size_type len = size()>>6;
    for (sdsl::bit_vector::size_type p = 0; p < len; ++p) {
        *b1_data = (*b1_data) ^ (*b2_data++);
        b1_data++;
    }
}

/*
Split Bloom Tree Compressed
*/
//====================================================================//
SBF::SBF(const std::string & f_sim, const std::string & f_dif, HashPair hp, int nh, uint64_t size) :
    BF(f_sim, f_dif, hp, nh),
    sim(nullptr),
    dif(nullptr)
{
    if (size > 0) {
        sim = new sdsl::bit_vector(size);
        dif = new sdsl::bit_vector(size);
    }
}


SBF::~SBF() {
    // call to base destructor happens automatically
    if (sim != nullptr) {
        delete sim;
    }

    if (dif != nullptr) {
        delete dif;
    }

}

//Before:
// split_filename() function added to address complications with one filename to two files.
//After:
//Now there's two filenames.
void SBF::load() {
    // read the actual bits
    assert(sim == nullptr);
    assert(dif == nullptr);
    sim = new sdsl::bit_vector();
    dif = new sdsl::bit_vector();
    sdsl::load_from_file(*sim, filename);
    sdsl::load_from_file(*dif, f2);
}

void SBF::save() {
    std::cerr << "Saving BF to " << filename << std::endl;
    sdsl::store_to_file(*sim, filename);
    sdsl::store_to_file(*dif, f2);
}


uint64_t SBF::size() const {
    uint64_t sim_size = sim->size();
    uint64_t dif_size = dif->size();
    assert(sim_size == dif_size);
    return sim_size;
}

int SBF::operator[](uint64_t pos) const {
    return (*sim)[pos] | (*dif)[pos];
}


//Set bit sets 'sim' filter
//This adheres to base BF add().
// XXX: MAKE SURE THIS IS ONLY USED IN LEAF BUILDING.
void SBF::set_bit(uint64_t p) {
    (*sim)[p] = 1;
}

void SBF::set_difbit(uint64_t p){
    (*dif)[p] = 1;
}

void SBF::unset_bit(uint64_t p) {
    (*sim)[p] = 0;
}

void SBF::unset_difbit(uint64_t p){
    (*dif)[p] = 0;
}


// When SBF unions, similar elements must be removed from both previous filters.
// That is going to be handled separately. 
BF* SBF::union_with(const std::string & new_sim, std::string & new_dif, const BF* f2) const {
    assert(size() == f2->size());
    const SBF* b = dynamic_cast<const SBF*>(f2);
    if (b == nullptr) {
        DIE("Can only union two uncompressed BF");
    }
    SBF* out = new SBF(new_sim, new_dif, hashes, num_hash);
    out->sim = sim_bv_fast(*this->sim, *b->sim); //sim filter is sim of two sim filters
    sdsl::bit_vector* new_diff = dif_bv_fast(*this->sim, *b->sim); //Dif filter is union of dif filters and new differences
    sdsl::bit_vector* temp = union_bv_fast(*this->dif, *b->dif); //Dif filters can be directly unioned
    out->dif = union_bv_fast(*temp, *new_diff);
    return out;
}

void SBF::union_into(const BF* f2) {
    assert(size() == f2->size());

    const SBF* b = dynamic_cast<const SBF*>(f2);
    if (b == nullptr) {
        DIE("Can only union two uncompressed BF");
    }

    uint64_t* b1_sim_data = this->sim->data();
    uint64_t* b1_dif_data = this->dif->data();
    const uint64_t* b2_sim_data = b->sim->data();
    const uint64_t* b2_dif_data = b->dif->data();

    sdsl::bit_vector::size_type len = size()>>6;
    for (sdsl::bit_vector::size_type p = 0; p < len; ++p) {
        *b1_sim_data = (*b1_sim_data) & (*b2_sim_data);
        *b1_dif_data = ((*b1_dif_data) | (*b2_dif_data)) | ((*b1_sim_data) ^ (*b2_sim_data));
        b1_sim_data++;
        b1_dif_data++;
    }
}

// Similarity could mean something different for SBF versus BF but right now it doesnt.
// This function simply is adapted to compare xor and or by combining dif and sim
uint64_t SBF::similarity(const BF* other, int type) const {
    assert(other->size() == size());

    const uint64_t* b1_sim_data = this->sim->data();
    const uint64_t* b1_dif_data = this->dif->data();

    const SBF* o = dynamic_cast<const SBF*>(other);
    if (o == nullptr) {
        DIE("Can only compute similarity on same type of BF.");
    }

    const uint64_t* b2_sim_data = o->sim->data();
    const uint64_t* b2_dif_data = o->dif->data();

    if (type == 1) {
		uint64_t xor_count = 0;
	   	uint64_t or_count = 0;
    
		sdsl::bit_vector::size_type len = size()>>6;
		for (sdsl::bit_vector::size_type p = 0; p < len; ++p) {
        		xor_count += __builtin_popcountl( ((*b1_sim_data) ^ (*b2_sim_data)) | ((*b1_dif_data) ^ (*b2_dif_data)) );
	        	or_count += __builtin_popcountl( ((*b1_sim_data) | (*b2_sim_data)) | ((*b1_dif_data) | (*b2_dif_data)) );
                b1_sim_data++;
                b1_dif_data++;
                b2_sim_data++;
                b2_dif_data++;
			//count += __builtin_popcountl((*b1_data++) ^ (*b2_data++));
  	  	}
   
    	return uint64_t(float(or_count - xor_count) / float(or_count) * size() );
	}
	else if (type == 0) {
		uint64_t count = 0;
		sdsl::bit_vector::size_type len = size()>>6;
                for (sdsl::bit_vector::size_type p = 0; p < len; ++p) {
                       count += __builtin_popcountl( ((*b1_sim_data) ^ (*b2_sim_data)) | ((*b1_dif_data) ^ (*b2_dif_data)) ); 
                }
		return size() - count;
	}
	
	DIE("ERROR: ONLY TWO TYPES IMPLEMENTED");
	return 0;
}


std::tuple<uint64_t, uint64_t> SBF::b_similarity(const BF* other) const {
    DIE("b_similarity has two bit vectors now. This method is not used");
    return std::make_tuple(0,0);
    /*
    assert(other->size() == size());

    const uint64_t* b1_sim_data = this->sim->data();
    const uint64_t* b1_dif_data = this->dif->data();

    const SBF* o = dynamic_cast<const SBF*>(other);
    if (o == nullptr) {
        DIE("Can only compute similarity on same type of BF.");
    }

    const uint64_t* b2_sim_data = o->sim->data();
    const uint64_t* b2_dif_data = o->dif->data();

    uint64_t and_count = 0;
    uint64_t or_count = 0;
    sdsl::bit_vector::size_type len = size()>>6;
    for (sdsl::bit_vector::size_type p = 0; p < len; ++p) {
        xor_count += __builtin_popcountl((*b1_sim_data) ^ (*b2_sim_data));
        or_count += __builtin_popcountl((*b1_data++) | (*b2_data++));
        //count += __builtin_popcountl((*b1_data++) ^ (*b2_data++));
    }
    return std::make_tuple(xor_count, or_count);
//size() - count;
    */
}

// return the # of 1s in the bitvector
uint64_t SBF::count_ones() const {
    uint64_t* sim_data = sim->data();
    uint64_t* dif_data = dif->data();
    sdsl::bit_vector::size_type len = size()>>6;
    uint64_t count_sim = 0;
    uint64_t count_dif = 0;
    for (sdsl::bit_vector::size_type p = 0; p < len; ++p) {
        int pc_sim = popcount(*sim_data);
        int pc_dif = popcount(*dif_data);
        DIE_IF(pc_sim != __builtin_popcountl(*sim_data), "popcountl and us disagree about popcount (sim)");
        DIE_IF(pc_dif != __builtin_popcountl(*dif_data), "popcountl and us disagree about popcount (dif)");
        sim_data++;
        dif_data++;
        count_sim += pc_sim;
        count_dif += pc_dif;
    }

    sdsl::rank_support_v<> rbv_sim(sim);
    sdsl::rank_support_v<> rbv_dif(dif);
    DIE_IF(rbv_sim.rank(size()) != count_sim, "SDSL and us disagree about number of 1s (sim)");
    DIE_IF(rbv_dif.rank(size()) != count_dif, "SDSL and us disagree about number of 1s (dif)");
    return count_sim+count_dif;
}

void SBF::compress() {
	sdsl::rrr_vector<255> rrr_sim(*sim);
    sdsl::rrr_vector<255> rrr_dif(*dif);
	std::cerr << "Compressed RRR sim vector is " << sdsl::size_in_mega_bytes(rrr_sim) << std::endl;
	std::cerr << "Compressed RRR diff vector is " << sdsl::size_in_mega_bytes(rrr_dif) << std::endl;
    sdsl::store_to_file(rrr_sim,filename+".rrr");
	sdsl::store_to_file(rrr_dif,f2+".rrr");
}


BF* SBF::sim_with(const std::string & new_name, const BF* f2) const{
    DIE("Not used - union is both sim and dif union");
    UncompressedBF* out = new UncompressedBF(new_name, hashes, num_hash);
    return out;
}


void SBF::sim_into(const BF* f2){
    DIE("Not used - union is both sim and dif union");
}

BF* SBF::dif_with(const std::string & new_name, const BF* f2) const{
    DIE("Not used - union is both sim and dif");
    UncompressedBF* out = new UncompressedBF(new_name, hashes, num_hash);
    return out;
}

void SBF::dif_into(const BF* f2){
    DIE("Not used - union is both sim and dif");
}




/*
Misc Code elements
*/
//====================================================================//
// union using 64bit integers
sdsl::bit_vector* union_bv_fast(const sdsl::bit_vector & b1, const sdsl::bit_vector& b2) {
    assert(b1.size() == b2.size());

    sdsl::bit_vector* out = new sdsl::bit_vector(b1.size(), 0);
    uint64_t* out_data = out->data();

    const uint64_t* b1_data = b1.data();
    const uint64_t* b2_data = b2.data();
    sdsl::bit_vector::size_type len = b1.size()>>6;
    for (sdsl::bit_vector::size_type p = 0; p < len; ++p) {
        (*out_data++) = (*b1_data++) | (*b2_data++);
    }
    return out;
}

BF* load_bf_from_file(const std::string & fn, HashPair hp, int nh) {
    if (fn.substr(fn.size()-4) == ".rrr") {
        return new BF(fn, hp, nh);
    } else if (fn.substr(fn.size()-3) == ".bv") {
        return new UncompressedBF(fn, hp, nh);
    } else {
        DIE("unknown bloom filter filetype (make sure extension is .rrr or .bv)");
        return nullptr;
    }
}

sdsl::bit_vector* sim_bv_fast(const sdsl::bit_vector & b1, const sdsl::bit_vector& b2){
    assert(b1.size() == b2.size());

    sdsl::bit_vector* out = new sdsl::bit_vector(b1.size(), 0);
    uint64_t* out_data = out->data();

    const uint64_t* b1_data = b1.data();
    const uint64_t* b2_data = b2.data();
    sdsl::bit_vector::size_type len = b1.size()>>6;
    for (sdsl::bit_vector::size_type p = 0; p < len; ++p) {
        (*out_data++) = (*b1_data++) & (*b2_data++);
    }
    return out;
}

sdsl::bit_vector* dif_bv_fast(const sdsl::bit_vector & b1, const sdsl::bit_vector& b2) {
    assert(b1.size() == b2.size());

    sdsl::bit_vector* out = new sdsl::bit_vector(b1.size(), 0);
    uint64_t* out_data = out->data();

    const uint64_t* b1_data = b1.data();
    const uint64_t* b2_data = b2.data();
    sdsl::bit_vector::size_type len = b1.size()>>6;
    for (sdsl::bit_vector::size_type p = 0; p < len; ++p) {
        (*out_data++) = (*b1_data++) ^ (*b2_data++);
    }
    return out;
}

std::string split_filename(const std::string & new_name){
    DIE("Not implemented");
    return "";
    /*
    std::vector<std::string> split_name;
    int count = SplitString(new_name, '.', split_name);
    for(std::vector<std::string>::const_iterator I = split_name.begin();
        I != split_name.end();
        ++I)
    {
        std::cerr << I << std::endl; 
    }
    */
}
