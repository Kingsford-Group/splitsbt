#include "BF.h"
#include "minhash.h"
#include "util.h"
#include "bit_pointer.hpp"
//#include "xxHash-dev/xxhash.c"

using HashPair = jellyfish::hash_pair<jellyfish::mer_dna>;

std::vector<uint64_t> minhash(sdsl::bit_vector & bv, uint64_t seed, int numhash){
    const uint64_t* bv_data = bv.data();
    sdsl::bit_vector::size_type len = bv.size()>>6;

    srand(seed);
        

    std::vector<uint64_t> minhash(numhash);
    uint64_t minindex[numhash];
    uint64_t i = 0;
    int first=numhash;
    uint64_t xorarray[numhash];
    for (int hash=0; hash < numhash; hash++){
        xorarray[hash]=rand();
    }

    for (uint64_t q =0; q < bv.size(); q++){
        bool temp = bv[q];
        if (temp){
            for (int hash=0; hash < numhash; hash++){
                uint64_t value = bitshift_hash(q, xorarray[hash], hash);
                if (first>0){
                    minhash[hash] = value;
                    minindex[hash] = q;
                    first--;
                } else if (value < minhash[hash]){
                    minhash[hash] = value;
                    minindex[hash] = q;
                }
            }
        }
    }
    for (int hash=0; hash < numhash; hash++){
        std::cerr << "Temphead: " << hash << " " << minhash[hash] << " " << minindex[hash] << " " <<  xorarray[hash] << std::endl;
    }
    return minhash;
}

std::vector<uint64_t> minhash_fast(sdsl::bit_vector & bv, uint64_t seed, int numhash){
    const uint64_t* bv_data = bv.data();
    sdsl::bit_vector::size_type len = bv.size()>>6;
    srand(seed);
    // *** XXX: TEMP DEBUG VALUE ***    
    bool DEBUG=0;
    std::vector<uint64_t> minhash(numhash);
    uint64_t minindex[numhash];
    uint64_t i = 0;
    int first=numhash;
    uint64_t xorarray[numhash];
    for (int hash=0; hash < numhash; hash++){
        xorarray[hash]=rand();
    }

    /*
    uint64_t temp = 4398046511104;
    uint64_t mask=1;
    for (uint64_t q=0; q<64; ++q){
        uint64_t bit = temp & mask;
            if (bit!=0){
                uint64_t index = 63-q;
                std::cerr << q << " " << index << std::endl;
            }
        mask = mask << 1;
    }
    return 0;
    */
    for (sdsl::bit_vector::size_type p =0; p < len; ++p){
        uint64_t mask=1;
        for (uint64_t q=0; q<64; ++q){
            uint64_t bit = (*bv_data) & mask;
            if (bit!=0){
                uint64_t index = p*64+q;
                if (DEBUG && bv[index]!=1){
                    std::cerr << "index: " << index << std::endl;
                    std::cerr << "p,q: " << p << " " << q << std::endl;
                    std::cerr << "bit: " << bit << std::endl;
                    std::cerr << "mask: " << mask << std::endl;
                    std::cerr << "uint: " << (*bv_data) << std::endl;
                    std::cerr << "bvval: " << bv[index] << std::endl;
                    break;
                }
                for (int hash=0; hash < numhash; hash++){
                    uint64_t value = bitshift_hash(index, xorarray[hash], hash);
                    if (first>0){
                        minhash[hash] = value;
                        minindex[hash] = index;
                        first--;
                    } else if (value < minhash[hash]){
                        minhash[hash] = value;
                        minindex[hash] = index;
                    }
                }
            } else {
                uint64_t index = p*64+q;
                if (DEBUG && bv[index]!=0){
                    std::cerr << "p,q: " << p << " " << q << std::endl;
                    std::cerr << "bit: " << bit << std::endl;
                    std::cerr << "mask: " << mask << std::endl;
                    std::cerr << "uint: " << (*bv_data) << std::endl;
                    std::cerr << "bvval: " << bv[index] << std::endl;

                }
            }
            mask = mask << 1;
        }
        bv_data++;
        //if (count64==10){
        //    break;
        //}
            //std::cerr << "Minhash value is: " << minhash << " at position " << minindex << " (Seed: " << seed << ")" << std::endl;
        //std::cerr << "Temphead: " << hash << " " << minhash << " " << minindex << " " <<  seed << std::endl;
    }
    //for (int hash=0; hash < numhash; hash++){
    //    std::cerr << "Temphead: " << hash << " " << minhash[hash] << " " << minindex[hash] << " " <<  xorarray[hash] << std::endl;
    //}
    return minhash;
}

int split_minhash(sdsl::bit_vector & bv, uint64_t seed, int numhash){
    const uint64_t* bv_data = bv.data();
    sdsl::bit_vector::size_type len = bv.size()>>6;

    srand(seed);
        

    std::size_t minhash[numhash];
    uint64_t minindex[numhash];
    uint64_t i = 0;
    int first=numhash;
    uint64_t xorarray[numhash];
    for (int hash=0; hash < numhash; hash++){
        xorarray[hash]=rand();
    }

    for (uint64_t q =0; q < bv.size(); q++){
        bool temp = bv[q];
        if (temp){
            for (int hash=0; hash < numhash; hash++){
                uint64_t value = bitshift_hash(q, xorarray[hash], hash);
                if (first>0){
                    minhash[hash] = value;
                    minindex[hash] = q;
                    first--;
                } else if (value < minhash[hash]){
                    minhash[hash] = value;
                    minindex[hash] = q;
                }
            }
        }
    }
    for (int hash=0; hash < numhash; hash++){
        std::cerr << "Temphead: " << hash << " " << minhash[hash] << " " << minindex[hash] << " " <<  xorarray[hash] << std::endl;
    }
    //return minhash;
}


void sbthash(sdsl::bit_vector & bv, uint64_t* top100){
    const uint64_t* bv_data = bv.data();
    sdsl::bit_vector::size_type len = bv.size()>>6;

    //top100=uint64_t[100];
    int count = 0;
    for (uint64_t q =0; q < bv.size(); q++){
        bool temp = bv[q];
        if (temp){
            top100[count]=q;
            count++;
        }
        if (count == 100){
            break;
        }
    }
    //uint64_t *count100;
    
}


// Bitshift table taken from Marsaglia_2003_xorshift
uint64_t bitshift_hash(uint64_t bit_index, uint64_t seed, int hash){
        int xorshift_array[275][3]={{1, 1,54},{1, 1,55},{1, 3,45},{1, 7, 9},{1, 7,44},{1, 7,46},{1, 9,50},{1,11,35},{1,11,50},
                                    {1,13,45},{1,15, 4},{1,15,63},{1,19, 6},{1,19,16},{1,23,14},{1,23,29},{1,29,34},{1,35, 5},
                                    {1,35,11},{1,35,34},{1,45,37},{1,51,13},{1,53, 3},{1,59,14},{2,13,23},{2,31,51},{2,31,53},
                                    {2,43,27},{2,47,49},{3, 1,11},{3, 5,21},{3,13,59},{3,21,31},{3,25,20},{3,25,31},{3,25,56},
                                    {3,29,40},{3,29,47},{3,29,49},{3,35,14},{3,37,17},{3,43, 4},{3,43, 6},{3,43,11},{3,51,16},
                                    {3,53, 7},{3,61,17},{3,61,26},{4, 7,19},{4, 9,13},{4,15,51},{4,15,53},{4,29,45},{4,29,49},
                                    {4,31,33},{4,35,15},{4,35,21},{4,37,11},{4,37,21},{4,41,19},{4,41,45},{4,43,21},{4,43,31},
                                    {4,53, 7},{5, 9,23},{5,11,54},{5,15,27},{5,17,11},{5,23,36},{5,33,29},{5,41,20},{5,45,16},
                                    {5,47,23},{5,53,20},{5,59,33},{5,59,35},{5,59,63},{6, 1,17},{6, 3,49},{6,17,47},{6,23,27},
                                    {6,27, 7},{6,43,21},{6,49,29},{6,55,17},{7, 5,41},{7, 5,47},{7, 5,55},{7, 7,20},{7, 9,38},
                                    {7,11,10},{7,11,35},{7,13,58},{7,19,17},{7,19,54},{7,23, 8},{7,25,58},{7,27,59},{7,33, 8},
                                    {7,41,40},{7,43,28},{7,51,24},{7,57,12},{8, 5,59},{8, 9,25},{8,13,25},{8,13,61},{8,15,21},
                                    {8,25,59},{8,29,19},{8,31,17},{8,37,21},{8,51,21},{9, 1,27},{9, 5,36},{9, 5,43},{9, 7,18},
                                    {9,19,18},{9,21,11},{9,21,20},{9,21,40},{9,23,57},{9,27,10},{9,29,12},{9,29,37},{9,37,31},
                                    {9,41,45},{10, 7,33},{10,27,59},{10,53,13},{11, 5,32},{11, 5,34},{11, 5,43},{11, 5,45},{11, 9,14},
                                    {11, 9,34},{11,13,40},{11,15,37},{11,23,42},{11,23,56},{11,25,48},{11,27,26},{11,29,14},{11,31,18},
                                    {11,53,23},{12, 1,31},{12, 3,13},{12, 3,49},{12, 7,13},{12,11,47},{12,25,27},{12,39,49},{12,43,19},
                                    {13, 3,40},{13, 3,53},{13, 7,17},{13, 9,15},{13, 9,50},{13,13,19},{13,17,43},{13,19,28},{13,19,47},
                                    {13,21,18},{13,21,49},{13,29,35},{13,35,30},{13,35,38},{13,47,23},{13,51,21},{14,13,17},{14,15,19},
                                    {14,23,33},{14,31,45},{14,47,15},{15, 1,19},{15, 5,37},{15,13,28},{15,13,52},{15,17,27},{15,19,63},
                                    {15,21,46},{15,23,23},{15,45,17},{15,47,16},{15,49,26},{16, 5,17},{16, 7,39},{16,11,19},{16,11,27},
                                    {16,13,55},{16,21,35},{16,25,43},{16,27,53},{16,47,17},{17,15,58},{17,23,29},{17,23,51},{17,23,52},
                                    {17,27,22},{17,45,22},{17,47,28},{17,47,29},{17,47,54},{18, 1,25},{18, 3,43},{18,19,19},{18,25,21},
                                    {18,41,23},{19, 7,36},{19, 7,55},{19,13,37},{19,15,46},{19,21,52},{19,25,20},{19,41,21},{19,43,27},
                                    {20, 1,31},{20, 5,29},{21, 1,27},{21, 9,29},{21,13,52},{21,15,28},{21,15,29},{21,17,24},{21,17,30},
                                    {21,17,48},{21,21,32},{21,21,34},{21,21,37},{21,21,38},{21,21,40},{21,21,41},{21,21,43},{21,41,23},
                                    {22, 3,39},{23, 9,38},{23, 9,48},{23, 9,57},{23,13,38},{23,13,58},{23,13,61},{23,17,25},{23,17,54},
                                    {23,17,56},{23,17,62},{23,41,34},{23,41,51},{24, 9,35},{24,11,29},{24,25,25},{24,31,35},{25, 7,46},
                                    {25, 7,49},{25, 9,39},{25,11,57},{25,13,29},{25,13,39},{25,13,62},{25,15,47},{25,21,44},{25,27,27},
                                    {25,27,53},{25,33,36},{25,39,54},{28, 9,55},{28,11,53},{29,27,37},{31, 1,51},{31,25,37},{31,27,35},
                                    {33,31,43},{33,31,55},{43,21,46},{49,15,61},{55, 9,56}};

    int hashvals[3]={xorshift_array[hash][0],xorshift_array[hash][1],xorshift_array[hash][2]};
    uint64_t outval = bit_index^seed;
    //std::cerr << bit_index << " " << outval << std::endl;
    outval^=(outval<<hashvals[0]); outval^=(outval>>hashvals[1]); outval^=(outval<<hashvals[2]);
    outval^=(outval<<hashvals[2]); outval^=(outval>>hashvals[1]); outval^=(outval<<hashvals[0]);
    outval^=(outval>>hashvals[0]); outval^=(outval<<hashvals[1]); outval^=(outval>>hashvals[2]);
    outval^=(outval>>hashvals[2]); outval^=(outval<<hashvals[1]); outval^=(outval>>hashvals[0]);
    outval^=(outval<<hashvals[0]); outval^=(outval<<hashvals[2]); outval^=(outval>>hashvals[1]);
    outval^=(outval<<hashvals[2]); outval^=(outval<<hashvals[0]); outval^=(outval>>hashvals[1]);
    outval^=(outval>>hashvals[0]); outval^=(outval>>hashvals[2]); outval^=(outval<<hashvals[1]);
    outval^=(outval>>hashvals[2]); outval^=(outval>>hashvals[0]); outval^=(outval<<hashvals[1]);
    return outval;
    //std::cerr<< outval<<std::endl;

}

void write_minhash(std::vector<uint64_t> myhash, std::string bvfile1){
        std::string parse = bvfile1.substr(0, bvfile1.find_first_of("."));
        std::cerr << parse+".minhash" << std::endl;
        std::ofstream out_file(parse+".minhash");
        std::ostream_iterator<uint64_t> output_iterator(out_file, "\n");
        std::copy(myhash.begin(), myhash.end(), output_iterator);
}

std::vector<uint64_t> read_minhash(std::string file, int numhash){
    std::vector<uint64_t> minhash(numhash);
    std::ifstream infile(file);
    std::string line;
    int i = 0;
    while (infile >> minhash[i]){
        i++;
    }
    //for (i=0; i<numhash; i++){    
    //    std::cerr << minhash[i] << " ";
    //}
    //std::cerr << std::endl;
    return minhash;
}
/*
int minhash_sim(std::vector<uint64_t> & mh1, std::vector<uint64_t> &mh2, int nh){
    int sim = 0;
    for(int i = 0; i < nh; i++){
        if( mh1[i] == mh2[i]){
            sim++;
        }
    }
    return sim;
}
*/

double minhash_sim(std::vector<uint64_t> & mh1, std::vector<uint64_t> &mh2, int nh){
    double sim = 0.0;
    for(int i = 0; i < nh; i++){
        if( mh1[i] == mh2[i]){
            sim++;
        }
    }
    return sim / nh;
}

mh_node::mh_node(std::string file, int numhash){
    std::string parse = file.substr(0, file.find_first_of("."));
    fname = parse+".sim.bf.bv";
    minhash = read_minhash(file, numhash);
}

MHcluster::MHcluster(std::string minlist, int numhash, std::string out_instruct){
    // read in minhashes, build array of them
    nh = numhash;
    outfile = out_instruct;

    std::ifstream list_file(minlist);
    std::string file;
    std::vector<uint64_t> minhash;

    num_elements = 0;
    while (list_file >> file){
        std::cerr << num_elements << std::endl;
        std::ifstream infile(file);
        mh_node* mynode = new mh_node(file, numhash);
        mhvector.push_back(mynode);
        num_elements++;
    }

    //build the NxN distance matrix
    std::cerr << "Calculating distance matrix..." << std::endl;
    calcDMatrix();
}

MHcluster::~MHcluster(){
    for(int i = 0; i < num_elements; i++){
        delete mhvector[i];
    }

    if (distance_matrix != NULL){
        for (int i = 1; i < num_elements; i++){
            if (distance_matrix[i] != NULL) free(distance_matrix[i]);
        }
        free(distance_matrix);
        distance_matrix = NULL;
    }
}

bool MHcluster::calcDMatrix(){
    int i,j;
    int n = num_elements;
    if (n < 2) return false;
    
    distance_matrix = (double**)malloc(n*sizeof(double*));
    if(distance_matrix==NULL) return NULL;
    distance_matrix[0]=NULL;
    for(i = 1; i < n; i++){
        distance_matrix[i]= (double*)malloc(i*sizeof(double));
        if (distance_matrix[i]==NULL) break;
    }

    if (i < n){
        j=1;
        for(int i = 1; i < j; i++) free(distance_matrix[i]);
        free(distance_matrix);
        distance_matrix=NULL;
        std::cerr << "Error. Failed to allocate memory." << std::endl;
        return false;
    }
    std::cerr << "Seeding distance table" << std::endl;
    for(i=1; i < num_elements; i++){
        for(j=0; j < i; j++){
            distance_matrix[i][j]=1.0-minhash_sim(mhvector[i]->minhash,mhvector[j]->minhash, nh);
            std::cerr << i << " " << j << " " << distance_matrix[i][j] << std::endl;
        }
    }
    return true;
}

//added idnex / matrix file to construct SSBT index while clustering
void MHcluster::gcluster(const std::string &indexfile, const std::string & matrix_file){
    // perform allsome greedy cluster
    //allsome_gclust(num_elements, dmatrix, outfile, nh);

// This is a re-implementation of the greedy clustering strategy implemented by ALLSOME SBT
// Instead of building the tree directly, it outputs the instructions to a file for later use
//void allsome_gclust(int nelements, double** distmatrix, std::string outfile, int nh){
    std::map<std::string, BloomTree*> topo;

    std::cerr << "Gclustering " << num_elements << " elements" << std::endl;
    //Initialize topology for leaves
    for(int i=0; i < num_elements; i++){
        topo[mhvector[i]->fname]=new BloomTree(mhvector[i]->fname);
    }

    int j;
    int n;
    std::ofstream out(outfile.c_str());

    BloomTree* newnode;
    for (n = num_elements; n> 1; n--){
        int is = 1;
        int js = 0;
        double dist = find_closest_pair(n, distance_matrix, &is, &js);
        std::cerr << "Best: " << is << ", " << js << " " << dist << std::endl;
        // record instruction (new node, child 1, child 2)
        std::string parse = mhvector[is]->fname.substr(0, mhvector[is]->fname.find_last_of("/")+1);
        std::string newname = parse+"internal_"+std::to_string(n)+"_union.sim.bf.bv";
        out << newname <<  " " << mhvector[is]->fname << " " << mhvector[js]->fname << std::endl; 


        //Build index here
        newnode = new BloomTree(newname);
        newnode->set_child(0, topo[mhvector[is]->fname]);
        newnode->set_child(1, topo[mhvector[js]->fname]);
        topo[newname]=newnode;

        //std::cerr << topo[newname]->name() << " " << topo[newname]->child(0)->name() << " " << topo[newname]->child(1)->name() << std::endl;


        if(js == n-1){
            int tmp = is;
            is = js;
            js = tmp;
        }

        if(is != n-1){
            unionAB2A_copyC2B(mhvector[js], mhvector[is], mhvector[n-1], newname);
        } else{
            unionAB2A(mhvector[js], mhvector[is], newname);
        }

        for (j = 0; j < js; j++){
            distance_matrix[js][j] = minhash_sim(mhvector[js]->minhash,mhvector[j]->minhash, nh);
        }
        for (j = js+1; j < is; j++){
            distance_matrix[j][js] = minhash_sim(mhvector[js]->minhash,mhvector[j]->minhash, nh);
        }
        for (j = is+1; j < n; j++){
            distance_matrix[j][js] = minhash_sim(mhvector[js]->minhash,mhvector[j]->minhash, nh);
        }

        if(is != n-1){
            for(j=0; j < is; j++) distance_matrix[is][j] = distance_matrix[n-1][j];
            for(j=is+1; j < n-1; j++) distance_matrix[j][is] = distance_matrix[n-1][j];
        }
    }
    write_bloom_tree(indexfile, newnode, matrix_file);

    // delete topology?
}

void MHcluster::unionAB2A_copyC2B(mh_node* a, mh_node* b, mh_node* c, std::string newn){
    a->fname=newn;
    b->fname=c->fname;

    for (int i = 0 ; i < nh; i++){
        uint64_t aval = a->minhash[i];
        uint64_t bval = b->minhash[i];
        if(aval < bval){
            a->minhash[i]=bval; //else a stays the same
        }
        b->minhash[i]=c->minhash[i]; 
    } 
}

void MHcluster::unionAB2A(mh_node* a, mh_node* b, std::string newn){
    a->fname=newn;
    
    for (int i = 0 ; i < nh; i++){
        uint64_t aval = a->minhash[i];
        uint64_t bval = b->minhash[i];
        if(aval < bval){
            a->minhash[i]=bval; //else a stays the same
        } 
    }
}

//Taken from cluster v 3.0
double find_closest_pair(int n, double** distmatrix, int* ip, int* jp)
{ int i, j;
  double temp;
  double distance = distmatrix[1][0];
  *ip = 1;
  *jp = 0;
  for (i = 1; i < n; i++)
  { for (j = 0; j < i; j++)
    { temp = distmatrix[i][j];
      if (temp<distance)
      { distance = temp;
        *ip = i;
        *jp = j;
      }
    }
  }
  return distance;
}

/*
    int count64=0;
    for (sdsl::bit_vector::size_type p = 0; p < len; ++p){
        count64+=1;
        uint64_t mask=1;
        for (uint64_t q=0; q<64; q++){
            uint64_t bit = (*bv_data) & mask;
            if (popcount(bit)==1) {
                std::cerr << 1; 
            } else {
                std::cerr << 0;
            }
            mask = mask << 1;
        }
        bv_data++;
        std::cerr << std::endl;
        if (count64==2){
            break;
        }
    } 
*/

