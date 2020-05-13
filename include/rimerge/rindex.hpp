//
//  rindex.hpp
//
//  Copyright 2020 Nicola Prezza. All rights reserved.
//

#ifndef r_index_hpp
#define r_index_hpp

#include <rimerge/utils.hpp>
#include <rimerge/rle_string.hpp>
#include <rimerge/bwtmerge.hpp>

namespace rimerge
{

//------------------------------------------------------------------------------

class rindex;

//------------------------------------------------------------------------------

struct Alphabet
{
public:
    
    std::vector<bool> used_chars = std::vector<bool>(256,false);
    bool initialized = false;
    
    bitvector alphabet;
    
    Alphabet() = default;
    
    size_type sigma() const;
    void update(byte_type i);
    void init();
    byte_type previous(byte_type c) const;
    byte_type following(byte_type c) const;
    
};

//------------------------------------------------------------------------------

struct SA_samples
{
public:
    
    using sample_type = std::pair<size_type, size_type>;
    
    std::vector<sample_type> starts;
    std::vector<sample_type>   ends;
    
    SA_samples() = default;
    SA_samples(const std::vector<sample_type>& s, const std::vector<sample_type>& e);
    
    size_type operator[](const size_type i) const;
    
    constexpr static size_type SAMPLE_BYTES = 5;
    static void write(sample_type& s, std::ofstream& out);
};

//------------------------------------------------------------------------------

class SA_updates
{
public:
    
    //                             RA value             SA[RA[j]]  SA[RA[j+1]]
    using left_map_type = std::map<size_type, std::pair<size_type, size_type>>;
    //                              RA value                       min j     j-th sa               max j     j-th sa
    using right_map_type = std::map<size_type, std::pair<std::pair<size_type,size_type>, std::pair<size_type, size_type>>>;
    
    std::vector<left_map_type> left_samples;
    std::vector<right_map_type> right_samples;
    
    SA_updates(size_type threads) :left_samples(threads), right_samples(threads) {}
    
    void update_left(const rindex& left, const rindex& rigth, size_type ra_i, size_type ra_j, size_type i, size_type j, size_type thread);
    void update_right(const rindex& left, const rindex& right, size_type ra_j, size_type j, size_type sa_value, size_type thread);
    
};

//------------------------------------------------------------------------------

class rindex{

public:
    
    byte_type terminator = DATA_TERMINATOR;
    
    // Constructors and Destructors
    rindex() = default;
    rindex(const rle_string& b, const SA_samples& sas);
    rindex(const rindex& source);
    rindex(const rindex& a, const rindex& b);
    
    // From disk
    rindex(byte_type* bwt_start, size_type bwt_size,
           byte_type* ssa_start, size_type ssa_size,
           byte_type* esa_start, size_type esa_size);
    
    // From disk, without samples
    rindex(byte_type* bwt_start, size_type bwt_size);
    
    ~rindex() = default;
    
    void init(byte_type* bwt_start, size_type bwt_size,
              byte_type* ssa_start, size_type ssa_size,
              byte_type* esa_start, size_type esa_size);
    
    // Statistics
    size_type size() const { return this->bwt_.size(); }
    bool empty() const { return (this->size() == 0); }
    size_type sequences() const { return this->sequences_; }
    size_type sigma() const { return this->alphabet.sigma(); }
    size_type runs() const { return this->bwt_.number_of_runs(); }
    size_type end_marker() const { return this->data_terminator_position; }
    
    // Metadata
    Alphabet alphabet;
    
    // Data
    const rle_string& bwt() const { return bwt_; }
    const SA_samples& samples() const {return samples_; }
    byte_type operator[](size_type i) const { return bwt_[i]; }
    size_type rank(size_type i, byte_type c) const { return bwt_.rank(i,c); }
    size_type select(size_type i, byte_type c) const {return bwt_.select(i,c); }
    size_type get_terminator_position() const { return data_terminator_position; }
    
    // get full BWT range
    range_type full_range();
    
    // access column F at position i
    byte_type F_at(size_type  i) const;
    
    // backward navigation of the BWT
    size_type LF(size_type i) const;
    
    // backward navigation of the BWT
    size_type LF(size_type i, byte_type c) const;
    
    // backward navigation of the BWt, range
    range_type LF(range_type rn, byte_type c) const;
    
    // forward navigation of the BWT
    size_type FL(size_type i) const;
    
    // forward navigation of the BWT, where for efficiency we give c=F[i] as input
    size_type FL(size_type i, byte_type c) const;
    
    // Return BWT range of pattern P
    range_type count(std::string &P) const;
    
    // serialize the structure to the ostream
    size_type serialize(std::ostream& out) const;
    
    // load the structure from the istream
    void load(const std::istream& in);
    
    size_type text_size() const;
    
    // Merge
    static void merge(const rindex& left, const rindex& right, MergeParameters& parameters);
    
    
private:
    
    void read_samples(byte_type* sa_start, size_type sa_size, std::vector<SA_samples::sample_type>& samples);
    void read_bwt(byte_type* bwt_start, size_type bwt_size);
    void samples_from_bwt(rle_string& bwt);
    
    /*
     * sparse RLBWT: r (log sigma + (1+epsilon) * log (n/r)) (1+o(1)) bits
     */
    
    // F column of the BWT (vector of 256 elements)
    std::array<size_type, 256> F = {0};
    
    // L column of the BWT, run-length compressed
    rle_string bwt_;
    size_type  data_terminator_position = 0;
    
    // SA samples
    SA_samples samples_;
    
    size_type sequences_ = 0;
    size_type alphabet_size_ = 0;
};


//------------------------------------------------------------------------------

namespace testing
{

class rindex_generator
{
public:
    
    rindex_generator()
    {
        original = {'G','A','T','T','A','C','A','T',
                    STRING_TERMINATOR,
                    'G','A','T','A','C','A','T',
                    STRING_TERMINATOR,
                    'G','A','T','T','A','G','A','T','A',
                    DATA_TERMINATOR};
        
        bwt = {'A','T','T','T','T','T','T','C','C','G','G','G','G','A','A','A',
               STRING_TERMINATOR,
               DATA_TERMINATOR,
               STRING_TERMINATOR,
               'A','A','A','T','A','T','A','A'};
        
        rle_bwt = rle_string(bwt);
        r = 14;
        
        SA = {26,8,16,25,4,12,21,6,14,23,10,1,18,5,13,22,9,0,17,7,15,24,3,11,20,2,19};
    
        std::vector<std::pair<size_type,size_type>> first;
        
        std::vector<std::pair<size_type,size_type>>  last;
        last.push_back(std::make_pair(0, 26));
        
        for (size_type j = 0; j < SA.size(); j++)
            if (bwt[j] != bwt[j + 1])
                first.push_back(std::make_pair(j, SA[j]));
    
        for (size_type j = 1; j < SA.size(); j++)
            if (bwt[j - 1] != bwt[j])
                last.push_back(std::make_pair(j, SA[j]));
            
    }
    
    rindex get()
    {
        return rindex(rle_bwt, samples);
    }
    
    string_type original;
    string_type      bwt;
    rle_string   rle_bwt;
    size_type          r;
    
    std::vector<size_type>       SA;
    SA_samples              samples;
    
};

}

}

#endif //r_index_hpp
