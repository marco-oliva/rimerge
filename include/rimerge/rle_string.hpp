//
//  rindex.hpp
//
//  Copyright 2020 Nicola Prezza. All rights reserved.
//

#ifndef rle_string_hpp
#define rle_string_hpp

#include <rimerge/utils.hpp>

namespace rimerge
{

//------------------------------------------------------------------------------

// Huffman string
class huff_string
{
public:
    
    huff_string() = default;
    
    huff_string(string_type &s)
    {
        s.push_back(IMPL_TERMINATOR);
        construct_im(wt, (const char*) s.data(), 1);
        assert(wt.size()==s.size() - 1);
    }
    
    byte_type operator[](size_type i) const
    {
        assert(i<wt.size());
        return wt[i];
    }
    
    size_type size() const
    {
        return wt.size();
    }
    
    size_type rank(size_type i, byte_type c) const
    {
        assert(i<=wt.size());
        return wt.rank(i,c);
    }
    
    // position of i-th character c. i starts from 0!
    size_type select(size_type i, byte_type c) const
    {
        return wt.select(i+1,c);
    }
    
    // serialize the structure to the ostream
    size_type serialize(std::ostream& out) const
    {
        return wt.serialize(out);
    }
    
    // load the structure from the istream
    void load(std::istream& in)
    {
        wt.load(in);
    }

private:
    
    sdsl::wt_huff<> wt;
};


//------------------------------------------------------------------------------

// Bitvector
class bitvector
{
public:
    
    bitvector() = default;
    
    /// Constructor, from a sdsl::bit_vector
    bitvector(const sdsl::bit_vector& in);
    
    /// Constructor, from a std::vector<bool>
    bitvector(const std::vector<bool>& in);
    
    bitvector& operator= (const bitvector& other);
    
    /// Subscript oprator, pos must be < size
    bool operator[](const size_type pos) const;
    
    /// At, alias for subscript operator
    bool at(const size_type pos) const;
    
    /// Rank_1, i must be < size
    size_type rank(const size_type i) const;
    
    /// Select_1, position og the i-th one in the vector
    size_type select(const size_type i) const;
    
    /// Predecessor, i must have a predecessor
    size_type predecessor(const size_type i) const;
    
    /// Predecessor of i (position i excluded) in rank pace. i must have a predecessor
    size_type predecessor_rank(const size_type i) const;
    
    size_type predecessor_rank_circular(const size_type i) const;
    
    /// Length of the i-th gap. Includes the leading 1
    size_type gap_at(const size_type i) const;
    
    /// Size
    size_type size() const;
    
    /// Number of 1s
    size_type number_of_1() const;
    
    /// Serialize to output
    size_type serialize(std::ostream& out) const;
    
    /// Load from input
    void load(std::istream& in);


private:
    
    void init_rank_and_select_structures()
    {
        rank1   = vector_type::rank_1_type(&vector);
        select1 = vector_type::select_1_type(&vector);
    }
    
    size_type size_ = 0;
    
    vector_type                  vector;
    vector_type::rank_1_type      rank1;
    vector_type::select_1_type  select1;
    
};

//------------------------------------------------------------------------------

class rle_string
{

public:
    
    void construct(byte_type* start, size_type size, size_type block_size);
    
    rle_string() = default;
    
    rle_string(mio::mmap_source& input, size_type B = 2);
    
    rle_string(string_type& input, size_type B = 2);
    
    byte_type operator[](size_type i) const;
    
    // position of i-th character c. i starts from 0!
    size_type select(size_type i, byte_type c) const;
    
    // number of c before position i
    size_type rank(size_type i, byte_type c) const;
    
    // text position i is inside this run
    size_type run_of_position(size_type i) const;
    
    //break range: given a range <l',r'> on the string and a character c, this function
    //breaks <l',r'> in maximal sub-ranges containing character c.
    //for simplicity and efficiency, we assume that characters at range extremities are both 'c'
    //thanks to the encoding (run-length), this function is quite efficient: O(|result|) ranks and selects
    std::vector<range_type> break_range(range_type rn, byte_type c) const;
    
    size_type size() const;
    
    // inclusive range of j-th run in the string
    std::pair<size_type,size_type> run_range(size_type j) const;
    
    // length of i-th run
    size_type run_at(size_type i) const;
    
    size_type number_of_runs() const;
    
    // serialize the structure to the ostream
    size_type serialize(std::ostream& out) const;
    
    // load the structure from the istream
    void load(std::istream& in);
    
    string_type to_string() const;
    
    /*
     * input: inclusive range rn, character c
     *
     * return the position j that is closest to rn.starts,
     * such that character in position j is c and that is
     * adjacent to a position j' inside rn that contains a
     * character != c
     *
     * rn must contain c and at least another character d!=c
     *
     */
    size_type closest_run_break(range_type rn, byte_type c) const;
    
    friend class rindex;

private:
    
    // <j=run of position i, ends position of j-th run>
    std::pair<size_type,size_type> run_of(size_type i) const;
    
    // block size: bitvector 'runs' has R/B bits set (R being number of runs)
    std::size_t B = 2;
    
    bitvector runs;
    
    // for each letter, its runs stored contiguously
    std::vector<bitvector> runs_per_letter;
    
    // store run heads in a compressed string supporting access/rank
    huff_string run_heads;
    
    // text length and number of runs
    std::size_t n = 0;
    std::size_t R = 0;
    
};
}

#endif //rle_string_hpp
