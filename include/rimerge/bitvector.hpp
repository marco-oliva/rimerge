//
//  bitvector.hpp
//
//  Copyright 2020 Marco Oliva. All rights reserved.
//

#ifndef bitvector_hpp
#define bitvector_hpp

#include <rimerge/utils.hpp>

namespace rimerge
{

class bitvector
{
public:
    
    bitvector() = default;
    
    /// Constructor, from a sdsl::bit_vector
    bitvector(const sdsl::bit_vector& in);
    
    /// Constructor, from a sdsl::sd_vector_builder
    bitvector(sdsl::sd_vector_builder& builder);
    
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

}



#endif //bitvector_hpp
