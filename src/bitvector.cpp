//
//  bitvector.cpp
//
//  Bitvector heavly inspired by Nicola Prezza's implementation.
//  https://github.com/nicolaprezza/r-index


#include <rimerge/bitvector.hpp>

namespace rimerge
{

//------------------------------------------------------------------------------

bitvector::bitvector(const sdsl::bit_vector& in)
{
    vector = vector_type(in);
    size_ = in.size();
    
    init_rank_and_select_structures();
}

bitvector::bitvector(sdsl::sd_vector_builder& builder) : vector(builder)
{
    size_ = vector.size();
    
    init_rank_and_select_structures();
}

bitvector::bitvector(const std::vector<bool>& in)
{
    sdsl::bit_vector bit_vector(in.size());
    
    for (std::size_t i = 0; i < in.size(); i++)
    {
        bit_vector[i] = in[i];
    }
    
    vector = vector_type(bit_vector);
    size_ = bit_vector.size();
    
    init_rank_and_select_structures();
}

bitvector&
bitvector::operator=(const bitvector& other)
{
    this->size_ = other.size();
    this->vector = vector_type(other.vector);
    
    init_rank_and_select_structures();
    
    return *this;
}

bool
bitvector::operator[](const size_type pos) const
{
    assert(pos < size_);
    return vector[pos];
}

bool
bitvector::at(size_type pos) const
{
    return operator[](pos);
}

size_type
bitvector::rank(size_type i) const
{
    return rank1(i);
}

size_type
bitvector::select(size_type i) const
{
    return select1(i + 1);
}

size_type
bitvector::predecessor(size_type i) const
{
    return select(rank(i) - 1);
}

size_type
bitvector::predecessor_rank(size_type i) const
{
    return rank(i) - 1;
}

size_type
bitvector::predecessor_rank_circular(size_type i) const
{
    return rank(i) == 0 ? number_of_1() - 1 : rank(i) - 1;
}

size_type
bitvector::gap_at(size_type i) const
{
    if (i == 0)
        return select(0) + 1;
    
    return select(i) - select(i - 1);
}

size_type
bitvector::size() const
{
    return size_;
}

size_type
bitvector::number_of_1() const
{
    if (size_ == 0)
        return 0;
    return rank1(size_);
}

size_type
bitvector::serialize(std::ostream& out) const
{
    size_type w_bytes = 0;
    
    out.write((char*) &size_, sizeof(size_));
    
    w_bytes += sizeof(size_);
    
    if (size_ == 0) return w_bytes;
    
    w_bytes += vector.serialize(out);
    
    return w_bytes;
}

void
bitvector::load(std::istream& in)
{
    in.read((char*) &size_, sizeof(size_));
    
    if (size_ == 0) return;
    
    vector.load(in);
    
    init_rank_and_select_structures();
}

//------------------------------------------------------------------------------

}
