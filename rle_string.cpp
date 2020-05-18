//
//  rle_string.cpp
//
//  Copyright 2020 Nicola Prezza. All rights reserved.
//

#include <rimerge/rle_string.hpp>


namespace rimerge
{

//------------------------------------------------------------------------------

bitvector::bitvector(const sdsl::bit_vector& in)
{
    vector = vector_type(in);
    size_  = in.size();
    
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
    size_  = bit_vector.size();
    
    init_rank_and_select_structures();
}

bitvector& bitvector::operator= (const bitvector& other)
{
    this->size_  = other.size();
    this->vector = vector_type(other.vector);
    
    init_rank_and_select_structures();
    
    return *this;
}

bool bitvector::operator[](const size_type pos) const
{
    assert(pos < size_);
    return vector[pos];
}

bool bitvector::at(size_type pos) const
{
    return operator[](pos);
}

size_type bitvector::rank(size_type i) const
{
    return rank1(i);
}

size_type bitvector::select(size_type i) const
{
    return select1(i + 1);
}

size_type bitvector::predecessor(size_type i) const
{
    return select(rank(i)-1);
}

size_type bitvector::predecessor_rank(size_type i) const
{
    return rank(i) - 1;
}

size_type bitvector::predecessor_rank_circular(size_type i) const
{
    return rank(i) == 0 ? number_of_1() - 1 : rank(i) - 1;
}

size_type bitvector::gap_at(size_type i) const
{
    if(i == 0)
        return select(0) + 1;
    
    return select(i) - select(i - 1);
}

size_type bitvector::size() const
{
    return size_;
}

size_type bitvector::number_of_1() const
{
    if (size_ == 0)
        return 0;
    return rank1(size_);
}

size_type bitvector::serialize(std::ostream& out) const
{
    size_type w_bytes = 0;
    
    out.write((char*)& size_, sizeof(size_));
    
    w_bytes += sizeof(size_);
    
    if(size_ == 0) return w_bytes;
    
    w_bytes += vector.serialize(out);
    
    return w_bytes;
}

void bitvector::load(std::istream& in)
{
    in.read((char*)& size_, sizeof(size_));
    
    if(size_ == 0) return;
    
    vector.load(in);
    
    init_rank_and_select_structures();
}

//------------------------------------------------------------------------------

void rle_string::construct(byte_type* start, size_type size, size_type block_size)
{
    this->B = block_size;
    this->n = 1;
    this->R = 0;
    
    if (size <= 1)
    {
        this->n = 0;
        return;
    }
    
    byte_type* iterator = start;
    byte_type* end = start + size;
    
    auto runs_per_letter_bv = std::vector<std::vector<bool>>(256);
    
    //runs in main bitvector
    std::vector<bool> runs_bv;
    runs_bv.reserve(size);
    string_type run_heads_s;
    
    byte_type last_c = *iterator;
    byte_type curr_c;
    
    while (++iterator != end)
    {
        curr_c = *iterator;
        if (last_c == 0) last_c = DATA_TERMINATOR;
        if (curr_c != last_c)
        {
            run_heads_s.push_back(last_c);
            runs_per_letter_bv[last_c].push_back(true);
            
            last_c = curr_c;
            
            //push back a bit set only at the end of a block
            runs_bv.push_back((R % B) == (B - 1));
            
            this->R++;
        }
        else
        {
            runs_bv.push_back(false);
            runs_per_letter_bv[last_c].push_back(false);
        }
        
        
        this->n++;
    }
    
    run_heads_s.push_back(last_c);
    runs_per_letter_bv[last_c].push_back(true);
    runs_bv.push_back(false);
    this->R++;
    
    runs = bitvector(runs_bv);
    
    //a fast direct array: char -> bitvector.
    runs_per_letter = std::vector<bitvector>(256);
    for(size_type i=0;i<256;++i)
        runs_per_letter[i] = bitvector(runs_per_letter_bv[i]);
    
    run_heads = huff_string(run_heads_s);
    
    assert(run_heads.size()==R);
}

rle_string::rle_string(mio::mmap_source& input, size_type B)
{
    construct((byte_type*) input.data(), input.size(), B);
}

rle_string::rle_string(string_type& input, size_type B)
{
    construct((byte_type*) input.data(), input.size(), B);
}

byte_type rle_string::operator[](size_type i) const
{
    return run_heads[run_of(i).first];
}

size_type rle_string::select(size_type i, byte_type c) const
{
    assert(i<runs_per_letter[c].size());
    
    //i-th c is inside j-th c-run (j starts from 0)
    assert(i<runs_per_letter[c].size());
    size_type j = runs_per_letter[c].rank(i);
    
    //starting position of i-th c inside its run
    assert(j==0 || i >= runs_per_letter[c].select(j-1) + 1);
    size_type before = (j==0 ? i : i - (runs_per_letter[c].select(j-1) + 1));
    
    //position in run_heads
    size_type r = run_heads.select(j,c);
    
    //k = number of bits before position of interest in the main string
    //here, k is initialized looking at the sampled runs
    assert(r/B==0 || r/B-1<runs.number_of_1());
    size_type k = (r/B==0?0 : runs.select(r/B-1)+1);
    
    //now add remaining run lengths to k
    for( size_type t = (r/B)*B; t<r; ++t ){
        
        k += run_at(t);
        
    }
    
    return k + before;
}

size_type rle_string::rank(size_type i, byte_type c) const
{
    assert(i<=n);
    
    //letter does not exist in the text
    if(runs_per_letter[c].size()==0) return 0;
    
    if(i==n) return runs_per_letter[c].size();
    
    size_type last_block = runs.rank(i);
    size_type current_run = last_block*B;
    
    //current position in the string: the starts of a block
    size_type pos = 0;
    if(last_block>0)
        pos = runs.select(last_block-1)+1;
    
    assert(pos <= i);
    
    size_type dist = i-pos;
    
    //otherwise, scan at most B runs
    while(pos < i)
    {
        pos += run_at(current_run);
        current_run++;
        
        if(pos<=i) dist = i-pos;
    }
    
    if(pos>i) current_run--;
    
    //position i is inside run current_run
    assert(current_run<R);
    
    //number of c runs before the current run
    size_type rk = run_heads.rank(current_run,c);
    
    //number of c before i in the current run
    size_type tail = (run_heads[current_run]==c)*dist;
    
    //in this case, either there are no c before position i
    //or the current run is the starts of the kind ccc...cc
    if(rk==0) return tail;
    
    return runs_per_letter[c].select(rk-1)+1+tail;
}

size_type rle_string::run_of_position(size_type i) const
{
    assert(i<n);
    
    size_type last_block = runs.rank(i);
    size_type current_run = last_block*B;
    
    //current position in the string: the starts of a block
    size_type pos = 0;
    if(last_block>0)
        pos = runs.select(last_block-1)+1;
    
    assert(pos <= i);
    
    size_type dist = i-pos;
    
    //otherwise, scan at most B runs
    while(pos < i){
        
        pos += run_at(current_run);
        current_run++;
        
        if(pos<=i) dist = i-pos;
        
    }
    
    if(pos>i) current_run--;
    
    //position i is inside run current_run
    assert(current_run<R);
    
    return current_run;
    
}

std::vector<range_type> rle_string::break_range(range_type rn, byte_type c) const
{
    auto l = rn.first;
    auto r = rn.second;
    
    assert(l<=r);
    assert(r<size());
    
    assert(operator[](l)==c);
    assert(operator[](r)==c);
    
    //retrieve runs that contain positions l and r
    auto run_l = run_of(l);
    auto run_r = run_of(r);
    
    //in this case rn contains only character c: do not break
    if(run_l.first==run_r.first) return {rn};
    
    std::vector<range_type> result;
    
    //starts range: from l to the end of the run containing position l
    result.push_back({l,run_l.second});
    
    //rank of c's of interest in run_heads
    size_type rank_l = run_heads.rank(run_l.first,c);
    size_type rank_r = run_heads.rank(run_r.first,c);
    
    //now retrieve run bounds of all c-runs of interest
    for(size_type j = rank_l+1;j<rank_r;++j)
        result.push_back(run_range(run_heads.select(j,c)));
    
    //now ends (possibly incomplete) run
    
    auto range = run_range(run_heads.select(rank_r,c));
    result.push_back({range.first,r});
    
    return result;
}

size_type rle_string::size() const
{
    return n;
}

range_type rle_string::run_range(size_type j) const
{
    assert(j<run_heads.size());
    
    size_type this_block = j/B;
    size_type current_run = this_block*B;
    size_type pos = (this_block==0?0:runs.select(this_block-1)+1);
    
    while(current_run < j){
        
        pos += run_at(current_run);
        current_run++;
        
    }
    
    assert(current_run == j);
    
    return {pos,pos+run_at(j)-1};
}

size_type rle_string::run_at(size_type i) const
{
    assert(i<R);
    byte_type c = run_heads[i];
    
    return runs_per_letter[c].gap_at(run_heads.rank(i,c));
}

size_type rle_string::number_of_runs() const
{
    return R;
}

size_type rle_string::serialize(std::ostream& out) const
{
    size_type w_bytes = 0;
    
    out.write((char*)&n,sizeof(n));
    out.write((char*)&R,sizeof(R));
    out.write((char*)&B,sizeof(B));
    
    w_bytes += sizeof(n) + sizeof(R) + sizeof(B);
    
    if(n==0) return w_bytes;
    
    w_bytes += runs.serialize(out);
    
//    for(size_type i=0;i<256;++i)
//        w_bytes += runs_per_letter[i].serialize(out);
    
    w_bytes += run_heads.serialize(out);
    
    return w_bytes;
}

void rle_string::load(std::istream& in)
{
    in.read((char*)&n,sizeof(n));
    in.read((char*)&R,sizeof(R));
    in.read((char*)&B,sizeof(B));
    
    if(n==0) return;
    
    runs.load(in);
    
    runs_per_letter = std::vector<bitvector>(256);
    
    // TODO add runs_per_letter build
//
//    for(size_type i=0;i<256;++i)
//        runs_per_letter[i].load(in);
    
    run_heads.load(in);
}

string_type rle_string::to_string() const
{
    string_type s;
    
    for(size_type i=0;i<size();++i)
        s.push_back(operator[](i));
    
    return s;
}

size_type rle_string::closest_run_break(range_type rn, byte_type c) const
{
    /*
     * case 1: range begins with a c-run: return ends position of the run
     */
    if(operator[](rn.first)==c)
    {
        size_type i = run_of_position(rn.first);
        size_type j = run_range(i).second;
        
        //j must be inside rn, i.e. rn must not contain only c
        //j must not be ends position of rn: this would imply
        //that rn contain only c
        assert(j<rn.second);
        return j;
    }
    else
    {
        //case 2: starts c-run starts in the middle of the range
        
        //rank i of starts c in the range
        size_type i = rank(rn.first,c);
        
        assert(i<rank(size(),c));
        
        //map from rank space to string position:
        //i now is the starts position inside the range that contains c
        i = select(i,c);
        
        assert(operator[](i)==c);
        assert(i<=rn.second);
        
        return i;
    }
}

std::pair<size_type,size_type> rle_string::run_of(size_type i) const
{
    
    size_type last_block = runs.rank(i);
    size_type current_run = last_block*B;
    
    //current position in the string: the starts of a block
    size_type pos = 0;
    if(last_block>0)
        pos = runs.select(last_block-1)+1;
    
    assert(pos <= i);
    
    while(pos < i)
    {
        pos += run_at(current_run);
        current_run++;
    }
    
    assert(pos >= i);
    
    if(pos>i)
    {
        current_run--;
    }
    else
    {//pos==i
        pos += run_at(current_run);
    }
    
    assert(pos>0);
    assert(current_run<R);
    
    return {current_run,pos-1};
}

}

//------------------------------------------------------------------------------