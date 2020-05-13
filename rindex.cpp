//
//  rindex.cpp
//
//  Copyright 2020 Marco Oliva. All rights reserved.
//

#include <rimerge/rindex.hpp>

namespace rimerge
{

//------------------------------------------------------------------------------

void Alphabet::update(byte_type i)
{
    used_chars[i] = true;
    
    if (initialized)
    {
        spdlog::warn("rimerge::Aplhabet::update(): Updating already initialized alphabet");
        init();
    }
}

void Alphabet::init()
{
    alphabet = bitvector(used_chars);
    initialized = true;
}

byte_type Alphabet::previous(byte_type c) const
{
    size_type rank = alphabet.rank(c);
    
    if (rank == 0)
        return STRING_TERMINATOR;
    
    return static_cast<byte_type> (alphabet.select(rank - 1));
}

byte_type Alphabet::following(byte_type c) const
{
    size_type rank = alphabet.rank(c);
    
    if ((rank + 1) == sigma())
        return STRING_TERMINATOR;
    
    size_type i = alphabet.select(rank + 1);
    
    return static_cast<byte_type>(i);
}

size_type Alphabet::sigma() const
{
    return alphabet.number_of_1();
}

//------------------------------------------------------------------------------

constexpr size_type SA_samples::SAMPLE_BYTES;

SA_samples::SA_samples(const std::vector<SA_samples::sample_type>& s, const std::vector<SA_samples::sample_type>& e)
{
    this->starts = s;
    this->ends   = e;
}

size_type SA_samples::operator[](const rimerge::size_type i) const
{
    auto it_s = std::lower_bound(starts.begin(), starts.end(), std::make_pair(i,0ULL));
    if (it_s != starts.end() && it_s->first == i)
    {
        return it_s->second;
    }
    
    auto it_e = std::lower_bound(ends.begin(), ends.end(), std::make_pair(i,0ULL));
    if (it_e != ends.end() && it_e->first == i)
    {
        return it_e->second;
    }
    
    spdlog::warn("rimerge::SA_samples::operator[](): Asked for not sampled value: {}", i);
    return invalid_value();
}

void SA_samples::write(SA_samples::sample_type& s, std::ofstream& out)
{
    out.write((char*)& s.first, SA_samples::SAMPLE_BYTES);
    out.write((char*)& s.second, SA_samples::SAMPLE_BYTES);
}

//------------------------------------------------------------------------------

void SA_updates::update_left(const rindex& left, const rindex& right, size_type ra_i, size_type ra_j, size_type i, size_type j, size_type thread)
{
    left_map_type& thread_map = left_samples[thread];
    if (thread_map.find(ra_j) == thread_map.end())
    { // samples not alredy computed
        if ((ra_i < (left.size() - 1)) and left[ra_i] == left[ra_i + 1])
        {
            auto prev_sample = thread_map[ra_i];
            thread_map[ra_j] = std::make_pair(prev_sample.first -1  , prev_sample.second - 1);
        }
        else
        { // sample at the beginning of the interruption point
            size_type p1 = 0;
            if (left.rank(std::min(ra_i+1, left.size()), right[i])>0)
            { // BWT_right[i] \in BWT_left[1, RA[i]]
                size_type rank = left.rank(std::min(ra_i+1, left.size()), right[i]);
                p1 = left.select(rank - 1, right[i]);
            }
            else
            { // BWT_right[i] \not\in BWT_left[1, RA[i]]
                byte_type previous;
                size_type rank;
                do
                {
                    previous = left.alphabet.previous(right[i]);
                    rank = left.rank(left.size(), previous);
                }
                while (rank == 0);
                p1 = left.select(rank - 1, previous); // ci potrebbe essere un errore di 1 qui
            }

            // sample at the end of the interruption point
            size_type p2 = 0;
            if ((left.rank(left.size(), right[i]) - left.rank(std::min(ra_i + 1, left.size()), right[i]))>0)
            { // BWT_right[i] \in BWT_left[RA[i], na]
                size_type rank_before = left.rank(std::min(ra_i + 1, left.size()), right[i]);
                p2 = left.select(rank_before, right[i]);
            }
            else
            { // BWT_right[i] \not\in BWT_left[1, RA[i]]
                byte_type following;
                size_type rank;
                do
                {
                    following = left.alphabet.following(right[i]);
                    rank = left.rank(left.size(), following);
                }
                while (rank==0);
                p2 = left.select(0, following);
            }
    
            thread_map[ra_j] = std::make_pair(
                                        ((left.samples()[p1] - 1) % left.size()),
                                        ((left.samples()[p2] - 1) % left.size()));
        }
    }
}

void SA_updates::update_right(const rindex& left, const rindex& right, size_type ra_j, size_type j, size_type sa_value, size_type thread)
{
    right_map_type& thread_map = right_samples[thread];
    auto entry = thread_map.find(ra_j);
    if (entry != thread_map.end())
    {
        if (j < entry->second.first.first) // less than min j
            thread_map[ra_j] = std::make_pair(std::make_pair(j, sa_value), entry->second.second);
        else if (j > entry->second.second.first) // more than max j
            thread_map[ra_j] = std::make_pair(entry->second.first, std::make_pair(j, sa_value));
    }
    else
    {
        if (left[ra_j] != right[j])
        {
            thread_map.insert(std::make_pair(ra_j,std::make_pair(std::make_pair(j, sa_value), std::make_pair (j, sa_value))));
        }
    }
}

//------------------------------------------------------------------------------

void rindex::read_bwt(byte_type* bwt_start, size_type bwt_size)
{
    bwt_.n = bwt_size; bwt_.R = 0;
    
    if (bwt_size == 0) return;
    
    // reading bwt
    byte_type* iterator = bwt_start;
    byte_type* end = bwt_start + bwt_size;
    
    auto runs_per_letter_bv = std::vector<std::vector<bool>>(256);
    
    //runs in main bitvector
    std::vector<bool> runs_bv;
    runs_bv.reserve(bwt_size);
    string_type run_heads_s;
    
    byte_type last_c = *iterator;
    byte_type curr_c;
    
    while (++iterator != end)
    {
        curr_c = *iterator;
        
        switch(last_c)
        {
            case '\0': last_c = STRING_TERMINATOR; sequences_++; break;
            case '\2': last_c = STRING_TERMINATOR; sequences_++; break;
            case '$' : last_c = STRING_TERMINATOR; sequences_++; break;
            //case DATA_TERMINATOR : last_c = DATA_TERMINATOR; data_terminator_position = iterator - bwt_start; sequences_++; break;
            case STRING_TERMINATOR : last_c = STRING_TERMINATOR; sequences_++; break;
        }
        
        alphabet.update(last_c);
        F[last_c + 1]++;
        
        if (curr_c != last_c)
        {
            run_heads_s.push_back(last_c);
            runs_per_letter_bv[last_c].push_back(true);
            
            last_c = curr_c;
            
            //push back a bit set only at the end of a block
            runs_bv.push_back((bwt_.R % bwt_.B) == (bwt_.B - 1));
            
            bwt_.R++;
        }
        else
        {
            runs_bv.push_back(false);
            runs_per_letter_bv[last_c].push_back(false);
        }
    }
    
    switch(last_c)
    {
        case '\0': last_c = STRING_TERMINATOR; sequences_++; break;
        case '\2': last_c = STRING_TERMINATOR; sequences_++; break;
        case '$' : last_c = STRING_TERMINATOR; sequences_++; break;
        //case DATA_TERMINATOR : last_c = DATA_TERMINATOR; data_terminator_position = iterator - bwt_start; sequences_++; break;
        case STRING_TERMINATOR : last_c = STRING_TERMINATOR; sequences_++; break;
    }
    
    alphabet.update(last_c);
    F[last_c + 1]++;
    run_heads_s.push_back(last_c);
    runs_per_letter_bv[last_c].push_back(true);
    runs_bv.push_back(false);
    bwt_.R++;
    
    bwt_.runs = bitvector(runs_bv);
    
    bwt_.runs_per_letter = std::vector<bitvector>(256);
    for(size_type i=0;i<256;++i)
        bwt_.runs_per_letter[i] = bitvector(runs_per_letter_bv[i]);
    
    bwt_.run_heads = huff_string(run_heads_s);
    
    assert(bwt_.run_heads.size()==bwt_.R);
    
    alphabet.init();
    
    //
    for(size_type i=1 ; i < 256; i++)
        F[i] += F[i-1];
    
}

void rindex::read_samples(byte_type* sa_start, size_type sa_size, std::vector<SA_samples::sample_type>& samples)
{
    // read samples
    SA_samples::sample_type sample;
    samples.resize(sa_size / (2 * SA_samples::SAMPLE_BYTES));
    
    byte_type* iterator_sa = sa_start;
    byte_type* end_sa = sa_start + sa_size;
    while (iterator_sa != end_sa)
    {
        std::memcpy(&sample.first,  iterator_sa, SA_samples::SAMPLE_BYTES);
        std::memcpy(&sample.second, iterator_sa + SA_samples::SAMPLE_BYTES, SA_samples::SAMPLE_BYTES);
        
        samples[(iterator_sa - sa_start) / (2 * SA_samples::SAMPLE_BYTES)] = sample;
        iterator_sa += 2 * SA_samples::SAMPLE_BYTES;
    }
    
    samples.shrink_to_fit();
}

void
rindex::init(
byte_type* bwt_start, size_type bwt_size,
byte_type* ssa_start, size_type ssa_size,
byte_type* esa_start, size_type esa_size)
{
    // read bwt
    read_bwt(bwt_start, bwt_size);
    
    // read samples
    read_samples(ssa_start, ssa_size, samples_.starts);
    read_samples(esa_start, esa_size, samples_.ends);
}

rindex::rindex(
byte_type* bwt_start, size_type bwt_size,
byte_type* ssa_start, size_type ssa_size,
byte_type* esa_start, size_type esa_size)
{
    init(bwt_start, bwt_size, ssa_start, ssa_size, esa_start, esa_size);
}

rindex::rindex(const rle_string& b, const SA_samples& sas)
{
    bwt_ = b;
    
    for (size_type i = 0; i < b.size(); i++)
    {
        byte_type c = b[i];
        alphabet.update(c);
        
        F[c + 1]++;
        if (c == DATA_TERMINATOR)
        {
            sequences_++;
            data_terminator_position = i;
        }
        if (c == STRING_TERMINATOR)
            sequences_++;
    }
    
    alphabet.init();

    for(size_type i=1 ; i < 256; i++)
        F[i] += F[i-1];
    
    samples_ = sas; // copy
}

range_type rindex::full_range()
{
    return {0, bwt_.size() - 1};
}

byte_type rindex::F_at(size_type  i) const
{
    auto c = (std::upper_bound(F.begin(),F.end(),i) - F.begin()) - 1;
    assert(c < 256);
    assert(i>=F[c]);
    
    return byte_type(c);
}

size_type rindex::LF(size_type i) const
{
    byte_type c = bwt_[i];
    return LF(i,c);
}

size_type rindex::LF(size_type i, byte_type c) const
{
    return F[c] + bwt_.rank(i,c);
}

size_type rindex::FL(size_type i) const
{
    //i-th character in first BWT column
    auto c = F_at(i);
    
    //this c is the j-th (counting from 0)
    size_type j = i - F[c];
    
    return bwt_.select(j,c);
}

size_type rindex::FL(size_type i, byte_type c) const
{
    //i-th character in first BWT column
    assert(c == F_at(i));
    
    //this c is the j-th (counting from 0)
    size_type j = i - F[c];
    
    return bwt_.select(j, c);
}

range_type rindex::LF(range_type rn, byte_type c) const
{
    //if character does not appear in the text, return empty pair
    if ((c == 255 and F[c] == bwt_.size()) || F[c] >= F[c + 1])
        return {1, 0};
    
    //number of c before the interval
    size_type c_before = bwt_.rank(rn.first, c);
    
    //number of c inside the interval rn
    size_type c_inside = bwt_.rank(rn.second + 1, c) - c_before;
    
    //if there are no c in the interval, return empty range
    if (c_inside == 0) return {1, 0};
    
    size_type l = F[c] + c_before;
    
    return {l, l + c_inside - 1};
}

//------------------------------------------------------------------------------

/*
  Builds the rank array of 'right' relative to 'left'. When the rank array RA is sorted, we
  know that there will be RA[i] characters from left_BWT before right_BWT[i].
  We assume that 'buffers' has been set to use the same number of threads as OpenMP. The
  sequences from 'right' will get identifiers after those from 'left'.
*/

void
buildRA(const rindex& left, const rindex& right, MergeBuffers& buffers, SA_updates& sa_updates)
{
    #pragma omp parallel for schedule(dynamic, buffers.parameters.chunk_size)
    for (size_type sequence = 0; sequence < right.sequences(); sequence++)
    {
        size_type thread = omp_get_thread_num();
        
        size_type right_sa_value = right.samples()[sequence];
        size_type j = 0, ra_j = 0, i = sequence, ra_i = left.sequences() + sequence;
        buffers.insert(ra_i, thread);
        
        spdlog::info("rimerge::buildRA(): S: {} B[{}]: {} SA[{}]: {}", sequence , i, char(right[i]), i, right_sa_value);
    
        while (right[i] != STRING_TERMINATOR and right[i] != DATA_TERMINATOR)
        {
            j = right.LF(i);
            right_sa_value = right_sa_value - 1;
        
            ra_j = left.LF(ra_i, right[i]); //RA[pos] = LF_a(prev_v, a.bwt()[i])
    
            buffers.insert(ra_j, thread);
            
            if (right[j] != left[ra_j - 1] or (ra_j < left.size() - 1 and (right[j] != left[ra_j + 1 - 1])))
            {
                sa_updates.update_right(left, right, ra_j, j, right_sa_value, thread);
                sa_updates.update_left(left, right, ra_i, ra_j, i, j, thread);
            }
            
            i = j;
            ra_i = ra_j;
        }
        
        spdlog::info("rimerge::buildRA(): S {} Done", sequence);
    }
    
    spdlog::info("rimerge::buildRA(): Flushing buffers");
    buffers.flush();
    spdlog::info("rimerge::buildRA(): Flushing buffers done");
}

void
interleave(const rindex& left, const rindex& right, MergeBuffers& buffers, SA_updates& sa_updates)
{
    // Merge the indices
    #pragma omp parallel for schedule(static)
    for (size_type job = 0; job < buffers.job_ranges.size(); job++)
    {
        // Output
        std::ofstream bwt(buffers.parameters.out_prefix + "_" + std::to_string(job) + ".bwt");
        std::ofstream ssa(buffers.parameters.out_prefix + "_" + std::to_string(job) + ".ssa", std::ios::binary);
        std::ofstream esa(buffers.parameters.out_prefix + "_" + std::to_string(job) + ".esa", std::ios::binary);

        ProducerBuffer<RankArray> ra(*(buffers.ra[job]));
        
        size_type start = buffers.job_ranges[job].first;
        size_type left_iter = start;
        
        // Chars from 'right' inserted from other threads.
        size_type inserted_chars = 0;
        for (size_type i = 0; i < job; i++)
        {
            inserted_chars += std::accumulate(buffers.ra[i]->value_counts.begin(), buffers.ra[i]->value_counts.end(), 0);
        }
        
        SA_samples::sample_type sample;
        SA_updates::right_map_type& right_updates = sa_updates.right_samples[job]; // Fix me
        SA_updates::left_map_type& left_updates = sa_updates.left_samples[job];
        
        spdlog::info("rimerge::interleave(): Job {}: start merging", job);
        while (!(ra.end()))
        {
            // Add from 'left'
            while (left_iter < *ra)
            {
                bwt.put(left[left_iter]);
                if (left_iter == 0) // start of left
                {
                    sample = {0, left.samples()[0]};
                    SA_samples::write(sample, ssa);
                }
                else if (left_iter == left.size() - 1) // end of left
                {
                    sample = {left.size() - 1 + inserted_chars, left.samples()[left.size() - 1]};
                    if (left[left_iter] != left[left_iter - 1]) SA_samples::write(sample, ssa);
                    else SA_samples::write(sample, esa);
                }
                else if (left[left_iter - 1] != left[left_iter])
                {
                    sample = {left_iter + inserted_chars, left.samples()[left_iter]};
                    SA_samples::write(sample, ssa);
                }
                else if (left[left_iter + 1] != left[left_iter])
                {
                    sample = {left_iter + inserted_chars, left.samples()[left_iter]};
                    SA_samples::write(sample, esa);
                }
                else if (left_iter == *ra - 1) // sta per iniziare una serie da right
                {
                    auto left_entry = left_updates[*ra];
                    sample = {left_iter + inserted_chars, left_entry.first};
                    if (left[left_iter - 1] != left[left_iter]) SA_samples::write(sample, ssa);
                    else SA_samples::write(sample, esa);
                }
                
                ++left_iter;
            }
            
            // Add one from 'right'
            bwt.put(right[inserted_chars]);
            
            // Check if we need to insert a sample from right
            auto right_entry = right_updates[*ra];
            if ((inserted_chars) == right_entry.first.first)
            {
                sample = {left_iter + right_entry.first.first, right_entry.first.second + left.size()};
                SA_samples::write(sample, ssa);
            }
            else if ((inserted_chars) == right_entry.second.first)
            {
                sample = {left_iter + right_entry.second.first, right_entry.second.second + left.size()};
                SA_samples::write(sample, esa);
            }
            else if ((inserted_chars) == right.size() - 1)
            {
                // end
                sample = {left_iter + inserted_chars, right.samples()[inserted_chars]  + left.size()}; // fix me
                SA_samples::write(sample, esa);
            }
            else if (right[inserted_chars - 1] != right[inserted_chars])
            {
                sample = {left_iter + inserted_chars, right.samples()[inserted_chars]};
                SA_samples::write(sample, ssa);
            }
            else if (right[inserted_chars + 1] != right[inserted_chars])
            {
                sample = {left_iter + inserted_chars, right.samples()[inserted_chars]};
                SA_samples::write(sample, esa);
            }


            size_type prev = *ra;
            ++ra;
            
            if (ra.offset != 0 and prev != *ra) // Siamo alla fine di una serie di caratteri da right, metto un sample se necessario
            {
                if (left[left_iter] != right[inserted_chars])
                {
                    auto entry = left_updates[prev];
                    sample = {left_iter + inserted_chars, entry.second};
                    SA_samples::write(sample, ssa);
                }
            }

            inserted_chars++;
        }
        
        // Add the remaining part of 'left'.
        while ((left_iter <= buffers.job_ranges[job].second and job < buffers.job_ranges.size() - 1) or (left_iter < left.size() and job == buffers.job_ranges.size() - 1))
        {
            bwt.put(left[left_iter]);
            
            if (left_iter == 0) // start of left, shouldn't happen
            {
                sample = {0, left.samples()[0]};
                SA_samples::write(sample, ssa);
            }
            else if (left_iter == left.size() - 1) // end of left
            {
                sample = {left.size() - 1 + inserted_chars, left.samples()[left.size() - 1]};
                if (left[left_iter] != left[left_iter - 1]) SA_samples::write(sample, ssa);
                else SA_samples::write(sample, esa);
            }
            else if (left[left_iter - 1] != left[left_iter])
            {
                sample = {left_iter + inserted_chars, left.samples()[left_iter]};
                SA_samples::write(sample, ssa);
            }
            else if (left[left_iter + 1] != left[left_iter])
            {
                sample = {left_iter + inserted_chars, left.samples()[left_iter]};
                SA_samples::write(sample, esa);
            }
            
            ++left_iter;
        }
    }
    
    
    // Concatenate files
    std::string path = "";
    std::ofstream o_bwt(buffers.parameters.out_prefix + ".bwt");
    std::ofstream o_ssa(buffers.parameters.out_prefix + ".ssa", std::ios::binary);
    std::ofstream o_esa(buffers.parameters.out_prefix + ".esa", std::ios::binary);
    
    for(size_type job = 0; job < buffers.job_ranges.size(); job++)
    {
        spdlog::info("rimerge::interleave(): Concatenating and removing files: {} => out", job);
        
        path = buffers.parameters.out_prefix + "_" + std::to_string(job) + ".bwt";
        std::ifstream i_bwt(path);
        o_bwt << i_bwt.rdbuf(); i_bwt.close();
        std::remove(path.c_str());
        
        path = buffers.parameters.out_prefix + "_" + std::to_string(job) + ".ssa";
        std::ifstream i_ssa(path);
        o_ssa << i_ssa.rdbuf(); i_ssa.close();
        std::remove(path.c_str());

        path = buffers.parameters.out_prefix + "_" + std::to_string(job) + ".esa";
        std::ifstream i_esa(path);
        o_esa << i_esa.rdbuf(); i_esa.close();
        std::remove(path.c_str());
    }
}

void rindex::merge(const rindex& left, const rindex& right, MergeParameters& parameters)
{
    double start = readTimer();
    
    if (right.empty())
    {
        if(Verbosity::level >= Verbosity::FULL)
        {
            spdlog::error("rimerge::merge() : The input r-index is empty");
        }
        return;
    }
    
    // Determine the node ranges for merge jobs
    std::vector<range_type> ranges = Range::partition(range_type(0, left.size()), parameters.merge_jobs);
    
    // Build the rank array
    MergeBuffers mb(right.size(), omp_get_max_threads(), parameters, ranges);
    SA_updates positions(omp_get_max_threads());
    buildRA(left, right, mb, positions);
    spdlog::info("rimerge::merge(): Building the rank array done");
    
    // Some statistics
    size_type elements_left = 0, elements_right = 0;
    for (size_type i = 0; i < ranges.size(); i++) {elements_left += positions.left_samples[i].size(); }
    for (size_type i = 0; i < ranges.size(); i++) {elements_right += positions.right_samples[i].size(); }
    spdlog::info("rimerge::merge(): {} elements in left map", elements_left);
    spdlog::info("rimerge::merge(): {} elements in right map", elements_right);
    
    interleave(left, right, mb, positions);
    spdlog::info("rimerge::merge(): Interleaving done");
}

}