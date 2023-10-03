//
//  RIndexRLE.cpp
//
// Copyright (c) Boucher Lab. All rights reserved.
// Licensed under the GNU license. See LICENSE file in the repository root for full license information.

#include <rimerge/bwtmerge.hpp>
#include <rimerge/support.hpp>
#include <rimerge/r-index-rle.hpp>

namespace rimerge
{

//------------------------------------------------------------------------------

const char RIndexRLE::name_[] = "RIndexRLE";
const char* RIndexRLE::name() const { return name_; }

//------------------------------------------------------------------------------

void
RIndexRLE::read_bwt(std::string& bwt_path)
{
    bwt_.load(bwt_path);
    
    // initialize F column and alphabet
    for (size_type i = 0; i < 256; i++)
    {
        size_type size = bwt_.runs_per_letter[i].size();
        F[i] = size;
        if (size > 0) { alphabet.update((byte_type) i); }
        if ((byte_type) i == STRING_TERMINATOR) { sequences_ = size; }
    }
    for (size_type i = F.size() - 1; i > 0; i--) F[i] = F[i-1];
    F[0] = 0;
    for(size_type i = 1; i < F.size(); i++) F[i] += F[i - 1];
    
    if (bwt_.size() != 0 and sequences_ == 0) { sequences_ = 1; }
    
    alphabet.init();
}

sample_genre
RIndexRLE::its(size_type i) const
{
    sample_genre out = sample_genre::NOT;
    if (i == 0) out = sample_genre::START;
    else if (i < sequences_) out = static_cast<sample_genre>(out | sample_genre::START_END);
    else if (i == bwt_.size() - 1) out = sample_genre::END;
    else
    {
        if ((i > 0) and (bwt_.runs[i - 1])) out = static_cast<sample_genre>(out | sample_genre::START);
        if (bwt_.runs[i]) out = static_cast<sample_genre>(out | sample_genre::END);
    }
    return out;
}

range_type
RIndexRLE::full_range() const
{
    return {0, bwt_.size() - 1};
}

byte_type
RIndexRLE::F_at(size_type  i) const
{
    auto c = (std::upper_bound(F.begin(),F.end(),i) - F.begin()) - 1;
    assert(c < 256);
    assert(i>=F[c]);
    
    return byte_type(c);
}

size_type
RIndexRLE::LF(size_type i) const
{
    byte_type c = bwt_[i];
    return LF(i,c);
}

size_type
RIndexRLE::LF(size_type i, byte_type c) const
{
    return F[c] + bwt_.rank(i,c);
}

size_type
RIndexRLE::FL(size_type i) const
{
    //i-th character in first BWT column
    auto c = F_at(i);
    
    //this c is the j-th (counting from 0)
    size_type j = i - F[c];
    
    return bwt_.select(j,c);
}

size_type
RIndexRLE::FL(size_type i, byte_type c) const
{
    //i-th character in first BWT column
    assert(c == F_at(i));
    
    //this c is the j-th (counting from 0)
    size_type j = i - F[c];
    
    return bwt_.select(j, c);
}

range_type
RIndexRLE::LF(range_type rn, byte_type c) const
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

string_type
RIndexRLE::get_sequence(size_type i) const
{
    string_type out; out.push_back(bwt_[i]);
    size_type j = LF(i);
    while (bwt_[j] != DATA_TERMINATOR and bwt_[j] != STRING_TERMINATOR)
    {
        out.push_back(bwt_[j]);
        j = LF(j);
    }
    std::reverse(out.begin(), out.end());
    return out;
}

//------------------------------------------------------------------------------
void
RIndexRLE::SamplesMergerRLE::operator()(size_type index, RLEString::RunCache& right_cache, RLEString::RunCache& left_cache, bool FL, size_type inserting_index, size_type ra_value, size_type prev_ra_value, size_type next_ra_value)
{
    using SG = sample_genre;
    int thread_id = omp_get_thread_num();
    
    SA_samples::sample_type sample;
    
    if (inserting_index == (left.sequences() + right.sequences()))
    {
        if (FL) { sample = {inserting_index, left.samples()[index]}; LFL = true; LLI = index; }
        else { sample = {inserting_index, right.samples()[index] + left.size()}; LFL = false; LRI = index; }
        SA_samples::write(sample, saes);
    }
    if (FL and index < left.sequences()) // first samples from left
    {
        sample = {inserting_index, left.samples()[index]};
        SA_samples::write(sample, saes);
        LFL = true; LLI = index;
    }
    else if (not FL and index < right.sequences()) // first samples from right
    {
        sample = {inserting_index, right.samples()[index] + left.size()};
        SA_samples::write(sample, saes);
        LFL = false; LRI = index;
    }
    else if (FL and index == (left.size() - 1))
    {
        sample = {inserting_index, left.samples()[index]};
        SA_samples::write(sample, saes);
        LFL = true; LLI = index;
    }
    else if (not FL and index == (right.size() - 1))
    {
        sample = {inserting_index, right.samples()[index] + left.size()};
        SA_samples::write(sample, saes);
        LFL = false; LRI = index;
    }
    else if (FL and LFL)
    {
        auto its = left.its(index);
        if (its & SG::START_END)
        {
            if ( (its & SG::START) or ((its & SG::END) and (index != (ra_value - 1))) )
            {
                sample = {inserting_index, left.samples()[index]};
                SA_samples::write(sample, saes);
            }
            else if ( (index == (ra_value - 1)) and (left_cache[index] != right_cache[LRI + 1]) )
            {
                sample = {inserting_index, left.samples()[index]};
                SA_samples::write(sample, saes);
            }
        }
        else if ( (index == (ra_value - 1)) and (left_cache[index] != right_cache[LRI + 1]) )
        {
            auto down = updates.find_left(ra_value - 1);
            counter_left++;
            if (down != updates.end_left())
            {
                sample = {inserting_index, down};
                SA_samples::write(sample, saes);
            }
            else
            {
                spdlog::error("Job: {} Sample missing in left_updates! {} 1 {}", thread_id, ra_value, inserting_index);  std::exit(EXIT_FAILURE);
            }
        }
        LFL = true; LLI = index;
    }
    else if (FL and not LFL)
    {
        auto its = left.its(index);
        if (left_cache[index] != right_cache[LRI])
        {
            auto up = updates.find_left(prev_ra_value);
            if (up != updates.end_left())
            {
                counter_left++;
                sample = {inserting_index, up};
                SA_samples::write(sample, saes);
            }
            else if (its & SG::START_END)
            {
                sample = {inserting_index, left.samples()[index]};
                SA_samples::write(sample, saes);
            }
            else
            {
                spdlog::error("Job: {} Sample missing in left_updates! {} 2 {}", thread_id, prev_ra_value, inserting_index);  std::exit(EXIT_FAILURE);
            }
        }
        else if ((index == ra_value - 1) and left_cache[index] != right_cache[LRI + 1])
        {
            auto down = updates.find_left(ra_value - 1);
            if (down != updates.end_left())
            {
                counter_left++;
                sample = {inserting_index, down};
                SA_samples::write(sample, saes);
            }
            else if (its & SG::START_END)
            {
                sample = {inserting_index, left.samples()[index]};
                SA_samples::write(sample, saes);
            }
            else
            {
                spdlog::error("Job: {} Sample missing in left_updates! {} 3 {}", thread_id, ra_value, inserting_index);  std::exit(EXIT_FAILURE);
            }
        }
        else if ( (index != ra_value - 1) and (its & SG::END) )
        {
            sample = {inserting_index, left.samples()[index]};
            SA_samples::write(sample, saes);
        }
        else if ( (index != (ra_value - 1)) and (its & SG::START) and (left_cache[index] != right_cache[LRI]) )
        {
            sample = {inserting_index, left.samples()[index]};
            SA_samples::write(sample, saes);
        }
        LFL = true; LLI = index;
    }
    else if (not FL and not LFL)
    {
        auto its = right.its(index);
        if (its & SG::START_END)
        {
            if ( (its & SG::START) or ((its & SG::END) and (ra_value == next_ra_value)) )
            {
                sample = {inserting_index, right.samples()[index] + left.size()};
                SA_samples::write(sample, saes);
            }
            else if ( (ra_value != next_ra_value) and (right_cache[index] != left_cache[LLI + 1]) )
            {
                sample = {inserting_index, right.samples()[index] + left.size()};
                SA_samples::write(sample, saes);
            }
        }
        else if ( (ra_value != next_ra_value) and (right_cache[index] != left_cache[LLI + 1]) )
        {
            auto down = updates.find_right_max(ra_value);
            counter_right++;
            if (down != updates.end_right())
            {
                sample = {inserting_index, down.second + left.size()};
                SA_samples::write(sample, saes);
            }
            else if (down == updates.end_right())
            {
                spdlog::error("Job: {} Sample missing in right_updates! {} 4 {}", thread_id, ra_value, inserting_index);  std::exit(EXIT_FAILURE);
            }
        }
        LFL = false; LRI = index;
    }
    else if (not FL and LFL)
    {
        auto its = right.its(index);
        if (its != SG::NOT)
        {
            if (right_cache[index] != left_cache[LLI])
            {
                sample = {inserting_index, right.samples()[index] + left.size()};
                SA_samples::write(sample, saes);
            }
            else if ( (ra_value != next_ra_value) and right_cache[index] != left_cache[LLI + 1] )
            {
                sample = {inserting_index, right.samples()[index] + left.size()};
                SA_samples::write(sample, saes);
            }
            else if ( (ra_value == next_ra_value) and (its & SG::END) )
            {
                sample = {inserting_index, right.samples()[index] + left.size()};
                SA_samples::write(sample, saes);
            }
        }
        else if (right_cache[index] != left_cache[LLI])
        {
            auto up = updates.find_right_min(ra_value);
            counter_right++;
            if (up != updates.end_right())
            {
                sample = {inserting_index, up.second + left.size()};
                SA_samples::write(sample, saes);
            }
            else
            {
                spdlog::error("Job: {} Sample missing in right_updates! {} 5 {}", thread_id, ra_value, inserting_index);  std::exit(EXIT_FAILURE);
            }
        }
        else if ( (ra_value != next_ra_value) and (right_cache[index] != left_cache[LLI + 1]) )
        {
            auto down = updates.find_right_max(ra_value);
            counter_right++;
            if (down != updates.end_right())
            {
                sample = {inserting_index, down.second + left.size()};
                SA_samples::write(sample, saes);
            }
            else
            {
                spdlog::error("Job: {} Sample missing in right_updates! {} 6 {}", thread_id, ra_value, inserting_index);  std::exit(EXIT_FAILURE);
            }
        }
        LFL = false; LRI = index;
    }
}


//------------------------------------------------------------------------------

std::pair<size_type, size_type>
RIndexRLE::SAUpdatesRLE::update_left(const RIndex<RIndexRLE, RLEString>& left, const RIndex<RIndexRLE, RLEString>& right, size_type ra_i, size_type ra_j, std::pair<size_type, size_type> prev_samples, size_type i, size_type j, size_type thread, RLEString::Accessor& right_accessor, RLEString::Accessor& left_accessor)
{
    left_map_type& thread_map = left_samples[thread];
    
    if (thread_map.find(ra_j) == thread_map.end() or thread_map.find(ra_j - 1) == thread_map.end())
    { // samples not alredy computed
        if (((ra_i < (left.size() /*- 1*/)) and left_accessor[ra_i - 1] == left_accessor[ra_i])
        and right_accessor[i] == left_accessor[ra_i])
        {
            return  std::make_pair(prev_samples.first - 1, prev_samples.second -1);
        }
        else
        { // sample at the beginning of the interruption point
            size_type p1 = 0;
            if (left.rank(std::min(ra_i, left.size()), right_accessor[i])>0)
            { // BWT_right[i] \in BWT_left[1, RA[i]]
                size_type rank = left.rank(std::min(ra_i, left.size()), right_accessor[i]);
                p1 = left.select(rank - 1, right_accessor[i]);
            }
            else
            { // BWT_right[i] \not\in BWT_left[1, RA[i]]
                byte_type previous;
                size_type rank;
                do
                {
                    previous = left.alphabet.previous(right_accessor[i]);
                    rank = left.rank(left.size(), previous);
                }
                while (rank == 0);
                p1 = left.select(rank - 1, previous); // ci potrebbe essere un errore di 1 qui
            }
            
            // sample at the end of the interruption point
            size_type p2 = 0;
            if ((left.rank(left.size(), right_accessor[i]) - left.rank(std::min(ra_i, left.size()), right_accessor[i]))>0)
            { // BWT_right[i] \in BWT_left[RA[i], na]
                size_type rank_before = left.rank(std::min(ra_i, left.size()), right_accessor[i]);
                p2 = left.select(rank_before, right_accessor[i]);
            }
            else
            { // BWT_right[i] \not\in BWT_left[1, RA[i]]
                byte_type following;
                size_type rank;
                do
                {
                    following = left.alphabet.following(right_accessor[i]);
                    rank = left.rank(left.size(), following);
                }
                while (rank==0);
                p2 = left.select(0, following);
            }
            return std::make_pair(((left.samples()[p1] - 1) % left.size()), ((left.samples()[p2] - 1) % left.size()));
        }
    }
    else //sample in map
    {
        return std::make_pair(thread_map.find(ra_j - 1)->second, thread_map.find(ra_j)->second);
    }
}

void
RIndexRLE::SAUpdatesRLE::update_right_min(const RIndex<RIndexRLE, RLEString>& left, const RIndex<RIndexRLE, RLEString>& right, size_type ra_j, size_type j, size_type sa_value, size_type thread, RLEString::Accessor& right_accessor, RLEString::Accessor& left_accessor)
{
    auto& thread_map_min = right_samples[thread].first;
    
    // Min
    auto entry_min = thread_map_min.find(ra_j);
    if (entry_min != thread_map_min.end())
    {
        if ( (j < entry_min->second.first) and (left_accessor[ra_j - 1] == right_accessor[j]) )
        { thread_map_min.erase(entry_min); }
        else if (j < entry_min->second.first)
        {thread_map_min[ra_j] = std::make_pair(j, sa_value);}
    }
    else if (right.its(j) == sample_genre::NOT and left_accessor[ra_j - 1] != right_accessor[j])
    { thread_map_min.insert(std::make_pair(ra_j, std::make_pair(j, sa_value))); }
}

void
RIndexRLE::SAUpdatesRLE::update_right_max(const RIndex<RIndexRLE, RLEString>& left, const RIndex<RIndexRLE, RLEString>& right, size_type ra_j, size_type j, size_type sa_value, size_type thread, RLEString::Accessor& right_accessor, RLEString::Accessor& left_accessor)
{
    auto& thread_map_max = right_samples[thread].second;
    
    // Max
    auto entry_max = thread_map_max.find(ra_j);
    if (entry_max != thread_map_max.end())
    {
        if ( (j > entry_max->second.first) and (ra_j < left.size() and left_accessor[ra_j] == right_accessor[j]) )
        { thread_map_max.erase(entry_max); }
        else if (j > entry_max->second.first)
        { thread_map_max[ra_j] = std::make_pair(j, sa_value); }
    }
    else if (ra_j >= left.size()) { thread_map_max.insert(std::make_pair(ra_j, std::make_pair(j, sa_value))); }
    else if (right.its(j) == sample_genre::NOT and left_accessor[ra_j] != right_accessor[j])
    { thread_map_max.insert(std::make_pair(ra_j, std::make_pair(j, sa_value))); }
}

//------------------------------------------------------------------------------

/*
  Builds the rank array of 'right' relative to 'left'. When the rank array RA is sorted, we
  know that there will be RA[i] characters from left_BWT before right_BWT[i].
  We assume that 'buffers' has been set to use the same number of threads as OpenMP. The
  sequences from 'right' will get identifiers after those from 'left'.
*/

void
buildRA(const RIndex<RIndexRLE, RLEString>& left, const RIndex<RIndexRLE, RLEString>& right, MergeBuffers& buffers, RIndexRLE::SAUpdatesRLE& sa_updates)
{
    double start = readTimer();
    
    #pragma omp parallel for schedule(dynamic, buffers.parameters.chunk_size)
    for (size_type sequence = 0; sequence < right.sequences(); sequence++)
    {
        size_type thread = omp_get_thread_num();
        
        RLEString::Accessor right_accessor(right.bwt());
        RLEString::Accessor left_accessor(left.bwt());
        
        size_type right_sa_value = right.samples()[sequence];
        size_type j = 0, ra_j = 0, i = sequence, ra_i = left.sequences();
        buffers.insert(ra_i, thread);
        
        std::pair<size_type, size_type> prev_samples = std::make_pair(left.samples()[ra_i -1], left.samples()[ra_i]);
        sa_updates.left_samples[thread].insert(std::make_pair((ra_i - 1), prev_samples.first));
        sa_updates.left_samples[thread].insert(std::make_pair((ra_i), prev_samples.second));
        spdlog::info("rimerge::buildRA(): S: {} B[{}]: {} SA[{}]: {}", sequence , i, char(right_accessor[i]), i, right_sa_value);
        
        while (right_accessor[i] != STRING_TERMINATOR and right_accessor[i] != DATA_TERMINATOR)
        {
            j = right.LF(i, right_accessor[i]);
            right_sa_value = right_sa_value - 1;
            
            ra_j = left.LF(ra_i, right_accessor[i]); //RA[pos] = LF_a(prev_v, a.bwt()[i])
            
            buffers.insert(ra_j, thread);
            
            prev_samples = sa_updates.update_left(left, right, ra_i, ra_j, prev_samples, i, j, thread, right_accessor, left_accessor);
            
            // Breaking a run
            if ((right_accessor[j] != left_accessor[ra_j - 1]) and (left.its(ra_j - 1) == sample_genre::NOT))
            {
                sa_updates.left_samples[thread].insert(std::make_pair(ra_j - 1, prev_samples.first));
                sa_updates.left_samples[thread].insert(std::make_pair(ra_j, prev_samples.second));
            }
            
            if (((ra_j <= left.size() - 1) and (right_accessor[j] != left_accessor[ra_j - 1 + 1]))
            and (left.its(ra_j - 1 + 1) == sample_genre::NOT))
            {
                sa_updates.left_samples[thread].insert(std::make_pair(ra_j - 1, prev_samples.first));
                sa_updates.left_samples[thread].insert(std::make_pair(ra_j, prev_samples.second));
            }
            
            sa_updates.update_right_min(left, right, ra_j, j, right_sa_value, thread, right_accessor, left_accessor);
            sa_updates.update_right_max(left, right, ra_j, j, right_sa_value, thread, right_accessor, left_accessor);
            
            i = j;
            ra_i = ra_j;
        }
        
        spdlog::debug("rimerge::buildRA(): S {} Done", sequence);
    }
    
    buffers.flush();
    
    if(Verbosity::level >= Verbosity::BASIC)
    {
        double seconds = readTimer() - start;
        spdlog::info("Built Rank Array in {} seconds", seconds);
    }
}

void
interleave(const RIndex<RIndexRLE, RLEString>& left, const RIndex<RIndexRLE, RLEString>& right, MergeBuffers& buffers, RIndexRLE::SAUpdatesRLE& sa_updates)
{
    double start = readTimer();
    
    std::string out_path = buffers.parameters.out_prefix + "/bwt.rle";
    RLEString::RLEncoderMerger encoders(out_path, buffers.job_ranges.size());
    
    // Merge the indices
    #pragma omp parallel for schedule(static)
    for (size_type job = 0; job < buffers.job_ranges.size(); job++)
    {
        // Output
        std::string base_path(buffers.parameters.out_prefix + "/tmp_" + std::to_string(job));
        std::ofstream saes(base_path + ".saes", std::ios::binary);
        if (not saes.is_open()) { spdlog::error("Can't open tmp file: {}", base_path + ".saes"); std::exit(EXIT_FAILURE); }
        
        ProducerBuffer<RankArray> ra(*(buffers.ra[job]));
        
        // Set up left iterator
        size_type left_iter = job == 0 ? 0 : buffers.job_ranges[job - 1].second;
        
        // Chars from 'right' inserted from other threads. Also find the last non empty ra
        size_type right_iter = 0, last_non_empty_ra = 0;
        for (size_type i = 0; i < job; i++)
        {
            size_type value_counts = buffers.range_values[i];
            right_iter += value_counts;
            
            if (value_counts != 0) { last_non_empty_ra = i; }
        }
        
        // Check how many ra values in current range,
        size_type current_ra_values = buffers.range_values[job];;
        spdlog::info("Job {} RA values in current range: {}", job, current_ra_values);
        
        // Set up this job
        size_type prev_ra = 0, next_ra = 0, curr_ra = 0;
        
        RIndexRLE::SamplesMergerRLE sample_merger(right, left, &saes, sa_updates);
        sample_merger.set_LLI(left_iter == 0 ? 0 : left_iter - 1);
        sample_merger.set_LRI(right_iter == 0 ? 0 : right_iter - 1);
    
        if (job != 0)
        {
            prev_ra = buffers.max_values[last_non_empty_ra];
            if (prev_ra == buffers.job_ranges[job - 1].second) { sample_merger.set_LFL(false); }
            else { sample_merger.set_LFL(true); }
        }
    
        RLEString::RunCache right_cache(right.bwt());
        RLEString::RunCache left_cache(left.bwt());
    
        RLEString::RLEncoder& rle_encoder = encoders.get_encoder(job);
        
        bool tok = true;
        curr_ra = *ra; ++ra; next_ra = *ra;
    
        spdlog::info("Job: {} Left Range: [{},{}] left_iter: {} right_iter: {}, curr_ra: {}, prev_ra: {}",
                     job, buffers.job_ranges[job].first,
                     buffers.job_ranges[job].second, left_iter, right_iter, curr_ra, prev_ra);
        
        while (curr_ra != invalid_value())
        {
            // Add from 'left'
            while (left_iter < curr_ra)
            {
                rle_encoder(left_cache[left_iter]);
                sample_merger(left_iter, right_cache, left_cache, true, left_iter + right_iter, curr_ra, prev_ra, next_ra);
                left_iter++;
            }
            
            // Add one from 'right'
            rle_encoder(right_cache[right_iter]);
            sample_merger(right_iter, right_cache, left_cache, false, left_iter + right_iter, curr_ra, prev_ra, next_ra);
            
            prev_ra = curr_ra;
            if (not tok) { curr_ra = invalid_value(); }
            else { curr_ra = next_ra; ++ra; }
            
            size_type next_range = job + 1;
            while (ra.end() and next_range < buffers.job_ranges.size())
            {
                if ( buffers.range_values[next_range] != 0 ) { break; }
                else { next_range++; }
            }
            
            if (ra.end() and next_range < buffers.job_ranges.size() and tok) { next_ra = buffers.min_values[next_range]; tok = false; }
            else next_ra = *ra;
            
            right_iter++;
        }
        
        // Add the remaining part of 'left'.
        while ((left_iter < buffers.job_ranges[job].second and job < buffers.job_ranges.size() - 1) or (left_iter < left.size() and job == buffers.job_ranges.size() - 1))
        {
            rle_encoder(left_cache[left_iter]);
            sample_merger(left_iter, right_cache, left_cache, true, left_iter + right_iter, curr_ra, prev_ra, next_ra);
            ++left_iter;
        }
        
        spdlog::info("Last left inserted by {}: {}", job, left_iter - 1);
        
        // Close files
        saes.close();
    }
    
    if(Verbosity::level >= Verbosity::BASIC)
    {
        double seconds = readTimer() - start;
        spdlog::info("Completed interleave in {} seconds", seconds);
    }
    
    // Concatenate files
    // Samples
    start = readTimer();
    
    std::string o_saes(buffers.parameters.out_prefix + "/samples.saes");
    std::remove(o_saes.c_str());
    std::ofstream out_saes_stream(o_saes);
    
    for (size_type job = 0; job < buffers.job_ranges.size(); job++)
    {
        std::string path = buffers.parameters.out_prefix + "/tmp_" + std::to_string(job) + ".saes";
        concatenate_file<byte_type>(out_saes_stream, path);
        std::remove(path.c_str());
    }
    
    if (Verbosity::level >= Verbosity::BASIC)
    {
        double seconds = readTimer() - start;
        spdlog::info("Completed sample merge in {} seconds", seconds);
    }
}

void
RIndexRLE::merge(const RIndex<RIndexRLE, RLEString>& left, const RIndex<RIndexRLE, RLEString>& right, MergeParameters& parameters)
{
    double start = readTimer();
    
    if (right.empty())
    {
        if (Verbosity::level >= Verbosity::FULL)
        {
            spdlog::error("rimerge::merge() : The input r-index is empty");
        }
        return;
    }
    
    // Compute node ranges for merge jobs
    std::vector<range_type> ranges = Range::partition(range_type(0, left.size()), parameters.merge_jobs);
    
    // Build the rank array
    omp_set_num_threads(parameters.search_jobs);
    MergeBuffers mb(right.size(), parameters.search_jobs, parameters, ranges);
    RIndexRLE::SAUpdatesRLE positions(parameters.search_jobs);
    buildRA(left, right, mb, positions);
    spdlog::info("rimerge::merge(): Building the rank array done");
    
    interleave(left, right, mb, positions);
    spdlog::info("rimerge::merge(): Interleaving done");
    
    if(Verbosity::level >= Verbosity::BASIC)
    {
        double seconds = readTimer() - start;
        spdlog::info("Completed merge of {} sequences in {} seconds", right.sequences(), seconds);
    }
}

bool
RIndexRLE::check(const RIndex<RIndexRLE, RLEString>& index, std::tuple<std::vector<size_type>, std::vector<size_type>, std::vector<size_type>>& errors)
{
    size_type mask = ~0ULL; mask = mask >> ((8 - SA_samples::SAMPLE_BYTES) * 8);
    
    std::vector<range_type> ranges = Range::partition(range_type(0, index.size() - 1), omp_get_max_threads());
    std::vector<std::vector<size_type>> missing(ranges.size()) , unnecessary(ranges.size()), invalid(ranges.size());
    
    #pragma omp parallel for
    for (size_type job = 0; job < ranges.size(); job++)
    {
        size_type it  = ranges[job].first, end = ranges[job].second;
        while (it <= end)
        {
            if (index.its(it) != sample_genre::NOT)
            {
                if (index.samples()[it] == rimerge::invalid_value()) { missing[job].push_back(it); }
                if (index.samples()[it] >= (mask - 1000)) { invalid[job].push_back(it); }
            }
            else if (index.samples()[it] != rimerge::invalid_value())
            {
                if (it > index.sequences()) { unnecessary[job].push_back(it); }
            }
            it++;
        }
    }
    
    for (auto& job_vector : missing)
        std::get<0>(errors).insert(std::get<0>(errors).end(), job_vector.begin(), job_vector.end());
    for (auto& job_vector : unnecessary)
        std::get<1>(errors).insert(std::get<1>(errors).end(), job_vector.begin(), job_vector.end());
    for (auto& job_vector : invalid)
        std::get<2>(errors).insert(std::get<2>(errors).end(), job_vector.begin(), job_vector.end());
    
    return std::get<0>(errors).size() == 0 and std::get<2>(errors).size() == 0;
}

std::size_t
RIndexRLE::check_sa_values(const RIndex<RIndexRLE, RLEString>& index)
{
    std::mutex num_of_errors_mtx;
    std::size_t num_of_errors;
    
    #pragma omp parallel for
    for (rimerge::size_type seq = 0; seq < index.sequences(); seq++)
    {
        rimerge::size_type pos = seq, sa_value = index.samples()[seq];
        spdlog::info("SA[{}] = {}", seq, sa_value);
        std::size_t per_thread_errors = 0;
        while (index.bwt()[pos] != rimerge::STRING_TERMINATOR)
        {
            if (index.its(pos) and index.samples()[pos] != sa_value) { per_thread_errors++; }
            
            pos = index.LF(pos);
            sa_value -= 1;
        }
        
        // update global errors
        {
            std::lock_guard<std::mutex> lock(num_of_errors_mtx);
            num_of_errors += per_thread_errors;
        }
    }
    
    return num_of_errors;
}

}