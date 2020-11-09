//
//  r-index.hpp
//
//  Copyright 2020 Marco Oliva. All rights reserved.
//

#ifndef r_index_hpp
#define r_index_hpp

#include <rimerge/inherit.hpp>
#include <rimerge/utils.hpp>
#include <rimerge/alphabet.hpp>
#include <rimerge/sa-samples.hpp>
#include <rimerge/support.hpp>

namespace rimerge
{

template <typename Derived, typename BWT>
class RIndex : public inherit<Derived>
{

public:

//------------------------------------------------------------------------------
    
    class SAUpdates
    {
    public:
        
        //                 RA[j]
        using left_type  = size_type;
        //                           j           j-th sa
        using right_type = std::pair<size_type, size_type>;
        
        //                              RA value
        using left_map_type = std::map<size_type, left_type>;
        //                              RA value
        using right_map_type = std::map<size_type, right_type>;
        
        std::vector<left_map_type> left_samples;
        //                      min j       max j
        std::vector<std::pair<right_map_type,right_map_type>> right_samples;
        
        right_type end_right_;
        left_type   end_left_;
        
    protected:
    
        SAUpdates(size_type threads) : left_samples(threads), right_samples(threads), end_left_(invalid_value()), end_right_(invalid_value(), invalid_value())
        {
            spdlog::info("Creating {} left maps and {} right maps", left_samples.size(), right_samples.size());
        }
        
    public:
    
        const right_type& end_right() const {return this->end_right_; }
        const left_type& end_left() const {return this->end_left_; }
        
        right_type find_right_min(size_type ra_value) const
        {
            std::vector<right_type> results;
    
            // Get results computed by all threads
            for (size_type i = 0; i < right_samples.size(); i++)
            {
                auto entry = right_samples[i].first.find(ra_value);
                if ( entry != right_samples[i].first.end())
                    results.push_back(entry->second);
            }
    
            // Combine results
            if (results.size() != 0)
            {
                right_type out = results[0];
                for (size_type i = 1; i < results.size(); i++)
                {
                    if (out.first > results[i].first) { out = results[i]; }
                }
                return out;
            }
            else
            {
                return end_right();
            }
        }
        
        right_type find_right_max(size_type ra_value) const
        {
            std::vector<right_type> results;
    
            // Get results computed by all threads
            for (size_type i = 0; i < right_samples.size(); i++)
            {
                auto entry = right_samples[i].second.find(ra_value);
                if ( entry != right_samples[i].second.end())
                    results.push_back(entry->second);
            }
    
            // Combine results
            if (results.size() != 0)
            {
                right_type out = results[0];
                for (size_type i = 1; i < results.size(); i++)
                {
                    if (out.first < results[i].first) { out = results[i]; }
                }
                return out;
            }
            else
            {
                return end_right();
            }
        }
        
        const left_type& find_left(size_type ra_value) const
        {
            for (size_type i = 0; i < left_samples.size(); i++)
            {
                auto entry = left_samples[i].find(ra_value);
                if ( entry != left_samples[i].end())
                    return entry->second;
            }
            return end_left();
        }
    
        ////////////////////////////////////////////////////////////
        // Generic methods
        
        std::pair<size_type, size_type> update_left(const RIndex<Derived, BWT>& left, const RIndex<Derived, BWT>& right, size_type ra_i, size_type ra_j, std::pair<size_type, size_type> prev_samples, size_type i, size_type j, size_type thread)
        {
            return this->self().update_left(left, right, ra_i, ra_j, prev_samples, i, j, thread);
        }
    
        void update_right_min(const RIndex<Derived, BWT>& left, const RIndex<Derived, BWT>& right, size_type ra_j, size_type j, size_type sa_value, size_type thread)
        {
            return this->self().update_right_min(left, right, ra_j, j, sa_value, thread);
        }
    
        void update_right_max(const RIndex<Derived, BWT>& left, const RIndex<Derived, BWT>& right, size_type ra_j, size_type j, size_type sa_value, size_type thread)
        {
            return this->self().update_right_max(left, right, ra_j, j, sa_value, thread);
        }
    };

//------------------------------------------------------------------------------

    class SamplesMerger
    {
    protected:
        const RIndex<Derived, BWT>& left; const RIndex<Derived, BWT>& right;
        std::ofstream* saes;
        const SAUpdates& updates;
    
        //////////////////////////////
        size_type counter_left = 0;
        size_type counter_right = 0;
        //////////////////////////////
    
        bool LFL = false; // last from left
        size_type LRI = 0; size_type LLI = 0; // last right index and last left index
        
        SamplesMerger(const RIndex<Derived, BWT>& r, const RIndex<Derived, BWT>& l, std::ofstream* ss, const SAUpdates& u)
        : left(l), right(r), updates(u), saes(ss), LFL(false), LRI(0), LLI(0) {}
    
        ~SamplesMerger()
        {
            if (Verbosity::level >= Verbosity::FULL)
                spdlog::info("Counter Left: {}  Counter Right: {}", counter_left, counter_right);
        }
        
    public:
    
        void set_LFL(bool value) { this->LFL = value; }
        void set_LLI(size_type value) { this->LLI = value; }
        void set_LRI(size_type value) { this->LRI = value; }
    
        ////////////////////////////////////////////////////////////
        // Generic methods
        
        void operator()(size_type index, bool FL, size_type inserting_index, size_type ra_value, size_type prev_ra_value, size_type next_ra_value)
        {
            return this->self().operator()(index, FL, inserting_index, ra_value, prev_ra_value, next_ra_value);
        }
    };
    
//------------------------------------------------------------------------------

    using bwt_type = BWT;

    byte_type terminator = DATA_TERMINATOR;
    
    // Constructors and Destructors
    RIndex() = default;
    ~RIndex() = default;
    
    // Statistics
    size_type size() const { return this->bwt_.size(); }
    bool empty() const { return (this->size() == 0); }
    size_type sequences() const { return this->sequences_; }
    size_type sigma() const { return this->alphabet.sigma(); }
    size_type runs() const { return this->bwt_.number_of_runs(); }
    size_type end_marker() const { return this->data_terminator_position; }
    
    const char* name() const { return this->self().name(); }
    
    // Metadata
    Alphabet alphabet;
    
    // Data
    const bwt_type& bwt() const { return bwt_; }
    const SA_samples& samples() const {return samples_; }
    byte_type operator[](size_type i) const { return bwt_[i]; }
    size_type get_terminator_position() const { return data_terminator_position; }
    
    size_type rank(size_type i, byte_type c) const { return this->bwt_.rank(i,c); }
    size_type select(size_type i, byte_type c) const {return this->bwt_.select(i,c); }
    
    ////////////////////////////////////////////////////////////
    // Generic methods
    
    void init(std::string bwt_path, std::string saes_path)
    {
        this->self().read_bwt(bwt_path);
        
        mio::mmap_source saes(saes_path);
        read_samples((rimerge::byte_type*) saes.data(), saes.size(), samples_.samples);
        samples_.init();
        saes.unmap();
    }
    
    // check if there is a sample in position i
    sample_genre its(size_type i) const { return this->self().its(i); }

    // get full BWT range
    range_type full_range() const { return this->self().full_range; }
    
    // access column F at position i
    byte_type F_at(size_type  i) const { return this->self().F_at(i); }
    
    // backward navigation of the BWT
    size_type LF(size_type i) const { return this->self().LF(i); }
    
    // backward navigation of the BWT
    size_type LF(size_type i, byte_type c) const { return this->self().LF(i, c); }
    
    // backward navigation of the BWt, range
    range_type LF(range_type rn, byte_type c) const { return this->self().LF(rn, c); }
    
    // forward navigation of the BWT
    size_type FL(size_type i) const { return this->self().FL(i); }
    
    // forward navigation of the BWT, where for efficiency we give c=F[i] as input
    size_type FL(size_type i, byte_type c) const { return this->self().FL(i, c); }

    // Extract i-th sequence from the index
    string_type get_sequence(size_type i) const {return this->self().get_sequence(i); }

    // Merge
    static void merge(const RIndex<Derived, BWT>& left, const RIndex<Derived, BWT>& right, MergeParameters& parameters)
    {
        return Derived::merge(left, right, parameters);
    }
    
    // Check if the structure of the index is correct
    static bool check(const RIndex<Derived, BWT>& index, std::tuple<std::vector<size_type>, std::vector<size_type>, std::vector<size_type>>& errors)
    {
        return Derived::check(index, errors);
    }
    
protected:
    
    void read_samples(byte_type* sa_start, size_type sa_size, std::vector<SA_samples::sample_type>& samples)
    {
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
    
    void read_bwt(std::string& bwt_path) { return this->read_bwt(bwt_path); }
    
    /*
     * sparse RLBWT: r (log sigma + (1+epsilon) * log (n/r)) (1+o(1)) bits
     */
    
    // F column of the BWT (vector of 256 elements)
    std::array<size_type, 256> F = {0};
    
    // L column of the BWT, run-length compressed
    bwt_type bwt_;
    size_type  data_terminator_position = 0;
    
    // SA samples
    SA_samples samples_;
    
    size_type sequences_ = 0;
};

//------------------------------------------------------------------------------

}

#endif //r_index_hpp
