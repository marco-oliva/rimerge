//
//  r-index-rled.hpp
//
//  Inspired from https://github.com/nicolaprezza/r-index
//  Gagie T, Navarro G, Prezza N. Optimal-time text indexing in BWT-runs bounded space. In Proceedings of the Twenty-Ninth Annual ACM-SIAM Symposium on Discrete Algorithms, SODA 2018, New Orleans, NA, USA, January 7-10 2017.
//


#ifndef r_index_rle_hpp
#define r_index_rle_hpp

#include <rimerge/r-index.hpp>
#include <rimerge/rle_string.hpp>

namespace rimerge
{

class RIndexRLE : public RIndex<RIndexRLE, RLEString>
{

public:

//------------------------------------------------------------------------------
    
    class SAUpdatesRLE : public RIndex<RIndexRLE, RLEString>::SAUpdates
    {
    public:
    
        SAUpdatesRLE(size_type threads) :RIndex<RIndexRLE, RLEString>::SAUpdates(threads) {}
    
        std::pair<size_type, size_type> update_left(const RIndex<RIndexRLE, RLEString>& left, const RIndex<RIndexRLE, RLEString>& right, size_type ra_i, size_type ra_j, std::pair<size_type, size_type> prev_samples, size_type i, size_type j, size_type thread, RLEString::Accessor& right_accessor, RLEString::Accessor& left_accessor);
        void update_right_min(const RIndex<RIndexRLE, RLEString>& left, const RIndex<RIndexRLE, RLEString>& right, size_type ra_j, size_type j, size_type sa_value, size_type thread, RLEString::Accessor& right_accessor, RLEString::Accessor& left_accessor);
        void update_right_max(const RIndex<RIndexRLE, RLEString>& left, const RIndex<RIndexRLE, RLEString>& right, size_type ra_j, size_type j, size_type sa_value, size_type thread, RLEString::Accessor& right_accessor, RLEString::Accessor& left_accessor);
    
        void print_size()
        {
            size_type elements_left = 0, elements_right = 0;
            for (size_type i = 0; i < left_samples.size(); i++) {elements_left += left_samples[i].size(); }
            for (size_type i = 0; i < right_samples.size(); i++) {elements_right += right_samples[i].first.size(); elements_right += right_samples[i].second.size(); }
            if(Verbosity::level >= Verbosity::FULL)
            {
                spdlog::info("rimerge::SAUpdatesRLE()::SAUpdatesRLE: {} elements in left map", elements_left);
                spdlog::info("rimerge::SAUpdatesRLE::SAUpdatesRLE): {} elements in right map", elements_right);
            }
        }
        
        ~SAUpdatesRLE(){ print_size(); }
    };

//------------------------------------------------------------------------------
    
    class SamplesMergerRLE : public RIndex<RIndexRLE, RLEString>::SamplesMerger
    {
    public:
    
        SamplesMergerRLE(const RIndex<RIndexRLE, RLEString>& r, const RIndex<RIndexRLE, RLEString>& l, std::ofstream* ss, const SAUpdatesRLE& u) : RIndex<RIndexRLE, RLEString>::SamplesMerger(r, l, ss, u) {}
        
        void operator()(size_type index, RLEString::RunCache& right_cache, RLEString::RunCache& left_cache, bool FL, size_type inserting_index, size_type ra_value, size_type prev_ra_value, size_type next_ra_value);
    };

//------------------------------------------------------------------------------
    
    
    // Constructors and Destructors
    RIndexRLE() : RIndex<RIndexRLE, RLEString>() {};
    ~RIndexRLE() = default;
    
    static const char name_[];
    const char* name() const;
    
    void read_bwt(std::string& bwt_path);

    // get full BWT range
    range_type full_range() const;
    
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

    // Extract i-th sequence from the index
    string_type get_sequence(size_type i) const;
    
    // check if there is a sample in position i
    sample_genre its(size_type i) const;
    
    // Merge
    static void merge(const RIndex<RIndexRLE, RLEString>& left, const RIndex<RIndexRLE, RLEString>& right, MergeParameters& parameters);
    
    // Check if the structure of the index is correct
    static bool check(const RIndex<RIndexRLE, RLEString>& index, std::tuple<std::vector<size_type>, std::vector<size_type>, std::vector<size_type>>& errors);
    
    // Check samples value
    static std::size_t check_sa_values(const RIndex<RIndexRLE, RLEString>& index);
  
};


}


#endif //r_index_rle_hpp
