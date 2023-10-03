//
//  sa-samples.hpp
//
// Copyright (c) Boucher Lab. All rights reserved.
// Licensed under the GNU license. See LICENSE file in the repository root for full license information.

#ifndef sa_samples_hpp
#define sa_samples_hpp

#include <rimerge/utils.hpp>
#include <rimerge/bitvector.hpp>

namespace rimerge
{

//------------------------------------------------------------------------------

struct SA_samples
{

private:

    // Rank and select to get the right sample and avoi log(n) complexity!!
    bitvector markers_bitvector;
    bool initialized = false;
    
public:
    
    using sample_type = std::pair<size_type, size_type>;
    
    std::vector<sample_type> samples;
    std::vector<size_type> compressed_samples;
    
    SA_samples() = default;
    
    // Empty the samples vector and build the bitvector and the compressed_samples
    void init();
    
    // Accessor
    size_type operator[](const size_type i) const;
    std::vector<sample_type> get_samples() const;
    
    // Output
    constexpr static size_type SAMPLE_BYTES = 5;
    static void write(sample_type& s, std::ofstream* out);
    
    static sample_type invalid_sample() {return {invalid_value(), 0}; }
    
    ~SA_samples()
    {
        if(Verbosity::level >= Verbosity::EXTENDED)
        {
            spdlog::info("rimerge::SA_samples::~SA_samples size {}GB", inGigabytes(8 * samples.size()));
        }
    };
};

//------------------------------------------------------------------------------

}

#endif //sa_samples_hpp
