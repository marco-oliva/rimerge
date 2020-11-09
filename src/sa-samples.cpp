//
//  sa-samples.cpp
//
//  Copyright 2020 Marco Oliva. All rights reserved.
//


#include <rimerge/sa-samples.hpp>

namespace rimerge
{

//------------------------------------------------------------------------------

constexpr size_type SA_samples::SAMPLE_BYTES;

void
SA_samples::init()
{
    if (initialized) { return; }
    this->initialized = true;
    
    size_type num_of_markers = this->samples.size();
    size_type tot_length = std::get<0>(*(std::max_element(this->samples.begin(), this->samples.end())));
    
    sdsl::sd_vector_builder markers_bv_builder(tot_length, num_of_markers);
    size_type runs_bv_it = 0;
    
    for (auto& entry : this->samples)
    {
        markers_bv_builder.set(entry.first);
        this->compressed_samples.push_back(entry.second);
    }
    this->markers_bitvector = bitvector(markers_bv_builder);
    this->samples.empty();
}

size_type
SA_samples::operator[](const rimerge::size_type i) const
{
    if (not initialized) { spdlog::error("Initialize SA_samples before access"); exit(EXIT_FAILURE);}
    if (i > this->markers_bitvector.size()) {spdlog::error("SA_samples::operator[]: invalid access"); exit(EXIT_FAILURE); }
    if (not markers_bitvector[i])
    {
        if (Verbosity::level >= Verbosity::EXTENDED) {spdlog::warn("rimerge::SA_samples::operator[](): Asked for not sampled value: {}", i);}
        return invalid_value();
    }
    
    size_type pos = markers_bitvector.rank(i);
    return this->compressed_samples[pos];
}

void
SA_samples::write(SA_samples::sample_type& s, std::ofstream* out)
{
    out->write(reinterpret_cast<const char *>(&s.first), SA_samples::SAMPLE_BYTES);
    out->write(reinterpret_cast<const char *>(&s.second), SA_samples::SAMPLE_BYTES);
}

//------------------------------------------------------------------------------

}