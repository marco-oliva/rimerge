//
//  benchmarks.cpp
//
//  Copyright 2020 Marco Oliva. All rights reserved.
//

#include <benchmark/benchmark.h>

#include <rimerge/rle_string.hpp>
#include <rimerge/rle_string.hpp>


//------------------------------------------------------------------------------

using namespace rimerge;

std::string test_bwt = "gcttcattcaatgtggcaatgatgatttaaggggttagcacactctataaaggtcaagcctcaggccgtcaggctaaaatcggtagttgtaatttttcaacgtccaactgggatgtcggtcctgaaacgtcagagtcgagtctagtcgcaggttggcatgcgatcacgcacacgtagtacaacccgcatagaagacccgtgttttattcgcctcattatgtattttttgtacccacgttctcagattagtggtgcatccaaaaatcccaacatgtcactgtgctgacagtatctcccaaacgttagggaggcacgcctctggtcccaaaaacacaggtgtttcctccgtacgatccttgccgcaggaactagtccaagtcactgggtccgggcgcagatggagcgtcagaagtctaatagctgtcggctacggcttggttactatcagaatgtccacgcctaaaacccttatgtactgacattttagatagctagacaaccggtcctccccggtactttgataggcgctatctacatggcgtagctatgccagactcccacttcttctaaacccctccctcagggtatggctatctgcttttttaaagagattcgccacaaaaaactccagtaaccactactcccagagcctgtagagtttcagatgccctaactcagctactggcactcggttagctatcaatttcactaaatgtccagttaaaactgttccacatttcgatattgtcgaaagggcgagccccagatacagtcacgagtcgaaaactttaggctccaccccacccgactaactcaagtcgttgacgtgcataattttagcgtacgcaagttccggttgtgtggttccctgtatgagtcgtagattgtctccgcgactg$tcaagaccttcacccgttgtgtaatctcgtacaaccatcactgactatgacatcacactcgggcctgaaggcaatgcaagcctggctgcgtaaaggcggaa";

//------------------------------------------------------------------------------

template <class BWT>
void full_LF(BWT& bwt, std::array<size_type, 256>& F)
{
    // Complete LF
    size_type i = 0; byte_type c = bwt[i];
    while (c != '$')
    {
        i = F[c] + bwt.rank(i,c);
        c = bwt[i];
    }
}
//------------------------------------------------------------------------------
// Random Access

static void BMRLE_RA(benchmark::State& state)
{
    // Setup
    Verbosity::set(Verbosity::SILENT);
    rimerge::string_type bwt;
    for (auto& c : test_bwt) { bwt.push_back(c); }
    rimerge::RLEString rle_string(bwt);
    
    for (auto _ : state)
    {
        for (size_type i = 0; i < rle_string.size(); i++)
        {
            // Multiple times at same position
            for (size_type j = 0; j < 10; j++) { byte_type c = rle_string[i]; }
        }
    }
}

static void BMRLE_RA_ACCESSOR(benchmark::State& state)
{
    // Setup
    Verbosity::set(Verbosity::SILENT);
    rimerge::string_type bwt;
    for (auto& c : test_bwt) { bwt.push_back(c); }
    rimerge::RLEString rle_string(bwt);
    rimerge::RLEString::Accessor rle_string_accessor(rle_string);
    
    for (auto _ : state)
    {
        for (size_type i = 0; i < rle_string.size(); i++)
        {
            // Multiple times at same position
            for (size_type j = 0; j < 10; j++) { byte_type c = rle_string_accessor[i]; }
        }
    }
}

// Register the function as a benchmark
BENCHMARK(BMRLE_RA);
BENCHMARK(BMRLE_RA_ACCESSOR);

// Run the benchmark
BENCHMARK_MAIN();