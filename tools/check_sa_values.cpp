//
//  check_correctness.cpp.c
//
// Copyright (c) Boucher Lab. All rights reserved.
// Licensed under the GNU license. See LICENSE file in the repository root for full license information.


#include <CLI/CLI.hpp>
#include <spdlog/spdlog.h>
#include <mio/mio.hpp>

#include <rimerge/utils.hpp>
#include <rimerge/r-index-rle.hpp>
#include <rimerge/rle_string.hpp>
#include <rimerge/version.hpp>


template<typename I, typename BWT>
void read_from_prefix(std::string& prefix, rimerge::RIndex<I, BWT>& index)
{
    index.init(prefix + "/bwt.rle", prefix + "/samples.saes");
}

int main(int argc, char **argv)
{
    rimerge::Verbosity::set(rimerge::Verbosity::BASIC);
    
    CLI::App app{"check"};
    
    std::string i_prefix = "";
    
    app.add_option("-i,--input", i_prefix, "Input index prefix")->required();
    app.add_flag_callback("-v,--version",rimerge::Version::print,"Version");
    
    CLI11_PARSE(app, argc, argv);
    
    spdlog::info("Input: {}", i_prefix);
    rimerge::RIndex<rimerge::RIndexRLE, rimerge::RLEString> I;
    
    // Read index A bwt and sufix samples, open the file as a mmap
    spdlog::info("Reading index disk");
    read_from_prefix(i_prefix, I);
    spdlog::info("Reading I completed\nA\n\tsize: {}\n\tsequences: {}\n\tr: {}", I.size(), I.sequences(), I.runs());
    
    spdlog::info("Start checking SA values");
    
    rimerge::size_type errors = 0;
    errors = rimerge::RIndex<rimerge::RIndexRLE, rimerge::RLEString>::check_sa_values(I);
    
    if (errors != 0) { spdlog::error("Errors: {}", errors); }
    else { spdlog::info("Errors: {}", errors); }
}