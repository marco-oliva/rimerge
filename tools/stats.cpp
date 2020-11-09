//
//  stats.cpp
//
//  Copyright 2020 Marco Oliva. All rights reserved.
//

#include <CLI/CLI.hpp>
#include <spdlog/spdlog.h>
#include <mio/mio.hpp>

#include <rimerge/utils.hpp>
#include <rimerge/r-index-rle.hpp>
#include <rimerge/version.hpp>


template<typename I, typename BWT>
void read_from_prefix(std::string& prefix, rimerge::RIndex<I, BWT>& index)
{
    index.init(prefix + "/bwt.rle", prefix + "/samples.saes");
}

int main(int argc, char **argv)
{
    // r-index type
    typedef rimerge::RIndexRLE RIndexType;
    typedef rimerge::RLEString BWTType;
    
    
    rimerge::Verbosity::set(rimerge::Verbosity::FULL);
    
    CLI::App app("rimerge");
    
    std::string i_prefix = "";
    
    app.add_option("-i,--input", i_prefix, "Index Prefix")->required();
    app.add_flag_callback("-v,--version",rimerge::Version::print,"Version");
    
    CLI11_PARSE(app, argc, argv);
    
    spdlog::info("Index: {}", i_prefix);
    
    rimerge::RIndex<RIndexType, BWTType> I;
    
    double reading = rimerge::readTimer();
    
    // Read index A bwt and sufix samples, open the file as a mmap
    spdlog::info("Reading A from disk");
    read_from_prefix(i_prefix, I);
    spdlog::info("Reading A completed\nA\n\tsize: {}\n\tsequences: {}\n\tr: {}\n\tn/r: {:03.2f}", I.size(), I.sequences(), I.runs(), float(I.size())/float(I.runs()));
    
    spdlog::info("Done");
}