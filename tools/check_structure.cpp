//
//  check_structure.cpp
//
//  Copyright 2020 Marco Oliva. All rights reserved.
//

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
    std::string o_prefix = "";
    
    app.add_option("-i,--input", i_prefix, "Input index prefix")->required();
    app.add_option("-o,--output", o_prefix, "Output files prefix")->required();
    app.add_flag_callback("-v,--version",rimerge::Version::print,"Version");
    
    CLI11_PARSE(app, argc, argv);
    
    spdlog::info("Input: {}", i_prefix);
    spdlog::info("Output: {}", o_prefix);
    
    rimerge::RIndex<rimerge::RIndexRLE, rimerge::RLEString> I;
    
    // Read index A bwt and sufix samples, open the file as a mmap
    spdlog::info("Reading index disk");
    read_from_prefix(i_prefix, I);
    spdlog::info("Reading I completed\nA\n\tsize: {}\n\tsequences: {}\n\tr: {}", I.size(), I.sequences(), I.runs());
    
    spdlog::info("Start checking the structure");
    
    rimerge::Verbosity::set(rimerge::Verbosity::BASIC);
    std::tuple<std::vector<rimerge::size_type>, std::vector<rimerge::size_type>, std::vector<rimerge::size_type>> errors;
    bool passed = rimerge::RIndex<rimerge::RIndexRLE, rimerge::RLEString>::check(I, errors);
    if (!passed)
    {
        spdlog::error("Missing Samples: {}", std::get<0>(errors).size());
        spdlog::error("Unnecessary Samples: {}", std::get<1>(errors).size());
        spdlog::error("Invalid Values: {}", std::get<2>(errors).size());
        spdlog::error("Total samples: {}", I.samples().samples.size());
        
        if (std::get<0>(errors).size())
        {
            std::ofstream m_f(o_prefix + ".m", std::ios::binary);
            m_f.write((char*) std::get<0>(errors).data(), std::get<0>(errors).size() * sizeof(rimerge::size_type));
        }
        if (std::get<1>(errors).size())
        {
            std::ofstream u_f(o_prefix + ".u", std::ios::binary);
            u_f.write((char*) std::get<1>(errors).data(), std::get<1>(errors).size() * sizeof(rimerge::size_type));
        }
        if (std::get<2>(errors).size())
        {
            std::ofstream i_f(o_prefix + ".i", std::ios::binary);
            i_f.write((char*) std::get<2>(errors).data(), std::get<2>(errors).size() * sizeof(rimerge::size_type));
        }
    }
    else
    {
        spdlog::info("No errors in index structure");
    }
    
    spdlog::info("SA[0]: {} Size: {}", I.samples()[0], I.size());
}