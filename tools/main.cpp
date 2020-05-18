//
//  main.cpp
//
//  Copyright 2020 Marco Oliva. All rights reserved.
//

#include <CLI/CLI.hpp>
#include <spdlog/spdlog.h>
#include <mio/mio.hpp>

#include <rimerge/utils.hpp>
#include <rimerge/rindex.hpp>
#include <rimerge/rle_string.hpp>
#include <rimerge/version.hpp>

static std::string    bwt_ext = ".bwt";
static std::string starts_ext = ".ssa";
static std::string   ends_ext = ".esa";
static std::string    nsa_ext = ".nsa";

int main(int argc, char **argv)
{
    rimerge::Verbosity::set(rimerge::Verbosity::EXTENDED);
    
    CLI::App app{"rimerge"};
    
    std::string a_prefix = "";
    std::string b_prefix = "";
    std::string o_prefix = "";
    
    app.add_option("-a,--aprefix", a_prefix, "Prefix of A")->required();
    app.add_option("-b,--bprefix", b_prefix, "Prefix of B")->required();
    app.add_option("-o, --output", o_prefix, "Prefix of the output")->required();
    app.add_flag_callback("-v,--version",rimerge::Version::print,"Version");
    
    CLI11_PARSE(app, argc, argv);
    
    spdlog::info("A: {}", a_prefix);
    spdlog::info("B: {}", b_prefix);
    
    rimerge::rindex A;
    rimerge::rindex B;
    
    // Read index A bwt and sufix samples, open the file as a mmap
    spdlog::info("Reading A from disk");
    mio::mmap_source bwt_a(a_prefix + bwt_ext);
    mio::mmap_source ssa_a(a_prefix + starts_ext);
    mio::mmap_source esa_a(a_prefix + ends_ext);
    mio::mmap_source nsa_a(a_prefix + nsa_ext);
    A.init((rimerge::byte_type*) bwt_a.data(), bwt_a.size(),
           (rimerge::byte_type*) ssa_a.data(), ssa_a.size(),
           (rimerge::byte_type*) esa_a.data(), esa_a.size(),
           (rimerge::byte_type*) nsa_a.data(), nsa_a.size());
    bwt_a.unmap(); ssa_a.unmap(); esa_a.unmap(); nsa_a.unmap();
    spdlog::info("Reading A completed\nA\n\tsize: {}\n\tsequences: {}\n\tr: {}", A.size(), A.sequences(), A.runs());
    
    spdlog::info("Reading B from disk");
    mio::mmap_source bwt_b(b_prefix + bwt_ext);
    mio::mmap_source ssa_b(b_prefix + starts_ext);
    mio::mmap_source esa_b(b_prefix + ends_ext);
    mio::mmap_source nsa_b(b_prefix + nsa_ext);
    B.init((rimerge::byte_type*) bwt_b.data(), bwt_b.size(),
           (rimerge::byte_type*) ssa_b.data(), ssa_b.size(),
           (rimerge::byte_type*) esa_b.data(), esa_b.size(),
           (rimerge::byte_type*) nsa_b.data(), nsa_b.size());
    bwt_b.unmap(); ssa_b.unmap(); esa_b.unmap(); nsa_b.unmap();
    spdlog::info("Reading B completed\nB\n\tsize: {}\n\tsequences: {}\n\tr: {}", B.size(), B.sequences(), B.runs());
    
    rimerge::MergeParameters merge_parameters;
    merge_parameters.out_prefix = o_prefix;
    rimerge::rindex::merge(A, B, merge_parameters);
    
    spdlog::info("Done");
}