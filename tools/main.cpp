//
//  main.cpp
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
    
    std::string a_prefix = "";
    std::string b_prefix = "";
    std::string o_prefix = "";
    rimerge::size_type merge_jobs = rimerge::MergeParameters::MERGE_JOBS;
    rimerge::size_type pos_size = rimerge::MergeParameters::POS_BUFFER_SIZE;
    rimerge::size_type thr_size = rimerge::MergeParameters::THREAD_BUFFER_SIZE;
    
    bool check_flag = false;
    
    app.add_option("-a,--aprefix", a_prefix, "Prefix of A")->required()->configurable();
    app.add_option("-b,--bprefix", b_prefix, "Prefix of B")->required()->configurable();
    app.add_option("-o, --output", o_prefix, "Prefix of the output")->required()->configurable();
    app.add_option("--pos-buffer", pos_size, "Positions buffer size")->configurable();
    app.add_option("--thread-buffer", thr_size, "Thread buffer size")->configurable();
    app.add_option("-j, --merge-jobs", merge_jobs, "Number of merge jobs")->required()->configurable();
    app.add_flag("-c, --check", check_flag, "Check index astructure after merging")->configurable();
    app.add_flag_callback("-v,--version",rimerge::Version::print,"Version");
    app.set_config("--configure");
    app.allow_windows_style_options();
    
    CLI11_PARSE(app, argc, argv);
    
    spdlog::info("A: {}", a_prefix);
    spdlog::info("B: {}", b_prefix);
    
    rimerge::RIndex<RIndexType, BWTType> A;
    rimerge::RIndex<RIndexType, BWTType> B;
    
    double reading = rimerge::readTimer();
    
    // Read index A bwt and sufix samples, open the file as a mmap
    spdlog::info("Reading A from disk");
    read_from_prefix(a_prefix, A);
    spdlog::info("Reading A completed\nA\n\tsize: {}\n\tsequences: {}\n\tr: {}\n\tn/r: {:03.2f}", A.size(), A.sequences(), A.runs(), float(A.size())/float(A.runs()));
    
    spdlog::info("Reading B from disk");
    read_from_prefix(b_prefix, B);
    spdlog::info("Reading B completed\nB\n\tsize: {}\n\tsequences: {}\n\tr: {}\n\tn/r: {:03.2f}", B.size(), B.sequences(), B.runs(), float(B.size())/float(B.runs()));
    
    double reading_end = rimerge::readTimer();
    spdlog::info("Wall time for reading: {} seconds", reading_end - reading);
    
    double merging = rimerge::readTimer();
    
    rimerge::MergeParameters merge_parameters;
    merge_parameters.out_prefix = o_prefix;
    merge_parameters.setPosBufferSize(pos_size);
    merge_parameters.setThreadBufferSize(thr_size);
    merge_parameters.setMergeJobs(merge_jobs);
    merge_parameters.setMergeBuffers(merge_jobs);
    rimerge::RIndex<RIndexType, BWTType>::merge(A, B, merge_parameters);
    
    double merging_end = rimerge::readTimer();
    spdlog::info("Wall time for merging: {} seconds", merging_end - merging);
    
    spdlog::info("Done");
}