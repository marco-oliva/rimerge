//
//  concat_fa.cpp
//
// Copyright (c) Boucher Lab. All rights reserved.
// Licensed under the GNU license. See LICENSE file in the repository root for full license information.


#include <string_view>
#include <fstream>

#include <mio/mio.hpp>
#include <spdlog/spdlog.h>
#include <CLI/CLI.hpp>
#include <rimerge/utils.hpp>
#include <rimerge/version.hpp>
#include <kseq.h>
#include <zlib.h>

KSEQ_INIT(gzFile, gzread)

typedef rimerge::byte_type byte_type;


int main(int argc, char **argv)
{
    std::string i_file = "";
    std::string o_file = "";
    byte_type terminator = 36;
    bool not_insert_endmarker = false;
    
    CLI::App app{"cfa"};
    app.add_option("-i,--input", i_file, "File name")->required()->check(CLI::ExistingFile)->configurable();
    app.add_option("-o,--output", o_file, "File name")->required()->configurable();
    app.add_option("-t,--terminator", terminator, "String terminator")->check(CLI::Range(0,255))->configurable();
    app.add_flag("--no-terminator", not_insert_endmarker, "Not place the end-marker between sequences")->configurable();
    app.add_flag_callback("-v,--version",rimerge::Version::print,"Version")->configurable();
    app.set_config("--configure");
    app.allow_windows_style_options();
    
    
    CLI11_PARSE(app, argc, argv);
    
    spdlog::info("Input file: {}", i_file);
    spdlog::info("Output file: {}", o_file);
    spdlog::info("Terminator: {} [{}]", terminator, char(terminator));

//------------------------------------------------------------------------------
    
    // output file
    std::fstream out(o_file, std::ios::out);
    if (!out.is_open())
    {
        spdlog::error("Failed to open {}", o_file);
        exit(EXIT_FAILURE);
    }
    
    // metadata
    std::fstream out_meta(o_file + ".meta", std::ios::out);
    if (!out_meta.is_open())
    {
        spdlog::error("Failed to open {}", o_file + ".meta");
        exit(EXIT_FAILURE);
    }
    
    // parse
    gzFile fp; kseq_t *seq; int l; std::size_t sequences = 0, tot_lenght = 0;
    fp = gzopen(i_file.c_str(), "r");
    seq = kseq_init(fp);
    while ((l = kseq_read(seq)) >= 0)
    {
        if (seq->seq.l > 0)
        {
            out.write(seq->seq.s, seq->seq.l);
            if (not not_insert_endmarker) { out.put(terminator); }
            sequences++; tot_lenght += seq->seq.l;
        }
    }
    kseq_destroy(seq); gzclose(fp);
    out.close();
    
    out_meta.write((char*) (&sequences), sizeof(sequences));
    out_meta.write((char*) (&tot_lenght), sizeof(sequences));
    out_meta.close();
    
    spdlog::info("End parsing");
}