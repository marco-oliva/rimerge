//
//  rle.cpp
//
//  Copyright 2020 Marco Oliva. All rights reserved.
//

#include <CLI/CLI.hpp>
#include <spdlog/spdlog.h>
#include <mio/mio.hpp>

#include <rimerge/utils.hpp>
#include <rimerge/rle_string.hpp>
#include <rimerge/version.hpp>


int main(int argc, char **argv)
{
    rimerge::Verbosity::set(rimerge::Verbosity::BASIC);
    
    CLI::App app("rle");
    
    std::string e_prefix = "";
    std::string d_prefix = "";

    app.add_option("-i,--input", e_prefix, "Input bwt")->configurable();
    app.add_option("-d, --decode", d_prefix, "Input rle")->configurable();
    app.add_flag_callback("-v,--version", rimerge::Version::print,"Version");
    app.set_config("--configure");
    app.allow_windows_style_options();

    CLI11_PARSE(app, argc, argv);

    if (e_prefix != "")
    {
        spdlog::info("Input: {}", e_prefix);
        spdlog::info("Encoding RLE");
        mio::mmap_source bwt(e_prefix + ".bwt");
        rimerge::byte_type* start = (rimerge::byte_type*) bwt.data();
        rimerge::byte_type* end   = start + bwt.size();

        rimerge::RLEString::RLEncoder encoder(e_prefix + ".rle");
        rimerge::byte_type* it = start; rimerge::byte_type to_ins;
        while (it != end)
        {
            if (*it == static_cast<rimerge::byte_type>(0)) { to_ins = rimerge::STRING_TERMINATOR; }
            else { to_ins = *it; }
            encoder(to_ins);
            it++;
        }

        return 0;
    }
    else if (d_prefix != "")
    {
    
        std::string d_prefix = argv[1];
        spdlog::info("Input: {}", d_prefix);
        spdlog::info("Decoding RLE");
        rimerge::RLEString::RLEDecoder decoder(d_prefix + ".rle");

        std::ofstream out_file(d_prefix + ".bwt");
        while (not decoder.end())
        {
            rimerge::RunType run = decoder.next();
            rimerge::byte_type c = rimerge::RunTraits::charachter(run);
            rimerge::size_type len = rimerge::RunTraits::length(run);

            if (len != 0)
            {
                std::fill_n(std::ostream_iterator<char>(out_file), len, c);
            }
        }
        out_file.close();

        return 0;
    }
    else
    {
        spdlog::error("Specify either --input or --decode");
    }
}