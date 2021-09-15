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
#include <rimerge/r-index-rle.hpp>


int main(int argc, char **argv)
{
    rimerge::Verbosity::set(rimerge::Verbosity::BASIC);
    
    CLI::App app("rle");
    
    std::string  e_path = "";
    std::string  d_path = "";
    std::string ex_path = "";

    app.add_option("-i,--input", e_path, "Input bwt")->configurable();
    app.add_option("-d, --decode", d_path, "Input rle")->configurable();
    app.add_flag_callback("-v,--version", rimerge::Version::print,"Version");
    app.set_config("--configure");
    app.allow_windows_style_options();

    CLI11_PARSE(app, argc, argv);

    if (e_path != "")
    {
        spdlog::info("Input: {}", e_path);
        spdlog::info("Encoding RLE");
        mio::mmap_source bwt(e_path);
        rimerge::byte_type* start = (rimerge::byte_type*) bwt.data();
        rimerge::byte_type* end   = start + bwt.size();

        rimerge::RLEString::RLEncoder encoder(e_path + ".rle");
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
    else if (d_path != "")
    {
        spdlog::info("Input: {}", d_path);
        spdlog::info("Decoding RLE");
        rimerge::RLEString::RLEDecoder decoder(d_path);

        std::ofstream out_file(d_path + ".ext");
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
        spdlog::error("Specify --input, --decode or --extract-sequences");
    }
}