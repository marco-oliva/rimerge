//
//  concat_fa.cpp
//
//  Copyright 2020 Marco Oliva. All rights reserved.
//

#include <string_view>
#include <fstream>

#include <mio/mio.hpp>
#include <spdlog/spdlog.h>
#include <CLI/CLI.hpp>

typedef std::uint8_t  byte_type;

int main(int argc, char **argv)
{
    std::string i_file = "";
    std::string o_file = "";
    byte_type terminator = 36;
    
    CLI::App app{"cfa"};
    app.add_option("-i,--input", i_file, "File name")->required()->check(CLI::ExistingFile);
    app.add_option("-o,--output", o_file, "File name")->required();
    app.add_option("-t,--terminator", terminator, "String terminator")->check(CLI::Range(0,255));
    
    CLI11_PARSE(app, argc, argv);
    
    spdlog::info("Input file: {}", i_file);
    spdlog::info("Output file: {}", o_file);
    spdlog::info("Terminator: {} [{}]", terminator, char(terminator));

//------------------------------------------------------------------------------
    
    // input file
    std::error_code error;
    mio::mmap_source fa_map;
    fa_map.map(i_file, error);
    if (error)
    {
        const auto& errmsg = error.message();
        spdlog::error(errmsg);
        exit(1);
    }
    
    // output file
    std::fstream out(o_file, std::ios::out);
    if (!out.is_open())
    {
        spdlog::error("Failed to open {}", o_file);
        exit(1);
    }
    
    
    const char* end = fa_map.data() + fa_map.size();
    
    std::string_view read;
    
    // start parsing
    spdlog::info("Start parsing");
    const char *l_1, *l_2, *e_l;
    l_1 = fa_map.data();
    l_2 = static_cast<const char*>(memchr(l_1, '\n', end - l_1)) + 1;
    e_l = static_cast<const char*>(memchr(l_2, '\n', end - l_2));
    read = std::string_view(l_2, e_l - l_2 - 1);
    
    out.write(read.data(), read.size());
    out.put(terminator);
    
    // keep parsing
    while (e_l + 1 < end)
    {
        // Next read
        l_1 = e_l + 1;
        l_2 = static_cast<const char*>(memchr(l_1, '\n', end - l_1)) + 1;
        e_l = static_cast<const char*>(memchr(l_2, '\n', end - l_2));
    
        read = std::string_view(l_2, e_l - l_2 - 1);
        if (read.size() > 0)
        {
            out.write(read.data(), read.size());
            out.put(terminator);
        }
        else
        {
            spdlog::error("Empty sequence found, skipping");
        }
        
    }
    
    spdlog::info("End parsing");
}

