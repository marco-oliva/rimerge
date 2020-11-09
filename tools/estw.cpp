//
//  estw.cpp
//
//  Copyright 2020 Marco Oliva. All rights reserved.
//

#include <string>
#include <vector>

#include <spdlog/spdlog.h>
#include <CLI/CLI.hpp>
#include <mio/mio.hpp>
#include <rimerge/r-index.hpp>

typedef rimerge::size_type size_type;
constexpr size_type SAMPLE_BYTES = rimerge::SA_samples::SAMPLE_BYTES;

std::pair<size_type, size_type>
get_index_and_value(char* pointer)
{
    size_type key = 0, value = 0;
    std::memcpy(&key,   pointer, SAMPLE_BYTES);
    std::memcpy(&value, pointer + SAMPLE_BYTES, SAMPLE_BYTES);
    return std::make_pair(key, value);
}

void
write_pair(std::pair<size_type, size_type>& pair, std::ofstream& out_file)
{
    out_file.write(reinterpret_cast<const char*>(&(pair.first)), SAMPLE_BYTES);
    out_file.write(reinterpret_cast<const char*>(&(pair.second)), SAMPLE_BYTES);
}

void
merge_samples(mio::mmap_source& file_a, mio::mmap_source& file_b, std::ofstream& out_file)
{
    // Iterators
    char *iterator_a, *iterator_b, *end_a, *end_b;
    iterator_a = const_cast<char*>(file_a.data()); end_a = const_cast<char*>(file_a.data()) + file_a.size();
    iterator_b = const_cast<char*>(file_b.data()); end_b = const_cast<char*>(file_b.data()) + file_b.size();

    // Merge the two files, no duplicates
    std::pair<size_type, size_type> last_wrote = {0,0};
    while (iterator_a != end_a and iterator_b != end_b)
    {
        std::pair<size_type, size_type> to_write;
        if (get_index_and_value(iterator_a) < get_index_and_value(iterator_b))
        {
            to_write = get_index_and_value(iterator_a);
            iterator_a += 2 * SAMPLE_BYTES;
        }
        else
        {
            to_write = get_index_and_value(iterator_b);
            iterator_b += 2 * SAMPLE_BYTES;
        }
        
        if (to_write != last_wrote) { write_pair(to_write, out_file); last_wrote = to_write; }
    }
    while (iterator_a != end_a)
    {
        std::pair<size_type, size_type> to_write = get_index_and_value(iterator_a);
        iterator_a += 2 * SAMPLE_BYTES;
        if (to_write != last_wrote) { write_pair(to_write, out_file); last_wrote = to_write; }
    }
    while (iterator_b != end_b)
    {
        std::pair<size_type, size_type> to_write = get_index_and_value(iterator_b);
        iterator_b += 2 * SAMPLE_BYTES;
        if (to_write != last_wrote) { write_pair(to_write, out_file); last_wrote = to_write; }
    }
    out_file.close();
}

int main(int argc, char **argv)
{
    std::string i_prefix = "";
    
    CLI::App app("estw");
    app.add_option("-i,--input", i_prefix, "Input Prefix")->required();
    CLI11_PARSE(app, argc, argv);
    
    spdlog::info("Input prefix: {}", i_prefix);
    spdlog::info("Output file: {}", i_prefix + ".saes");
    //------------------------------------------------------------------------------
    
    // merge samples, first pass
    mio::mmap_source ssa(i_prefix + ".ssa");
    mio::mmap_source esa(i_prefix + ".esa");
    
    std::string tmp_file_name = rimerge::TempFile::getName("estw");
    std::ofstream tmp_file(tmp_file_name, std::ios::binary);
    spdlog::info("Merging .ssa and .esa into {}", tmp_file_name);
    merge_samples(ssa, esa, tmp_file);
    tmp_file.close();
    
    // merge samples, second pass
    mio::mmap_source nsa(i_prefix + ".nsa");
    mio::mmap_source tmp_mmap(tmp_file_name);
    
    std::ofstream out_samples(i_prefix + ".saes", std::ios::binary);
    spdlog::info("Merging {} and .nsa into {}", tmp_file_name, i_prefix + ".saes");
    merge_samples(tmp_mmap, nsa, out_samples);
    out_samples.close();
    
    spdlog::info("Done");
}


