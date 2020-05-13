//
//  unit_tests.cpp
//
//  Copyright 2020 Marco Oliva. All rights reserved.
//

#include <gtest/gtest.h>

#include <rimerge/rle_string.hpp>
#include <rimerge/rindex.hpp>


namespace
{

std::string testfiles_dir = "/Users/marco/Repositories/papers/r-merge/tests/files";

//------------------------------------------------------------------------------

TEST(Bitvector, BitvectorEmptyConstructor)
{
    rimerge::bitvector bitvector;
    
    EXPECT_EQ(bitvector.size(), 0);
}

TEST(Bitvector, Rank)
{
    std::vector<bool> set = {0, 0, 1, 0, 0, 1, 1, 1, 0, 0, 1, 1};
    
    rimerge::bitvector bitvector(set);
    
    EXPECT_EQ(bitvector.size(), 12);
    
    EXPECT_EQ(0, bitvector.rank(1));
    EXPECT_EQ(0, bitvector.rank(2));
    EXPECT_EQ(1, bitvector.rank(3));
    EXPECT_EQ(1, bitvector.rank(5));
    EXPECT_EQ(3, bitvector.rank(7));
    EXPECT_EQ(5, bitvector.rank(11));
    EXPECT_EQ(6, bitvector.rank(12));
}

TEST(Bitvector, Select)
{
    std::vector<bool> set = {0, 0, 1, 0, 0, 1, 1, 1, 0, 0, 1, 1};
    
    rimerge::bitvector bitvector(set);
    
    EXPECT_EQ(bitvector.size(), 12);
    
    EXPECT_EQ(2, bitvector.select(0));
    EXPECT_EQ(7, bitvector.select(3));
    EXPECT_EQ(11, bitvector.select(5));
}

TEST(RLEString, EmptyConstructor)
{
    rimerge::rle_string rle_test;
    EXPECT_EQ(0, rle_test.size());
}

TEST(RLEString, ShortStringConstructor)
{
    rimerge::string_type short_string = {'A', 'A', 'A', 'A', 'A', 'B', 'C', 'B'};
    rimerge::rle_string rle_test(short_string);
    
    EXPECT_EQ(8, rle_test.size());
    EXPECT_EQ(short_string, rle_test.to_string());
    
    EXPECT_EQ('A', rle_test[2]);
    EXPECT_EQ('B', rle_test[7]);
    
    EXPECT_EQ(3, rle_test.rank(3, 'A'));
    EXPECT_EQ(0, rle_test.rank(3, 'B'));
    EXPECT_EQ(1, rle_test.rank(6, 'B'));
    
    EXPECT_EQ(7, rle_test.select(1, 'B'));
    EXPECT_EQ(0, rle_test.select(0, 'A'));
    EXPECT_EQ(5, rle_test.select(0, 'B'));
}

TEST(RLEString, CharInARange)
{
    rimerge::string_type short_string = {'A', 'A', 'A', 'A', 'A', 'B', 'C', 'C'};
    rimerge::rle_string rle_test(short_string);
    
    rimerge::size_type n = rle_test.rank(rle_test.size(), 'B') - rle_test.rank(4, 'B');
    EXPECT_EQ(1, n);
}

TEST(RLEString, FromFile)
{
    std::string path = testfiles_dir + "/rle_string_plain.txt";
    mio::mmap_source input(path);
    
    rimerge::rle_string from_file(input);
    
    rimerge::string_type short_string = {'A', 'A', 'A', 'A', 'A', 'B', 'B', 'B'};
    
    EXPECT_EQ(8, from_file.size());
    EXPECT_EQ(short_string, from_file.to_string());
    
    EXPECT_EQ('A', from_file[2]);
    EXPECT_EQ('B', from_file[7]);
    
    EXPECT_EQ(3, from_file.rank(3, 'A'));
    EXPECT_EQ(0, from_file.rank(3, 'B'));
    
    EXPECT_EQ(7, from_file.select(2, 'B'));
}

//------------------------------------------------------------------------------

TEST(Alphabet, ThreeCharacters)
{
    rimerge::Alphabet alphabet;
    alphabet.update('A');
    alphabet.update('C');
    alphabet.update('G');
    
    alphabet.init();
    
    EXPECT_EQ(3, alphabet.sigma());
    
    EXPECT_EQ('A', alphabet.previous('C'));
    EXPECT_EQ('C', alphabet.previous('G'));
    EXPECT_EQ(rimerge::STRING_TERMINATOR, alphabet.previous('A'));
    
    EXPECT_EQ('C', alphabet.following('A'));
    EXPECT_EQ('G', alphabet.following('C'));
    EXPECT_EQ(rimerge::STRING_TERMINATOR, alphabet.following('G'));
}

TEST(Rindex, EmptyConstructor)
{
    rimerge::rindex empty_index;
    
    EXPECT_EQ(0, empty_index.size());
    EXPECT_EQ(0, empty_index.sigma());
    EXPECT_EQ(0, empty_index.runs());
}

TEST(Rindex, FromFilesShort)
{
    mio::mmap_source bwt(testfiles_dir + "/test.bwt");
    mio::mmap_source ssa(testfiles_dir + "/test.ssa");
    mio::mmap_source esa(testfiles_dir + "/test.esa");
    
    rimerge::rindex index((rimerge::byte_type*) bwt.data(), bwt.size(),
                          (rimerge::byte_type*) ssa.data(), ssa.size(),
                          (rimerge::byte_type*) esa.data(), esa.size());
    
    EXPECT_EQ(27, index.size());
}

TEST(Rindex, FromFilesLong)
{
    mio::mmap_source bwt(testfiles_dir + "/long.bwt");
    mio::mmap_source ssa(testfiles_dir + "/long.ssa");
    mio::mmap_source esa(testfiles_dir + "/long.esa");
    
    rimerge::rindex index((rimerge::byte_type*) bwt.data(), bwt.size(),
                          (rimerge::byte_type*) ssa.data(), ssa.size(),
                          (rimerge::byte_type*) esa.data(), esa.size());
    
    EXPECT_EQ(755, index.size());
    EXPECT_EQ(5, index.sequences());
}

TEST(Rindex, RindexGeneratorConstructor)
{
    rimerge::testing::rindex_generator generator;
    rimerge::rindex index = generator.get();
    
    EXPECT_EQ(27, index.bwt().size());
    EXPECT_EQ(3, index.sequences());
}

TEST(Rindex, LFMapping)
{
    rimerge::testing::rindex_generator generator;
    rimerge::rindex index = generator.get();
    
    EXPECT_EQ(3, index.LF(0));
    EXPECT_EQ(12, index.LF(26));
    EXPECT_EQ(13, index.LF(7));
}

TEST(Rindex, LFMappingRange)
{
    rimerge::testing::rindex_generator generator;
    rimerge::rindex index = generator.get();
    
    rimerge::range_type range = index.full_range();
    rimerge::range_type new_range = index.LF(range, 'T');
    
    EXPECT_EQ(rimerge::range_type(19, 26), new_range); // add more
}

TEST(Rindex, FLMapping)
{
    rimerge::testing::rindex_generator generator;
    rimerge::rindex index = generator.get();
    
    EXPECT_EQ(17, index.FL(0));
    EXPECT_EQ(24, index.FL(26));
    EXPECT_EQ(7, index.FL(13));
    
    EXPECT_EQ(17, index.FL(0, rimerge::DATA_TERMINATOR));
    EXPECT_EQ(24, index.FL(26, 'T'));
    EXPECT_EQ(7, index.FL(13, 'C'));
}

//------------------------------------------------------------------------------

} // namespace