//
//  unit_tests.cpp
//
//  Copyright 2020 Marco Oliva. All rights reserved.
//

#include <gtest/gtest.h>

#include <rimerge/utils.hpp>
#include <rimerge/rle_string.hpp>
#include <rimerge/r-index-rle.hpp>


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

rimerge::string_type random_string(rimerge::size_type length)
{
    srand((rimerge::size_type) rimerge::readTimer());

    // See: https://stackoverflow.com/questions/440133/how-do-i-create-a-random-alpha-numeric-string-in-c
    auto randchar = []() -> char
    {
        const rimerge::byte_type charset[] = "ACGTN";
        const size_t max_index = (sizeof(charset) - 1);
        return charset[ rand() % max_index ];
    };
    rimerge::string_type str(length,0);
    std::generate_n( str.begin(), length, randchar );
    str.push_back('$');
    return str;
}

TEST(RLEString, FromEncoder)
{
    std::string path = testfiles_dir + "/test.rle";
    rimerge::string_type string = {'A','A','A','G','C','A','A','T','T','T','T','T','T','A','A','G','G','A','A','A'};
    rimerge::RLEString rlet_string(string);


    rimerge::RLEString::RLEncoder encoder(path);
    for (auto e : string) { encoder(e); }
    encoder.close();

    rimerge::RLEString rle_string;
    rle_string.load(path);

    EXPECT_EQ(string.size(), rle_string.size());
    for (std::size_t i = 0; i < string.size(); i++)
    {
        EXPECT_EQ(string[i], rle_string[i]);
    }
}

TEST(RLESEncoder, RandomSequence)
{
    for (rimerge::size_type r = 0; r < 10; r++)
    {
        std::string path = testfiles_dir + "/rnd_seq.seq";
        rimerge::string_type rnds = random_string(10000);

        rimerge::RLEString::RLEncoder encoder(path);

        for (auto& e : rnds) { encoder(e); }
        encoder.close();

        rimerge::RLEString rle_string;
        rle_string.load(path);

        for (rimerge::size_type i = 0; i < rnds.size(); i++)
        {
            EXPECT_EQ(rnds[i], rle_string[i]);
        }
    }
}

TEST(RLEString, FromConcatRunInterruption)
{
    std::string path = testfiles_dir + "/test_cnct.seq";
    rimerge::RLEString::RLEncoderMerger encoders(path, 2);

    rimerge::string_type string_1 = {'A','A','A','G','C','A','A','T','T','T','T','T','T','A','A','G','G','A','A','A'};
    rimerge::string_type string_2 = {'G','A','A','G','C','A','A','T','T','T','T','T','T','A','A','G','G','A','A','A'};

    rimerge::string_type string;
    string.reserve(string_1.size() + string_2.size());
    string.insert(string.end(), string_1.begin(), string_1.end());
    string.insert(string.end(), string_2.begin(), string_2.end());

    rimerge::RLEString::RLEncoder& encoder_1 = encoders.get_encoder(0);
    for (auto e : string_1)
    { encoder_1(e); }

    rimerge::RLEString::RLEncoder& encoder_2 = encoders.get_encoder(1);
    for (auto e : string_2)
    { encoder_2(e); }

    encoders.close();

    rimerge::RLEString rle_string;
    rle_string.load(testfiles_dir + "/test_cnct.seq");

    EXPECT_EQ(string.size(), rle_string.size());
    for (std::size_t i = 0; i < string.size(); i++)
    {
        EXPECT_EQ(string[i], rle_string[i]);
    }
}

TEST(RLEString, FromConcatNoRunInterruption)
{
    std::string path = testfiles_dir + "/test_cnct.seq";
    rimerge::RLEString::RLEncoderMerger encoders(path, 2);

    rimerge::string_type string_1 = {'A','A','A','G','C','A','A','T','T','T','T','T','T','A','A','G','G','A','A','A'};
    rimerge::string_type string_2 = {'A','A','A','G','C','A','A','T','T','T','T','T','T','A','A','G','G','A','A','A'};

    rimerge::string_type string;
    string.reserve(string_1.size() + string_2.size());
    string.insert(string.end(), string_1.begin(), string_1.end());
    string.insert(string.end(), string_2.begin(), string_2.end());

    rimerge::RLEString::RLEncoder& encoder_1 = encoders.get_encoder(0);
    for (auto e : string_1)
    { encoder_1(e); }

    rimerge::RLEString::RLEncoder& encoder_2 = encoders.get_encoder(1);
    for (auto e : string_2)
    { encoder_2(e); }

    encoders.close();

    rimerge::RLEString rle_string;
    rle_string.load(testfiles_dir + "/test_cnct.seq");

    EXPECT_EQ(string.size(), rle_string.size());
    for (std::size_t i = 0; i < string.size(); i++)
    {
        EXPECT_EQ(string[i], rle_string[i]);
    }
}

TEST(RLEString, RandomConcat)
{
    for (rimerge::size_type r = 0; r < 10; r++)
    {
        std::string path = testfiles_dir + "/rnd_cnct.seq";
        rimerge::RLEString::RLEncoderMerger encoders(path, 2);

        rimerge::string_type string_1 = random_string(1000);
        rimerge::string_type string_2 = random_string(100);

        rimerge::string_type string;
        string.reserve(string_1.size() + string_2.size());
        string.insert(string.end(), string_1.begin(), string_1.end());
        string.insert(string.end(), string_2.begin(), string_2.end());

        rimerge::RLEString::RLEncoder& encoder_1 = encoders.get_encoder(0);
        for (auto e : string_1)
        { encoder_1(e); }

        rimerge::RLEString::RLEncoder& encoder_2 = encoders.get_encoder(1);
        for (auto e : string_2)
        { encoder_2(e); }

        encoders.close();

        rimerge::RLEString rle_string;
        rle_string.load(testfiles_dir + "/rnd_cnct.seq");

        EXPECT_EQ(string.size(), rle_string.size());
        for (std::size_t i = 0; i < string.size(); i++)
        {
            EXPECT_EQ(string[i], rle_string[i]);
        }
    }

}

TEST(RLEString, RunCache_1)
{
    rimerge::string_type string = {'A','A','A','G','C','A','A','T','T','T','T','T','T','A','A','G','G','A','A','A'};
    rimerge::RLEString rle_string(string);

    rimerge::RLEString::RunCache cache(rle_string);

    for (std::size_t i = 0; i < string.size(); i++)
    {
        EXPECT_EQ(string[i], cache[i]);
    }
}

TEST(RLEString, RunCache_2)
{
    rimerge::string_type string = {'A', 'G', 'C', 'A', 'T', 'A', 'G', 'A'};
    rimerge::RLEString rle_string(string);

    rimerge::RLEString::RunCache cache(rle_string);

    for (std::size_t i = 0; i < string.size(); i++)
    {
        EXPECT_EQ(string[i], cache[i]);
    }
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

//------------------------------------------------------------------------------

TEST(RIndex, Inheritance)
{
    rimerge::RIndex<rimerge::RIndexRLE, rimerge::RLEString> index;

    EXPECT_EQ(std::string("RIndexRLE"), std::string(index.name()));
}


} // namespace