//
//  unit_tests.cpp
//
// Copyright (c) Boucher Lab. All rights reserved.
// Licensed under the GNU license. See LICENSE file in the repository root for full license information.

#include <random>

#include <gtest/gtest.h>

#include <rimerge/utils.hpp>
#include <rimerge/rle_string.hpp>
#include <rimerge/r-index-rle.hpp>



namespace
{

std::string testfiles_dir = "../../tests/files";

rimerge::string_type random_string(rimerge::size_type length, rimerge::string_type seed = {})
{
    std::random_device rd;
    std::mt19937 generator(rd());
    
    if (seed.empty())
    {
        std::string str;
        str.append(length / 5, 'A'); str.append(length / 5, 'C'); str.append(length / 5, 'G');
        str.append(length / 5, 'T'); str.append(length / 5, 'N');
        std::shuffle(str.begin(), str.end(), generator);
        str.push_back('$');
        rimerge::string_type out; for (auto& c : str) { out.push_back(c); }
        return out;
    }
    else
    {
        std::uniform_real_distribution<> dist(0, 1);
        std::uniform_real_distribution<> dist_char(0, 5);
        std::string dict = "ACGTN";
        
        double error_rate = 0.15;
        rimerge::string_type str;
        for (auto& c : seed)
        {
            double rand = dist(generator);
            if (rand < error_rate) { str.push_back(dict[std::floor(dist_char(generator))]); }
            else { str.push_back(c); }
        }
        str.push_back('$');
        return str;
    }
    
}

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
    EXPECT_EQ(rimerge::DATA_TERMINATOR, alphabet.previous('A'));

    EXPECT_EQ('C', alphabet.following('A'));
    EXPECT_EQ('G', alphabet.following('C'));
    EXPECT_EQ(rimerge::DATA_TERMINATOR, alphabet.following('G'));
}

//------------------------------------------------------------------------------

TEST(RIndex, Inheritance)
{
    rimerge::RIndex<rimerge::RIndexRLE, rimerge::RLEString> index;

    EXPECT_EQ(std::string("RIndexRLE"), std::string(index.name()));
}

//------------------------------------------------------------------------------

void
build_csa(sdsl::csa_bitcompressed<>& out_csa, rimerge::string_type s)
{
    std::string tmp;
    tmp.insert(tmp.end(), s.begin(), s.end());
    sdsl::construct_im(out_csa, tmp.c_str(), 1);
}

void
csa_to_ri(sdsl::csa_bitcompressed<>& in_csa, std::string out_prefix)
{
    rimerge::RLEString::RLEncoder encoder(testfiles_dir + "/" + out_prefix + ".rle");
    std::ofstream samples(testfiles_dir + "/" + out_prefix + ".saes");
    for (std::size_t i = 0; i < in_csa.size(); i++)
    {
        if ((i == 0 or (i == 1 or (i == in_csa.size() - 1))) or (in_csa.bwt[i] != in_csa.bwt[i - 1] or (in_csa.bwt[i] != in_csa.bwt[i + 1])))
        {
            rimerge::SA_samples::sample_type t = {i, in_csa[i]};
            rimerge::SA_samples::write(t, &samples);
        }
        // replace 0 with data terminator
        if (in_csa.bwt[i] != '\0') { encoder(in_csa.bwt[i]); }
        else { encoder(rimerge::DATA_TERMINATOR); }
        
    }
}

template<typename I, typename BWT>
void
read_from_prefix(std::string prefix, rimerge::RIndex<I, BWT>& index)
{
    if (prefix[prefix.size() - 1] == '/') { index.init(prefix + "bwt.rle", prefix + "samples.saes"); }
    else { index.init(prefix + ".rle", prefix + ".saes"); }
}

TEST(TwoStrings, Rimerge)
{
    // r-index type
    typedef rimerge::RIndexRLE RIndexType;
    typedef rimerge::RLEString BWTType;
    
    for (std::size_t i = 0; i < 1; i++)
    {
        // Generate two random r-indxes
    
        // sequences, DATA_TERMINATOR is pushed by sdsl. (It pushes 0 and we substitute it with DATA_TERMINATOR)
        rimerge::string_type s1, s2, sm, seed;
        seed = random_string(1000); seed.pop_back(); // remove terminator
        s1 = random_string(0, seed); s1.pop_back(); // remove terminator
        s2 = random_string(0, seed); s2.pop_back();
        sm.insert(sm.end(), s2.begin(), s2.end());
        sm.push_back(rimerge::STRING_TERMINATOR);
        sm.insert(sm.end(), s1.begin(), s1.end());
    
        // compute bwt and samples for both sequences and for merged sequence
        sdsl::csa_bitcompressed<> csa1, csa2, csam;
        build_csa(csa1, s1);
        build_csa(csa2, s2);
        build_csa(csam, sm);
    
        // r-indexes on disk
        csa_to_ri(csa1, "s1");
        csa_to_ri(csa2, "s2");
    
        // Merge s1 and s2
        rimerge::RIndex<RIndexType, BWTType> r1, r2, rm;
        read_from_prefix(testfiles_dir + "/s1", r1);
        read_from_prefix(testfiles_dir + "/s2", r2);
    
        rimerge::MergeParameters merge_parameters;
        merge_parameters.setMergeJobs(1);
        merge_parameters.setMergeBuffers(1);
        merge_parameters.setSearchJobs(1);
        merge_parameters.out_prefix = testfiles_dir + "/rm";
        rimerge::RIndex<RIndexType, BWTType>::merge(r1, r2, merge_parameters);
    
        read_from_prefix(testfiles_dir + "/rm/", rm);
    
        // Extract sequences from bwt
        rimerge::string_type o1 = rm.get_sequence(0);
        rimerge::string_type o2 = rm.get_sequence(1);
        EXPECT_EQ(s1, o1);
        EXPECT_EQ(s2, o2);
    
        // Check bwt
        rimerge::string_type bwt_string = rm.bwt().to_string();
        EXPECT_EQ(csam.bwt.size(), bwt_string.size());
        for (std::size_t i = 0; i < bwt_string.size(); i++)
        {
            EXPECT_TRUE((csam.bwt[i] == bwt_string[i]) or (rimerge::is_terminator(bwt_string[i])));
        }
    
        // Suffx Array Samples, structural correctness
        std::tuple<std::vector<rimerge::size_type>, std::vector<rimerge::size_type>, std::vector<rimerge::size_type>> errors;
        bool check = rimerge::RIndexRLE::check(rm, errors);
        EXPECT_TRUE(check);
    
        // Suffix Array Samples, values. Extract each sequence and check samples
        rimerge::size_type sa_value = s1.size() + s2.size() + 2 - 1;
        rimerge::size_type j = 1;
        while (not rimerge::is_terminator(rm.bwt()[j]))
        {
            if (rm.its(j) != rimerge::sample_genre::NOT) { EXPECT_EQ(sa_value, rm.samples()[j]); }
            j = rm.LF(j); sa_value -= 1;
        }
    
        sa_value = s1.size();
        j = 0;
        while (not rimerge::is_terminator(rm.bwt()[j]))
        {
            if (rm.its(j) != rimerge::sample_genre::NOT) { EXPECT_EQ(sa_value, rm.samples()[j]); }
            j = rm.LF(j); sa_value -= 1; // error when j = 206
        }
    };
}


} // namespace