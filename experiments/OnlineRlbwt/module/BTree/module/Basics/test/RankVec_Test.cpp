#include "RankVec.hpp"
#include <gtest/gtest.h>
#include <stdint.h>

namespace itmmti
{
  class RankVec_Test : public ::testing::Test
  {
  };

  TEST_F(RankVec_Test, Constructor)
  {
    size_t num = 8200;

    // {
    //   RankVec64<> v64;
    //   v64.printStatistics();
    //   RankVec32<> v32;
    //   v32.printStatistics();
    // }

    RankVec<> v1(num);
    {
      uint64_t sum = 0;
      for (uint64_t j = 0; j < num; ++j) {
        const auto diff0 = (j % 1000);
        if (v1.capacity() - v1.size() < diff0 + 1) {
          // v1.changeCapacity((v1.size() + diff0 + 1) * 2); // reserve "needed size * 2"
          v1.changeCapacity(((v1.size() + diff0 + 1) * 3) / 2); // reserve "needed size * 1.5"
        }
        for (uint64_t i = 0; i < diff0; ++i) {
          v1.appendBit(false);
        }
        v1.appendBit(true);
        sum += diff0 + 1;
      }
      ASSERT_EQ(sum, v1.size());
      ASSERT_LE(sum, v1.capacity());
    }
    // v1.printStatistics();

    {
      RankVec<> v2(v1);
      ASSERT_EQ(v1.size(), v2.size());
      for (uint64_t j = 0; j < num; ++j) {
        ASSERT_EQ(v1.readBit(j), v2.readBit(j));
        ASSERT_EQ(v1.rank_1(j), v2.rank_1(j));
        ASSERT_EQ(v1.rank_0(j), v2.rank_0(j));
      }
    }
    {
      RankVec<> v2(num/2);
      v2 = v1;
      ASSERT_EQ(v1.size(), v2.size());
      for (uint64_t j = 0; j < num; ++j) {
        ASSERT_EQ(v1.readBit(j), v2.readBit(j));
        ASSERT_EQ(v1.rank_1(j), v2.rank_1(j));
        ASSERT_EQ(v1.rank_0(j), v2.rank_0(j));
      }
    }
    {
      RankVec<> v_copy(v1);
      RankVec<> v2(std::move(v_copy));
      ASSERT_EQ(0, v_copy.size());
      ASSERT_EQ(0, v_copy.capacity());
      ASSERT_EQ(v1.size(), v2.size());
      for (uint64_t j = 0; j < num; ++j) {
        ASSERT_EQ(v1.readBit(j), v2.readBit(j));
        ASSERT_EQ(v1.rank_1(j), v2.rank_1(j));
        ASSERT_EQ(v1.rank_0(j), v2.rank_0(j));
      }
    }
    {
      RankVec<> v_copy(v1);
      RankVec<> v2(num/2);
      v2 = std::move(v_copy);
      ASSERT_EQ(0, v_copy.size());
      ASSERT_EQ(0, v_copy.capacity());
      ASSERT_EQ(v1.size(), v2.size());
      for (uint64_t j = 0; j < num; ++j) {
        ASSERT_EQ(v1.readBit(j), v2.readBit(j));
        ASSERT_EQ(v1.rank_1(j), v2.rank_1(j));
        ASSERT_EQ(v1.rank_0(j), v2.rank_0(j));
      }
    }
  }


  TEST_F(RankVec_Test, RankSelect0)
  {
    size_t num = 8200;

    RankVec<> v(num);
    ASSERT_EQ(0, v.size());
    for (uint64_t j = 0; j < num; ++j) {
      const auto diff0 = (j % 1000);
      if (v.capacity() - v.size() < diff0 + 1) {
        // v.changeCapacity((v.size() + diff0 + 1) * 2); // reserve "needed size * 2"
        v.changeCapacity(((v.size() + diff0 + 1) * 3) / 2); // reserve "needed size * 1.5"
      }
      for (uint64_t i = 0; i < diff0; ++i) {
        v.appendBit(false);
      }
      v.appendBit(true);
    }
    // v.printStatistics(true);
    uint64_t count = 0;
    for (uint64_t j = 0; j < num; ++j) {
      const auto ans = v.rank_0(j);
      ASSERT_LE(count, ans);
      ASSERT_EQ(ans - count, 1 - v.readBit(j));
      if (ans > count) {
        ASSERT_EQ(j, v.select_0(ans));
        count = ans;
      }
    }
  }


  TEST_F(RankVec_Test, RankSelect1)
  {
    size_t num = 8200;

    RankVec<> v(num);
    ASSERT_EQ(0, v.size());
    for (uint64_t j = 0; j < num; ++j) {
      const auto diff0 = (j % 1000);
      if (v.capacity() - v.size() < diff0 + 1) {
        // v.changeCapacity((v.size() + diff0 + 1) * 2); // reserve "needed size * 2"
        v.changeCapacity(((v.size() + diff0 + 1) * 3) / 2); // reserve "needed size * 1.5"
      }
      for (uint64_t i = 0; i < diff0; ++i) {
        v.appendBit(false);
      }
      v.appendBit(true);
    }
    // v.printStatistics(true);
    uint64_t count = 0;
    for (uint64_t j = 0; j < num; ++j) {
      const auto ans = v.rank_1(j);
      ASSERT_LE(count, ans);
      ASSERT_EQ(ans - count, v.readBit(j));
      if (ans > count) {
        ASSERT_EQ(j, v.select_1(ans));
        count = ans;
      }
    }
  }


  TEST_F(RankVec_Test, PredSucc0)
  {
    size_t num = 8200;

    RankVec<> v(num);
    ASSERT_EQ(0, v.size());
    for (uint64_t j = 0; j < num; ++j) {
      const auto diff0 = (j % 1000);
      if (v.capacity() - v.size() < diff0 + 1) {
        // v.changeCapacity((v.size() + diff0 + 1) * 2); // reserve "needed size * 2"
        v.changeCapacity(((v.size() + diff0 + 1) * 3) / 2); // reserve "needed size * 1.5"
      }
      for (uint64_t i = 0; i < diff0; ++i) {
        v.appendBit(false);
      }
      v.appendBit(true);
    }
    // v.printStatistics(true);
    const auto num_zeros = v.getNum_0();
    { // pred
      uint64_t pos0 = 0;
      uint64_t pos1 = (num_zeros > 0)? v.select_0(1) : num;
      for ( ; pos0 < pos1; ++pos0) {
        ASSERT_EQ(UINT64_MAX, v.pred_0(pos0));
      }
      for (uint64_t r = 1; r < num_zeros; ++r) {
        pos1 = v.select_0(r+1);
        for ( ; pos0 < pos1; ++pos0) {
          ASSERT_EQ(v.select_0(r), v.pred_0(pos0));
        }
      }
      for ( ; pos0 < 2 * num; ++pos0) {
        ASSERT_EQ(pos1, v.pred_0(pos0));
      }
    }
    {
      uint64_t pos0 = 0;
      uint64_t pos1 = 0;
      for (uint64_t r = 1; r <= num_zeros; ++r) {
        pos1 = v.select_0(r);
        // std::cout << r << ": " << pos1 << std::endl;
        for ( ; pos0 <= pos1; ++pos0) {
          // std::cout << "pos0 = " << pos0 << ", pos1 = " << pos1 << std::endl;
          ASSERT_EQ(pos1, v.succ_0(pos0));
        }
      }
      for ( ; pos0 < 2 * num; ++pos0) {
        ASSERT_EQ(UINT64_MAX, v.succ_0(pos0));
      }
    }
  }


  TEST_F(RankVec_Test, PredSucc1)
  {
    size_t num = 8200;

    RankVec<> v(num);
    ASSERT_EQ(0, v.size());
    for (uint64_t j = 0; j < num; ++j) {
      const auto diff0 = (j % 1000);
      if (v.capacity() - v.size() < diff0 + 1) {
        // v.changeCapacity((v.size() + diff0 + 1) * 2); // reserve "needed size * 2"
        v.changeCapacity(((v.size() + diff0 + 1) * 3) / 2); // reserve "needed size * 1.5"
      }
      for (uint64_t i = 0; i < diff0; ++i) {
        v.appendBit(false);
      }
      v.appendBit(true);
    }
    // v.printStatistics(true);
    const auto num_ones = v.getNum_1();
    { // pred
      uint64_t pos0 = 0;
      uint64_t pos1 = (num_ones > 0)? v.select_1(1) : num;
      for ( ; pos0 < pos1; ++pos0) {
        ASSERT_EQ(UINT64_MAX, v.pred_1(pos0));
      }
      for (uint64_t r = 1; r < num_ones; ++r) {
        pos1 = v.select_1(r+1);
        for ( ; pos0 < pos1; ++pos0) {
          ASSERT_EQ(v.select_1(r), v.pred_1(pos0));
        }
      }
      for ( ; pos0 < 2 * num; ++pos0) {
        ASSERT_EQ(pos1, v.pred_1(pos0));
      }
    }
    {
      uint64_t pos0 = 0;
      uint64_t pos1 = 0;
      for (uint64_t r = 1; r <= num_ones; ++r) {
        pos1 = v.select_1(r);
        // std::cout << r << ": " << pos1 << std::endl;
        for ( ; pos0 <= pos1; ++pos0) {
          // std::cout << "pos0 = " << pos0 << ", pos1 = " << pos1 << std::endl;
          ASSERT_EQ(pos1, v.succ_1(pos0));
        }
      }
      for ( ; pos0 < 2 * num; ++pos0) {
        ASSERT_EQ(UINT64_MAX, v.succ_1(pos0));
      }
    }
  }


  TEST_F(RankVec_Test, Shorten)
  {
    size_t num = 8200;

    RankVec<> v(num);
    ASSERT_EQ(0, v.size());
    for (uint64_t j = 0; j < num; ++j) {
      const auto diff0 = (j % 3);
      if (v.capacity() - v.size() < diff0 + 1) {
        // v.changeCapacity((v.size() + diff0 + 1) * 2); // reserve "needed size * 2"
        v.changeCapacity(((v.size() + diff0 + 1) * 3) / 2); // reserve "needed size * 1.5"
      }
      for (uint64_t i = 0; i < diff0; ++i) {
        v.appendBit(false);
      }
      v.appendBit(true);
    }
    // v.printStatistics(true);

    const uint64_t cutsize[] = {0, 1, 2, 3, 4, 32, 256, 1024, 8191, 8192, 8193};
    for (uint32_t e = 0; e < sizeof(cutsize) / sizeof(cutsize[0]); ++e) {
      RankVec<> v_copy(v);
      if (v.size() < cutsize[e]) {
        continue;
      }
      const auto shortenSize = v.size() - cutsize[e];
      v_copy.shorten(shortenSize);
      ASSERT_EQ(shortenSize, v_copy.size());
      // v_copy.printStatistics(true);
      // check RankSelect1
      uint64_t count = 0;
      for (uint64_t j = 0; j < shortenSize; ++j) {
        const auto ans = v.rank_1(j);
        ASSERT_EQ(ans, v_copy.rank_1(j));
        ASSERT_EQ(ans - count, v_copy.readBit(j));
        if (ans > count) {
          ASSERT_EQ(j, v_copy.select_1(ans));
          count = ans;
        }
      }
    }

    const uint64_t size[] = {0, 1, 2, 3, 4, 32, 256, 257, 1024, 4095, 4096, 4097};
    for (uint32_t e = 0; e < sizeof(cutsize) / sizeof(cutsize[0]); ++e) {
      RankVec<> v_copy(v);
      if (v.size() < size[e]) {
        continue;
      }
      const auto shortenSize = size[e];
      v_copy.shorten(shortenSize);
      ASSERT_EQ(shortenSize, v_copy.size());
      // v_copy.printStatistics(true);
      // check RankSelect1
      uint64_t count = 0;
      for (uint64_t j = 0; j < shortenSize; ++j) {
        const auto ans = v.rank_1(j);
        ASSERT_EQ(ans, v_copy.rank_1(j));
        ASSERT_EQ(ans - count, v_copy.readBit(j));
        if (ans > count) {
          ASSERT_EQ(j, v_copy.select_1(ans));
          count = ans;
        }
      }
    }
  }
} // namespace



