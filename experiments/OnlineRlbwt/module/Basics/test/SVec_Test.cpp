#include "SVec.hpp"
#include "RankVec.hpp"
#include <gtest/gtest.h>
#include <stdint.h>

namespace itmmti
{
  class SVec_Test : public ::testing::Test
  {
  };

  /*
  TEST_F(RankVec_Test, Constructor)
  {
    size_t num = 2000;

    RankVec bv1(num);
    bv1.resize(num);
    for (uint64_t j = 0; j < num; ++j) {
      bv1.writeBit((UINT64_C(1) << (j % 64)) & bits::UINTW_MAX(1), j);
    }
    {
      BitVec bv2(bv1);
      ASSERT_EQ(bv1.capacity(), bv2.capacity());
      ASSERT_EQ(bv1.size(), bv2.size());
      for (uint64_t j = 0; j < num; ++j) {
        ASSERT_EQ(bv1.readBit(j), bv2.readBit(j));
      }
    }
    {
      BitVec bv2(num/2);
      bv2 = bv1;
      ASSERT_EQ(bv1.capacity(), bv2.capacity());
      ASSERT_EQ(bv1.size(), bv2.size());
      for (uint64_t j = 0; j < num; ++j) {
        ASSERT_EQ(bv1.readBit(j), bv2.readBit(j));
      }
    }
    {
      BitVec bv_copy(bv1);
      BitVec bv2(std::move(bv_copy));
      ASSERT_EQ(0, bv_copy.size());
      ASSERT_EQ(0, bv_copy.capacity());
      ASSERT_EQ(bv1.capacity(), bv2.capacity());
      ASSERT_EQ(bv1.size(), bv2.size());
      for (uint64_t j = 0; j < num; ++j) {
        ASSERT_EQ(bv1.readBit(j), bv2.readBit(j));
      }
    }
    {
      BitVec bv_copy(bv1);
      BitVec bv2(num/2);
      bv2 = std::move(bv_copy);
      ASSERT_EQ(0, bv_copy.size());
      ASSERT_EQ(0, bv_copy.capacity());
      ASSERT_EQ(bv1.capacity(), bv2.capacity());
      ASSERT_EQ(bv1.size(), bv2.size());
      for (uint64_t j = 0; j < num; ++j) {
        ASSERT_EQ(bv1.readBit(j), bv2.readBit(j));
      }
    }
  }
  */

  TEST_F(SVec_Test, RankSelect0)
  {
    size_t num = 8200;

    SVec<RankVec<> > v(8, num);
    ASSERT_EQ(0, v.size());
    uint64_t val = 0;
    for (uint64_t j = 0; j < num; ++j) {
      // val += 1;
      val += (j % 1000) + 1;
      v.append(val);
      ASSERT_EQ(val, v.select_1(j+1));
    }
    ASSERT_EQ(num, v.size());
    // v.printDebugInfo();
    uint64_t count = 0;
    for (uint64_t j = 0; j < num; ++j) {
      const auto ans = v.rank_0(j);
      ASSERT_LE(count, ans);
      if (ans > count) {
        ASSERT_EQ(ans - count, 1);
        ASSERT_EQ(j, v.select_0(ans));
        count = ans;
      }
    }
  }


  TEST_F(SVec_Test, RankSelect1)
  {
    size_t num = 8200;

    SVec<RankVec<> > v(8, num);
    ASSERT_EQ(0, v.size());
    uint64_t val = 0;
    for (uint64_t j = 0; j < num; ++j) {
      // val += 1;
      val += (j % 1000) + 1;
      v.append(val);
      ASSERT_EQ(val, v.select_1(j+1));
    }
    ASSERT_EQ(num, v.size());
    // v.printDebugInfo();
    uint64_t count = 0;
    for (uint64_t j = 0; j < num; ++j) {
      const auto ans = v.rank_1(j);
      ASSERT_LE(count, ans);
      if (ans > count) {
        ASSERT_EQ(ans - count, 1);
        ASSERT_EQ(j, v.select_1(ans));
        count = ans;
      }
    }
  }


  TEST_F(SVec_Test, ConvertDec)
  {
    size_t num = 8200;

    SVec<RankVec<> > sv(12, num);
    ASSERT_EQ(0, sv.size());
    uint64_t val = 0;
    for (uint64_t j = 0; j < num; ++j) {
      // val += 1;
      val += (j % 1000) + 1;
      sv.append(val);
      ASSERT_EQ(val, sv.select_1(j+1));
    }
    ASSERT_EQ(num, sv.size());
    SVec<RankVec<> > sv1(sv);
    sv.printStatistics();
    std::cout << "shrink_to_fit() ->" << std::endl;
    sv.shrink_to_fit();
    sv.printStatistics();
    {
      const auto optLoW = sv1.calcOptimalLoW(sv1.getMax(), num);
      std::cout << "Optimize " << (uint64_t)(sv1.getLoW()) << " to " << (uint64_t)(optLoW) << std::endl;
      sv1.convert(optLoW, 0, 1.0, true); // Convert to optimal loW with shrinking.
    }
    sv1.printStatistics();

    const auto num_ones = sv.getNum_1();
    const auto num_zeros = sv.getNum_0();
    for (uint64_t i = 0; i < num; ++i) {
      ASSERT_EQ(sv.rank_1(i), sv1.rank_1(i));
      ASSERT_EQ(sv.rank_0(i), sv1.rank_0(i));
    }
    for (uint64_t i = 1; i <= num_ones; ++i) {
      ASSERT_EQ(sv.select_1(i), sv1.select_1(i));
    }
    for (uint64_t i = 1; i <= num_zeros; ++i) {
      ASSERT_EQ(sv.select_0(i), sv1.select_0(i));
    }
  }


  TEST_F(SVec_Test, ConvertInc)
  {
    size_t num = 8200;

    SVec<RankVec<> > sv(4, num);
    ASSERT_EQ(0, sv.size());
    uint64_t val = 0;
    for (uint64_t j = 0; j < num; ++j) {
      // val += 1;
      val += (j % 1000) + 1;
      sv.append(val);
      ASSERT_EQ(val, sv.select_1(j+1));
    }
    ASSERT_EQ(num, sv.size());
    SVec<RankVec<> > sv1(sv);
    sv.printStatistics();
    std::cout << "shrink_to_fit() ->" << std::endl;
    sv.shrink_to_fit();
    sv.printStatistics();
    {
      const auto optLoW = sv1.calcOptimalLoW(sv1.getMax(), num);
      std::cout << "Optimize " << (uint64_t)(sv1.getLoW()) << " to " << (uint64_t)(optLoW) << std::endl;
      sv1.convert(optLoW, 0, 1.0, true); // Convert to optimal loW with shrinking.
    }
    sv1.printStatistics();

    const auto num_ones = sv.getNum_1();
    const auto num_zeros = sv.getNum_0();
    for (uint64_t i = 0; i < num; ++i) {
      ASSERT_EQ(sv.rank_1(i), sv1.rank_1(i));
      ASSERT_EQ(sv.rank_0(i), sv1.rank_0(i));
    }
    for (uint64_t i = 1; i <= num_ones; ++i) {
      ASSERT_EQ(sv.select_1(i), sv1.select_1(i));
    }
    for (uint64_t i = 1; i <= num_zeros; ++i) {
      ASSERT_EQ(sv.select_0(i), sv1.select_0(i));
    }
  }


#if 0
  TEST_F(SVec_Test, PredSucc0)
  {
    size_t num = 16384;

    RankVec<> rv(num);
    // std::cout << "calcMemBytes: " << bv.calcMemBytes()
    //           << "bytes, size / capacity_: " << bv.size() << " / " << bv.capacity() << std::endl;
    ASSERT_EQ(0, rv.size());
    for (uint64_t j = 0; j < num; ++j) {
      rv.appendBit(!(j % 571) || !(j % 285));
    }
    ASSERT_EQ(num, rv.size());
    // rv.printDebugInfo();
    const auto num_zeros = rv.getNum_0();
    { // pred
      uint64_t pos0 = 0;
      uint64_t pos1 = (num_zeros > 0)? rv.select_0(1) : num;
      for ( ; pos0 < pos1; ++pos0) {
        ASSERT_EQ(UINT64_MAX, rv.pred_0(pos0));
      }
      for (uint64_t r = 1; r < num_zeros; ++r) {
        pos1 = rv.select_0(r+1);
        for ( ; pos0 < pos1; ++pos0) {
          ASSERT_EQ(rv.select_0(r), rv.pred_0(pos0));
        }
      }
      for ( ; pos0 < 2 * num; ++pos0) {
        ASSERT_EQ(pos1, rv.pred_0(pos0));
      }
    }
    {
      uint64_t pos0 = 0;
      uint64_t pos1 = 0;
      for (uint64_t r = 1; r <= num_zeros; ++r) {
        pos1 = rv.select_0(r);
        // std::cout << r << ": " << pos1 << std::endl;
        for ( ; pos0 <= pos1; ++pos0) {
          // std::cout << "pos0 = " << pos0 << ", pos1 = " << pos1 << std::endl;
          ASSERT_EQ(pos1, rv.succ_0(pos0));
        }
      }
      for ( ; pos0 < 2 * num; ++pos0) {
        ASSERT_EQ(UINT64_MAX, rv.succ_0(pos0));
      }
    }
  }


  TEST_F(RankVec_Test, PredSucc1)
  {
    size_t num = 16384;

    RankVec<> rv(num);
    // std::cout << "calcMemBytes: " << bv.calcMemBytes()
    //           << "bytes, size / capacity_: " << bv.size() << " / " << bv.capacity() << std::endl;
    ASSERT_EQ(0, rv.size());
    for (uint64_t j = 0; j < num; ++j) {
      rv.appendBit(!(j % 571) || !(j % 285));
    }
    ASSERT_EQ(num, rv.size());
    // rv.printDebugInfo();
    const auto num_ones = rv.getNum_1();
    { // pred
      uint64_t pos0 = 0;
      uint64_t pos1 = (num_ones > 0)? rv.select_1(1) : num;
      for ( ; pos0 < pos1; ++pos0) {
        ASSERT_EQ(UINT64_MAX, rv.pred_1(pos0));
      }
      for (uint64_t r = 1; r < num_ones; ++r) {
        pos1 = rv.select_1(r+1);
        for ( ; pos0 < pos1; ++pos0) {
          ASSERT_EQ(rv.select_1(r), rv.pred_1(pos0));
        }
      }
      for ( ; pos0 < 2 * num; ++pos0) {
        ASSERT_EQ(pos1, rv.pred_1(pos0));
      }
    }
    {
      uint64_t pos0 = 0;
      uint64_t pos1 = 0;
      for (uint64_t r = 1; r <= num_ones; ++r) {
        pos1 = rv.select_1(r);
        // std::cout << r << ": " << pos1 << std::endl;
        for ( ; pos0 <= pos1; ++pos0) {
          // std::cout << "pos0 = " << pos0 << ", pos1 = " << pos1 << std::endl;
          ASSERT_EQ(pos1, rv.succ_1(pos0));
        }
      }
      for ( ; pos0 < 2 * num; ++pos0) {
        ASSERT_EQ(UINT64_MAX, rv.succ_1(pos0));
      }
    }
  }
#endif
} // namespace



