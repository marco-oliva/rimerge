#include "WBitsBlockVec.hpp"
#include <gtest/gtest.h>
#include <stdint.h>

namespace itmmti
{
  using WBitsBlockVecT = WBitsBlockVec<1024>;

  class WBitsBlockVecTest : public ::testing::Test
  {
  };

  TEST_F(WBitsBlockVecTest, Constructor)
  {
    size_t num = 2000;
    uint8_t w = 61;

    WBitsBlockVecT wbv1(w, num);
    wbv1.resize(num);
    for (uint64_t j = 0; j < num; ++j) {
      wbv1.write((UINT64_C(1) << (j % 64)) & bits::UINTW_MAX(w), j);
    }
    {
      WBitsBlockVecT wbv2(wbv1);
      ASSERT_EQ(wbv1.capacity(), wbv2.capacity());
      ASSERT_EQ(wbv1.size(), wbv2.size());
      ASSERT_EQ(wbv1.getW(), wbv2.getW());
      for (uint64_t j = 0; j < num; ++j) {
        ASSERT_EQ(wbv1.read(j), wbv2.read(j));
      }
    }
    {
      WBitsBlockVecT wbv2(w/2, num/2);
      wbv2 = wbv1;
      ASSERT_EQ(wbv1.capacity(), wbv2.capacity());
      ASSERT_EQ(wbv1.size(), wbv2.size());
      ASSERT_EQ(wbv1.getW(), wbv2.getW());
      for (uint64_t j = 0; j < num; ++j) {
        ASSERT_EQ(wbv1.read(j), wbv2.read(j));
      }
    }
    {
      WBitsBlockVecT wbv_copy(wbv1);
      WBitsBlockVecT wbv2(std::move(wbv_copy));
      ASSERT_EQ(0, wbv_copy.size());
      ASSERT_EQ(0, wbv_copy.capacity());
      ASSERT_EQ(wbv1.capacity(), wbv2.capacity());
      ASSERT_EQ(wbv1.size(), wbv2.size());
      ASSERT_EQ(wbv1.getW(), wbv2.getW());
      for (uint64_t j = 0; j < num; ++j) {
        ASSERT_EQ(wbv1.read(j), wbv2.read(j));
      }
    }
    {
      WBitsBlockVecT wbv_copy(wbv1);
      WBitsBlockVecT wbv2(w/2, num/2);
      wbv2 = std::move(wbv_copy);
      ASSERT_EQ(0, wbv_copy.size());
      ASSERT_EQ(0, wbv_copy.capacity());
      ASSERT_EQ(wbv1.capacity(), wbv2.capacity());
      ASSERT_EQ(wbv1.size(), wbv2.size());
      ASSERT_EQ(wbv1.getW(), wbv2.getW());
      for (uint64_t j = 0; j < num; ++j) {
        ASSERT_EQ(wbv1.read(j), wbv2.read(j));
      }
    }
  }


  TEST_F(WBitsBlockVecTest, ReadWrite)
  {
    size_t num = 2000;

    for (uint8_t w = 1; w <= 64; ++w) {
      WBitsBlockVecT wbv(w, num);
      wbv.resize(num);
      // std::cout << "calcMemBytes [w = " << static_cast<uint32_t>(wbv.getW()) << "]: " << wbv.calcMemBytes()
      //           << "bytes, size / capacity_: " << wbv.size() << " / " << wbv.capacity() << std::endl;
      for (uint64_t j = 0; j < num; ++j) {
        wbv.write((UINT64_C(1) << (j % 64)) & bits::UINTW_MAX(w), j);
      }
      for (uint64_t j = 0; j < num; ++j) {
        ASSERT_EQ((UINT64_C(1) << (j % 64)) & bits::UINTW_MAX(w), wbv.read(j));
      }
    }
  }


  TEST_F(WBitsBlockVecTest, ReadWriteSubscript)
  {
    size_t num = 2000;

    for (uint8_t w = 1; w <= 64; ++w) {
      WBitsBlockVecT wbv(w, num);
      wbv.resize(num);
      // std::cout << "calcMemBytes [w = " << static_cast<uint32_t>(wbv.getW()) << "]: " << wbv.calcMemBytes()
      //           << "bytes, size / capacity_: " << wbv.size() << " / " << wbv.capacity() << std::endl;
      for (uint64_t j = 0; j < num; ++j) {
        wbv[j] = (UINT64_C(1) << (j % 64)) & bits::UINTW_MAX(w);
      }
      for (uint64_t j = 0; j < num; ++j) {
        ASSERT_EQ((UINT64_C(1) << (j % 64)) & bits::UINTW_MAX(w), wbv[j]);
      }
    }
  }


  TEST_F(WBitsBlockVecTest, ChangeW)
  {
    size_t num = 2000;

    for (uint8_t step = 1; step < 63; ++step) {
      for (uint8_t w = 1; w + step <= 64; ++w) {
        WBitsBlockVecT wbv(w, num);
        wbv.resize(num);
        for (uint64_t j = 0; j < num; ++j) {
          wbv.write((UINT64_C(1) << (j % 64)) & bits::UINTW_MAX(w), j);
        }
        // std::cout << "Convert " << static_cast<uint32_t>(wbv.getW()) << " (" << wbv.calcMemBytes() << " bytes) -> ";
        wbv.increaseW(w+step); // convert 'w' bits vec to 'w+step' bits vec
        // std::cout << static_cast<uint32_t>(wbv.getW()) << " (" << wbv.calcMemBytes() << " bytes)" << std::endl;
        for (uint64_t j = 0; j < num; ++j) {
          ASSERT_EQ((UINT64_C(1) << (j % 64)) & bits::UINTW_MAX(w), wbv.read(j));
        }
      }
    }

    for (uint8_t step = 1; step < 63; ++step) {
      for (uint8_t w = 64; w > step; --w) {
        WBitsBlockVecT wbv(w, num);
        wbv.resize(num);
        for (uint64_t j = 0; j < num; ++j) {
          wbv.write((UINT64_C(1) << (j % 64)) & bits::UINTW_MAX(w), j);
        }
        // std::cout << "Convert " << static_cast<uint32_t>(wbv.getW()) << " (" << wbv.calcMemBytes() << " bytes) -> ";
        wbv.decreaseW(w-step); // convert 'w' bits vec to 'w-step' bits vec
        // std::cout << static_cast<uint32_t>(wbv.getW()) << " (" << wbv.calcMemBytes() << " bytes)" << std::endl;
        for (uint64_t j = 0; j < num; ++j) {
          ASSERT_EQ((UINT64_C(1) << (j % 64)) & bits::UINTW_MAX(w-step), wbv.read(j));
        }
      }
    }
  }
} // anonymous namespace
