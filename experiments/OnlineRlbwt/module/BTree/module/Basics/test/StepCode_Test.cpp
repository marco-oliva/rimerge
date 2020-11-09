#include "StepCode.hpp"
#include <gtest/gtest.h>
#include <stdint.h>

namespace itmmti
{
  class StepCode_Test : public ::testing::Test
  {
  };


  TEST_F(StepCode_Test, Constructor)
  {
    const size_t num = 32;

    // { // Compare with different SizeT
    //   StepCode<16, uint64_t> code1;
    //   code1.printStatistics();
    //   StepCode<16, uint32_t> code2;
    //   code2.printStatistics();
    //   StepCode<16, uint16_t> code3;
    //   code3.printStatistics();
    // }

    using SC = StepCode<num, uint16_t>;
    SC code1(64);
    for (uint64_t j = 0; j < num; ++j) {
      const auto val = UINT64_C(1) << (j % 64);
      const auto newBvSize = code1.bitSize() + StepCodeUtil::calcSteppedW(val);
      if (newBvSize > code1.bitCapacity()) {
        code1.changeBitCapacity(static_cast<size_t>(newBvSize * 1.25));
      }
      code1.append(val);
    }
    // code1.printStatistics(true);

    {
      SC code2(code1);
      ASSERT_EQ(code1.size(), code2.size());
      ASSERT_EQ(code1.bitSize(), code2.bitSize());
      for (uint64_t j = 0; j < num; ++j) {
        ASSERT_EQ(code1.read(j), code2.read(j));
      }
    }
    {
      SC code2(num/2);
      code2 = code1;
      ASSERT_EQ(code1.size(), code2.size());
      ASSERT_EQ(code1.bitSize(), code2.bitSize());
      for (uint64_t j = 0; j < num; ++j) {
        ASSERT_EQ(code1.read(j), code2.read(j));
      }
    }
    {
      SC code_copy(code1);
      SC code2(std::move(code_copy));
      ASSERT_EQ(0, code_copy.size());
      ASSERT_EQ(0, code_copy.bitCapacity());
      ASSERT_EQ(code1.size(), code2.size());
      ASSERT_EQ(code1.bitSize(), code2.bitSize());
      for (uint64_t j = 0; j < num; ++j) {
        ASSERT_EQ(code1.read(j), code2.read(j));
      }
    }
    {
      SC code_copy(code1);
      SC code2(num/2);
      code2 = std::move(code_copy);
      ASSERT_EQ(0, code_copy.size());
      ASSERT_EQ(0, code_copy.bitCapacity());
      ASSERT_EQ(code1.size(), code2.size());
      ASSERT_EQ(code1.bitSize(), code2.bitSize());
      for (uint64_t j = 0; j < num; ++j) {
        ASSERT_EQ(code1.read(j), code2.read(j));
      }
    }
  }


  TEST_F(StepCode_Test, AppendRead)
  {
    const size_t num = 128;

    using SC = StepCode<num, uint16_t>;
    SC code(64);
    for (uint64_t j = 0; j < num; ++j) {
      const auto val = UINT64_C(1) << (j % 64);
      const auto newBvSize = code.bitSize() + StepCodeUtil::calcSteppedW(val);
      if (newBvSize > code.bitCapacity()) {
        code.changeBitCapacity(static_cast<size_t>(newBvSize * 1.25));
      }
      code.append(val);
    }
    // code.printStatistics(true);
    uint64_t pos = 0;
    for (uint64_t j = 0; j < num; ++j) {
      const auto w = code.readW(j);
      ASSERT_EQ(UINT64_C(1) << (j % 64), code.readWBits(pos, w));
      pos += w;
    }
    // code.changeBitCapacity();
    // code.printStatistics(true);
  }


  TEST_F(StepCode_Test, ChangeWCodes)
  {
    const size_t num = 32;

    using SC = StepCode<num, uint16_t>;
    SC code(64);
    for (uint64_t j = 0; j < num - 16; ++j) {
      const auto val = UINT64_C(1) << (j % 64);
      const auto newBvSize = code.bitSize() + StepCodeUtil::calcSteppedW(val);
      if (newBvSize > code.bitCapacity()) {
        code.changeBitCapacity(static_cast<size_t>(newBvSize * 1.25));
      }
      code.append(val);
    }
    code.printStatistics(true);
    uint64_t wCodeSrc = 0;
    for (uint8_t j = 0; j < 16; ++j) {
      StepCodeUtil::writeWCode(j, &wCodeSrc, j);
    }

    {
      SC cc(code);
      const uint64_t srcIdxBeg = 0;
      const uint64_t srcIdxLen = 5;
      const uint64_t tgtIdxBeg = 0;
      const uint64_t tgtIdxLen = 0;
      const uint64_t bitPos = cc.calcBitPos(tgtIdxBeg);
      const auto insBitLen = StepCodeUtil::sumW(&wCodeSrc, srcIdxBeg, srcIdxBeg + srcIdxLen);
      const auto delBitLen = StepCodeUtil::sumW(cc.getConstPtr_wCodes(), tgtIdxBeg, tgtIdxBeg + tgtIdxLen);
      const auto newBvSize = cc.bitSize() + insBitLen - delBitLen;
      if (insBitLen > delBitLen && newBvSize > cc.bitCapacity()) {
        cc.changeBitCapacity(newBvSize);
      }
      std::cout << "srcIdxBeg = " << srcIdxBeg << ", srcIdxLen = " << srcIdxLen << ", tgtIdxBeg = " << tgtIdxBeg << ", tgtIdxLen = "
                << tgtIdxLen << ", bitPos = " << bitPos << ", insBitLen = " << insBitLen << ", delBitLen = " << delBitLen << std::endl;
      cc.changeWCodesAndValPos(&wCodeSrc, srcIdxBeg, srcIdxLen, tgtIdxBeg, tgtIdxLen, bitPos, insBitLen, delBitLen);
      cc.printStatistics(true);
    }
    {
      SC cc(code);
      const uint64_t srcIdxBeg = 0;
      const uint64_t srcIdxLen = 16;
      const uint64_t tgtIdxBeg = 0;
      const uint64_t tgtIdxLen = 3;
      const uint64_t bitPos = cc.calcBitPos(tgtIdxBeg);
      const auto insBitLen = StepCodeUtil::sumW(&wCodeSrc, srcIdxBeg, srcIdxBeg + srcIdxLen);
      const auto delBitLen = StepCodeUtil::sumW(cc.getConstPtr_wCodes(), tgtIdxBeg, tgtIdxBeg + tgtIdxLen);
      const auto newBvSize = cc.bitSize() + insBitLen - delBitLen;
      if (insBitLen > delBitLen && newBvSize > cc.bitCapacity()) {
        cc.changeBitCapacity(newBvSize);
      }
      std::cout << "srcIdxBeg = " << srcIdxBeg << ", srcIdxLen = " << srcIdxLen << ", tgtIdxBeg = " << tgtIdxBeg << ", tgtIdxLen = "
                << tgtIdxLen << ", bitPos = " << bitPos << ", insBitLen = " << insBitLen << ", delBitLen = " << delBitLen << std::endl;
      cc.changeWCodesAndValPos(&wCodeSrc, srcIdxBeg, srcIdxLen, tgtIdxBeg, tgtIdxLen, bitPos, insBitLen, delBitLen);
      cc.printStatistics(true);
    }
    {
      SC cc(code);
      const uint64_t srcIdxBeg = 0;
      const uint64_t srcIdxLen = 8;
      const uint64_t tgtIdxBeg = 8;
      const uint64_t tgtIdxLen = 1;
      const uint64_t bitPos = cc.calcBitPos(tgtIdxBeg);
      const auto insBitLen = StepCodeUtil::sumW(&wCodeSrc, srcIdxBeg, srcIdxBeg + srcIdxLen);
      const auto delBitLen = StepCodeUtil::sumW(cc.getConstPtr_wCodes(), tgtIdxBeg, tgtIdxBeg + tgtIdxLen);
      const auto newBvSize = cc.bitSize() + insBitLen - delBitLen;
      if (insBitLen > delBitLen && newBvSize > cc.bitCapacity()) {
        cc.changeBitCapacity(newBvSize);
      }
      std::cout << "srcIdxBeg = " << srcIdxBeg << ", srcIdxLen = " << srcIdxLen << ", tgtIdxBeg = " << tgtIdxBeg << ", tgtIdxLen = "
                << tgtIdxLen << ", bitPos = " << bitPos << ", insBitLen = " << insBitLen << ", delBitLen = " << delBitLen << std::endl;
      cc.changeWCodesAndValPos(&wCodeSrc, srcIdxBeg, srcIdxLen, tgtIdxBeg, tgtIdxLen, bitPos, insBitLen, delBitLen);
      cc.printStatistics(true);
    }
    { // delete
      SC cc(code);
      const uint64_t srcIdxBeg = 0;
      const uint64_t srcIdxLen = 0;
      const uint64_t tgtIdxBeg = 7;
      const uint64_t tgtIdxLen = 4;
      const uint64_t bitPos = cc.calcBitPos(tgtIdxBeg);
      const auto insBitLen = StepCodeUtil::sumW(&wCodeSrc, srcIdxBeg, srcIdxBeg + srcIdxLen);
      const auto delBitLen = StepCodeUtil::sumW(cc.getConstPtr_wCodes(), tgtIdxBeg, tgtIdxBeg + tgtIdxLen);
      const auto newBvSize = cc.bitSize() + insBitLen - delBitLen;
      if (insBitLen > delBitLen && newBvSize > cc.bitCapacity()) {
        cc.changeBitCapacity(newBvSize);
      }
      std::cout << "srcIdxBeg = " << srcIdxBeg << ", srcIdxLen = " << srcIdxLen << ", tgtIdxBeg = " << tgtIdxBeg << ", tgtIdxLen = "
                << tgtIdxLen << ", bitPos = " << bitPos << ", insBitLen = " << insBitLen << ", delBitLen = " << delBitLen << std::endl;
      cc.changeWCodesAndValPos(&wCodeSrc, srcIdxBeg, srcIdxLen, tgtIdxBeg, tgtIdxLen, bitPos, insBitLen, delBitLen);
      cc.printStatistics(true);
    }
  }

} // namespace



