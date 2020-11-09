/*!
 * Copyright (c) 2017 Tomohiro I
 *
 * This program is released under the MIT License.
 * http://opensource.org/licenses/mit-license.php
 */
/*!
 * @file BitsUtil.hpp
 * @brief Useful functions to work on bits.
 * @author Tomohiro I
 * @date 2017-02-23
 * 
 * @note
 *   For a uint64_t array words[0..n) of length n, let words[[0..64*n)) denote a bit vector defined as follows:
 *   The i-th bit of words[[0..64*n)) is the (i%64)-th (low significant) bit of words[i/64].
 *   We use the following terminology:
 *   - bit-pos: position (0base) in a bit vector.
 *   - bit-region: sub-vector of a bit vector.
 *   - offset: bit-pos in a word, i.e., offset is in [0, 64).
 */
#ifndef INCLUDE_GUARD_BitsUtil
#define INCLUDE_GUARD_BitsUtil

#include <stdint.h> // include uint64_t etc.
#include <assert.h>
#include <iostream>
#include <cstring>

namespace itmmti
{
  /*!
   * @namespace ctcbits
   * @brief Bits functions specialized to run in compile time.
   */
  namespace ctcbits
  {
    /*!
     * @brief Compile time specialized version of bits::UINTW_MAX.
     */
    constexpr uint64_t UINTW_MAX
    (
     uint8_t w //!< Bit-width in [0, 64].
     ) {
      if (w >= 64) {
        return UINT64_MAX;
      } else {
        return (UINT64_C(1) << w) - 1;
      }
    }
  } // namespace ctcbits


  /*!
   * @namespace bits
   * @brief Useful functions to work on bits.
   */
  namespace bits
  {
    inline void printBits
    (
     const uint64_t * words,
     const uint64_t n
     ) noexcept {
      for (uint64_t i = 0; i < n; ++i) {
        for (uint64_t j = 0; j < 64; ++j) {
          std::cout << ((*words >> (63 - j)) & 1);
        }
        std::cout << " ";
        ++words;
      }
      std::cout << std::endl;
    }


    /*!
     * @brief Given w in [0, 64], TBL_LoSet[w] is the 64bit uint such that 'w' low significant bits are raised and the other bits not.
     */
    extern const uint64_t TBL_LoSet[65];

    /*!
     * @brief Table to implement select queries for 8bit uint.
     * @note
     *   Given 8bit uint 'v' and 'rank' in [1, 8],
     *   TBL_Sel8[(rank << 8) + v] returns the bit-pos where the 'rank'-th 1 appears in v[[0..8]].
     */
    extern const uint8_t TBL_Sel8[256 * 9];

    /*!
     * @brief Table to implement popcount for 8bit uint (the first 256 entries of TBL_Sel8 are used).
     * @note
     *   Given 8bit uint 'v', TBL_Cnt8[v] returns the popcount(v).
     */
    extern const uint8_t * const TBL_Cnt8;


    /*!
     * @brief Return the maximum 64bit uint with 'w' bits.
     */
    inline uint64_t UINTW_MAX
    (
     uint8_t w //!< Bit-width in [0, 64].
     ) {
      assert(w <= 64);

      return TBL_LoSet[w];
    }


    /*!
     * @brief Return if 'val' is power of two.
     */
    constexpr bool isPowerOfTwo
    (
     uint64_t val
     ) {
      return !(val & (val - 1));
    }


    /*!
     * @brief Return popcount(val).
     */
    constexpr uint8_t cnt64
    (
     uint64_t val
     ) {
      return static_cast<uint8_t>(__builtin_popcountll(val));
    }


    /*!
     * @brief Return Num of bits needed to represent 'val' (where we think that 0 needs 1 bit).
     */
    constexpr uint8_t bitSize
    (
     uint64_t val
     ) {
      if(val) {
        return static_cast<uint8_t>(64 - __builtin_clzll(val));
      }
      return 1;
    }


    /*!
     * @brief Return num of leading zero bits.
     */
    constexpr uint8_t clz
    (
     uint64_t val
     ) {
      if (val) {
        return static_cast<uint8_t>(__builtin_clzll(val));
      }
      return 64;
    }


    /*!
     * @brief Return num of trailing zero bits.
     */
    constexpr uint8_t ctz
    (
     uint64_t val
     ) {
      if (val) {
        return static_cast<uint8_t>(__builtin_ctzll(val));
      }
      return 64;
    }


    /*!
     * @brief Return bit-pos of the 'rank'-th 1 in 'word'.
     */
    inline uint8_t sel64
    (
     uint64_t word, //!< 64bit uint
     uint64_t rank //!< in [1, 64].
     ) {
      assert(1 <= rank && rank <= 64);
      assert(cnt64(word) >= rank); //! @pre 'word' must contain at least 'rank' 1's.

      for (int j = 0; j < 7; ++j) {
        const uint8_t x = word & 0xff;
        const uint64_t c = TBL_Cnt8[x];
        if (rank <= c) {
          return static_cast<uint8_t>(j*8 + TBL_Sel8[(rank<<8) + x]);
        }
        rank -= c;
        word >>= 8;
      }
      return static_cast<uint8_t>(56 + TBL_Sel8[(rank<<8) + word]);
    }


    /*!
     * @brief Return bit-pos of the 'rank'-th 1 in words[[0..]].
     * @pre 'rank'-th 1 must be found before going out of bounds.
     */
    template<typename ArrayT>
    inline uint64_t sel_0
    (
     const ArrayT & array,
     uint64_t i,
     uint64_t rank //!< 'rank' must be > 0.
     ) {
      uint64_t ret = 0;
      uint8_t cnt = cnt64(~(array[i]));
      while (rank > cnt) {
        ret += 64;
        rank -= cnt;
        cnt = cnt64(~(array[++i]));
      };
      return ret + sel64(~(array[i]), rank);
    }


    /*!
     * @brief Return bit-pos of the 'rank'-th 1 in words[[0..]].
     * @pre 'rank'-th 1 must be found before going out of bounds.
     */
    template<typename ArrayT>
    inline uint64_t sel_1
    (
     const ArrayT & array,
     uint64_t i,
     uint64_t rank //!< 'rank' must be > 0.
     ) {
      uint64_t ret = 0;
      uint8_t cnt = cnt64(array[i]);
      while (rank > cnt) {
        ret += 64;
        rank -= cnt;
        cnt = cnt64(array[++i]);
      };
      return ret + sel64(array[i], rank);
    }


    /*!
     * @brief Return # of 0's in words[[0..bitPos]].
     * @pre The bit-region must not be out of bounds.
     */
    template<typename ArrayT>
    inline uint64_t cnt_0
    (
     const ArrayT & array,
     uint64_t i,
     uint64_t bitPos //!< Bit-pos specifying the last position of the bit-region.
     ) {
      uint64_t ret = 0;
      while (bitPos > 63) {
        ret += cnt64(~(array[i++]));
        bitPos -= 64;
      }
      assert(bitPos+1 <= 64);
      ret += cnt64((~(array[i])) & bits::UINTW_MAX(static_cast<uint8_t>(bitPos+1)));
      return ret;
    }


    /*!
     * @brief Return # of 1's in words[[0..bitPos]].
     * @pre The bit-region must not be out of bounds.
     */
    template<typename ArrayT>
    inline uint64_t cnt_1
    (
     const ArrayT & array,
     uint64_t i,
     uint64_t bitPos //!< Bit-pos specifying the last position of the bit-region.
     ) {
      uint64_t ret = 0;
      while (bitPos > 63) {
        ret += cnt64(array[i++]);
        bitPos -= 64;
      }
      assert(bitPos+1 <= 64);
      ret += cnt64(array[i] & bits::UINTW_MAX(static_cast<uint8_t>(bitPos+1)));
      return ret;
    }


    /*!
     * @brief Predecessor query.
     * @ret
     *   Return the largest unset bit position smaller than or equal to a given bitPos.
     *   Return UINT64_MAX when the answer is not found after investigating 'numWords' words.
     */
    template<typename ArrayT>
    inline uint64_t pred_0
    (
     const ArrayT & array,
     const uint64_t bitPos, //!< BitPos specifying the beginning of the bit-region.
     const uint64_t numWords //!< # words to investigate in.
     ) {
      const uint64_t idx = bitPos >> 6;
      auto word = ~(array[idx]) & bits::UINTW_MAX((bitPos & 0x3f) + 1);
      uint64_t i = 0;
      while (true) {
        if (word) {
          return ((idx - i + 1) << 6) - static_cast<uint64_t>(__builtin_clzll(word)) - 1;
        }
        if (++i >= numWords) {
          break;
        }
        word = ~(array[idx - i]);
      }
      return UINT64_MAX;
    }


    /*!
     * @brief Predecessor query.
     * @ret
     *   Return the largest set bit position smaller than or equal to a given bitPos.
     *   Return UINT64_MAX when the answer is not found after investigating 'numWords' words.
     */
    template<typename ArrayT>
    inline uint64_t pred_1
    (
     const ArrayT & array,
     const uint64_t bitPos, //!< BitPos specifying the beginning of the bit-region.
     const uint64_t numWords //!< # words to investigate in.
     ) {
      const uint64_t idx = bitPos >> 6;
      auto word = (array[idx]) & bits::UINTW_MAX((bitPos & 0x3f) + 1);
      uint64_t i = 0;
      while (true) {
        if (word) {
          return ((idx - i + 1) << 6) - static_cast<uint64_t>(__builtin_clzll(word)) - 1;
        }
        if (++i >= numWords) {
          break;
        }
        word = array[idx - i];
      }
      return UINT64_MAX;
    }


    /*!
     * @brief Successor query.
     * @ret
     *   Return the smallest unset bit position larger than or equal to a given bitPos.
     *   Return UINT64_MAX when the answer is not found after investigating 'numWords' words.
     * @pre
     *   Investigating bit-region must not be out-of-bounds.
     */
    template<typename ArrayT>
    inline uint64_t succ_0
    (
     const ArrayT & array,
     const uint64_t bitPos, //!< BitPos specifying the beginning of the bit-region.
     uint64_t numWords //!< # words to investigate in.
     ) {
      const auto idx = bitPos >> 6;
      auto word = ~(array[idx]) & ~bits::UINTW_MAX((bitPos & 0x3f));
      uint64_t i = 0;
      while (true) {
        if (word) {
          return ((idx + i) << 6) + static_cast<uint64_t>(__builtin_ctzll(word));
        }
        if (++i >= numWords) {
          break;
        }
        word = ~(array[idx + i]);
      }
      return UINT64_MAX;
    }


    /*!
     * @brief Successor query.
     * @ret
     *   Return the smallest set bit position larger than or equal to a given bitPos.
     *   Return UINT64_MAX when the answer is not found after investigating 'numWords' words.
     * @pre
     *   Investigating bit-region must not be out-of-bounds.
     */
    template<typename ArrayT>
    inline uint64_t succ_1
    (
     const ArrayT & array,
     const uint64_t bitPos, //!< BitPos specifying the beginning of the bit-region.
     uint64_t numWords //!< # words to investigate in.
     ) {
      const auto idx = bitPos >> 6;
      auto word = array[idx] & ~bits::UINTW_MAX((bitPos & 0x3f));
      uint64_t i = 0;
      while (true) {
        if (word) {
          return ((idx + i) << 6) + static_cast<uint64_t>(__builtin_ctzll(word));
        }
        if (++i >= numWords) {
          break;
        }
        word = array[idx + i];
      }
      return UINT64_MAX;
    }


    /*!
     * @brief Read 'w'-bits written in the bit-region beginning at p[[bitPos..]].
     * @return Value represented by p[[bitPos..bitPos+w)).
     * @pre The bit-region must not be out of bounds.
     */
    template<typename ArrayT>
    inline uint64_t readWBits
    (
     const ArrayT & array,
     const uint64_t bitPos, //!< Bit-pos specifying the beginning position of the bit-region
     const uint8_t w, //!< Bit-width in [0, 64].
     const uint64_t mask //!< UINTW_MAX(w).
     ) {
      assert(w <= 64);

      const auto idx = (bitPos >> 6);
      const uint8_t offset = bitPos & 0x3f;
      const uint8_t loBits = 64 - offset;
      uint64_t ret = array[idx] >> offset;
      if (loBits < w) {
        ret |= array[idx + 1] << loBits;
      }
      return ret & mask;
    }


    /*!
     * @brief Simplified version of ::readWBits that can be used when reading bits in a single word.
     */
    template<typename ArrayT>
    inline uint64_t readWBits_S
    (
     const ArrayT & array, //!< Array type that returns uint64_t by [] operator
     const uint64_t bitPos, //!< Bit-pos specifying the beginning position of the bit-region
     const uint64_t mask //!< UINTW_MAX(w).
     ) {
      assert(bits::bitSize(mask) + (bitPos & 0x3f) <= 64);

      const uint8_t offset = bitPos & 0x3f;
      uint64_t ret = array[bitPos >> 6] >> offset;
      return ret & mask;
    }


    /*!
     * @brief Write 'w'-bit value 'val' to the bit-region beginning at p[[bitPos..]].
     * @pre The bit-region must not be out of bounds.
     */
    template<typename ArrayT>
    inline void writeWBits
    (
     const uint64_t val, //!< in [0, 2^w).
     ArrayT & array, //!< Array type that returns uint64_t & by [] operator
     const uint64_t bitPos, //!< Bit-pos.
     const uint8_t w, //!< Bit-width in [0, 64].
     const uint64_t mask //!< UINTW_MAX(w).
     ) {
      assert(w <= 64);
      assert(val == 0 || bits::bitSize(val) <= w);

      const auto idx = (bitPos >> 6);
      const uint8_t offset = bitPos & 0x3f;
      const uint8_t loBits = 64 - offset;
      array[idx] &= ~(mask << offset);
      array[idx] |= (val << offset);
      if (loBits < w) {
        array[idx + 1] &= ~(mask >> loBits);
        array[idx + 1] |= val >> loBits;
      }
    }


    /*!
     * @brief Simplified version of ::writeWBits that can be used when writing bits in a single word.
     */
    template<typename ArrayT>
    inline void writeWBits_S
    (
     const uint64_t val, //!< in [0, 2^w).
     ArrayT & array, //!< Array type that returns uint64_t & by [] operator
     const uint64_t bitPos, //!< Bit-pos.
     const uint64_t mask //!< UINTW_MAX(w).
     ) {
      assert(bits::bitSize(mask) + (bitPos & 0x3f) <= 64);
      assert(bits::bitSize(val) <= bits::bitSize(mask));

      const auto idx = (bitPos >> 6);
      const uint8_t offset = bitPos & 0x3f;
      array[idx] &= ~(mask << offset);
      array[idx] |= (val << offset);
    }


    template<typename SrcArrayT, typename TgtArrayT>
    inline void mvBitsLR_DiffOffs
    (
     const SrcArrayT & src,
     const uint64_t srcBitPos,
     TgtArrayT & tgt,
     const uint64_t tgtBitPos,
     uint64_t bitLen
     ) {
      uint64_t src_i = srcBitPos / 64;
      uint64_t tgt_i = tgtBitPos / 64;
      uint8_t srcOffset = srcBitPos % 64;
      uint8_t tgtOffset = tgtBitPos % 64;

      uint8_t diff1, diff2;
      uint64_t val = src[src_i] >> srcOffset;
      val <<= tgtOffset;
      val |= tgt[tgt_i] & bits::UINTW_MAX(tgtOffset);
      bitLen += tgtOffset;
      if (srcOffset > tgtOffset) {
        diff1 = srcOffset - tgtOffset;
        diff2 = 64 - diff1;
        if (bitLen <= diff2) {
          goto last;
        }
        val |= src[++src_i] << diff2;
      } else {
        diff2 = tgtOffset - srcOffset;
        diff1 = 64 - diff2;
      }
      while (bitLen > 64) {
        tgt[tgt_i++] = val;
        bitLen -= 64;
        val = src[src_i] >> diff1;
        if (bitLen <= diff2) { goto last; }
        val |= src[++src_i] << diff2;
      }
    last:
      assert(bitLen <= 64);
      const uint64_t mask = bits::UINTW_MAX(static_cast<uint8_t>(bitLen));
      tgt[tgt_i] &= ~mask;
      tgt[tgt_i] |= val & mask;
    }


    template<typename SrcArrayT, typename TgtArrayT>
    inline void mvBitsRL_DiffOffs
    (
     const SrcArrayT & src,
     const uint64_t srcBitPos,
     TgtArrayT & tgt,
     const uint64_t tgtBitPos,
     uint64_t bitLen
     ) {
      uint64_t src_i = (srcBitPos + bitLen) / 64;
      uint64_t tgt_i = (tgtBitPos + bitLen) / 64;
      uint8_t srcOffset = (srcBitPos + bitLen) % 64;
      uint8_t tgtOffset = (tgtBitPos + bitLen) % 64;

      if (srcOffset == 0) {
        srcOffset = 64;
        --src_i;
      } else if (tgtOffset == 0) {
        tgtOffset = 64;
        --tgt_i;
      }
      uint8_t diff1, diff2;
      uint64_t val = tgt[tgt_i] & ~bits::UINTW_MAX(tgtOffset);
      bitLen += 64 - tgtOffset;
      if (srcOffset > tgtOffset) {
        diff2 = srcOffset - tgtOffset;
        diff1 = 64 - diff2;
        val |= (src[src_i] & bits::UINTW_MAX(srcOffset)) >> diff2;
      } else {
        diff1 = tgtOffset - srcOffset;
        diff2 = 64 - diff1;
        val |= (src[src_i] & bits::UINTW_MAX(srcOffset)) << diff1;
        if (bitLen <= diff2) {
          goto last;
        }
        val |= src[--src_i] >> diff2;
      }
      while (bitLen > 64) {
        tgt[tgt_i--] = val;
        bitLen -= 64;
        val = src[src_i] << diff1;
        if (bitLen <= diff2) { goto last; }
        val |= (src[--src_i]) >> diff2;
      }
    last:
      assert(64 - bitLen <= 64);
      const uint64_t mask = bits::UINTW_MAX(static_cast<uint8_t>(64 - bitLen));
      tgt[tgt_i] &= mask;
      tgt[tgt_i] |= val & ~mask;
    }


    template<typename SrcArrayT, typename TgtArrayT>
    inline void mvBitsLR_SameOffs
    (
     const SrcArrayT & src,
     const uint64_t srcBitPos,
     TgtArrayT & tgt,
     const uint64_t tgtBitPos,
     uint64_t bitLen
     ) {
      assert(srcBitPos % 64 == tgtBitPos % 64);

      uint64_t src_i = srcBitPos / 64;
      uint64_t tgt_i = tgtBitPos / 64;
      const uint8_t offset = srcBitPos % 64;
      const uint64_t mask1 = bits::UINTW_MAX(offset);
      uint64_t val = tgt[tgt_i] & mask1;
      val |= src[src_i] & ~mask1;
      for (bitLen += offset; bitLen > 64; bitLen -= 64) {
        tgt[tgt_i++] = val;
        val = src[++src_i];
      }
      assert(bitLen <= 64);
      const uint64_t mask2 = bits::UINTW_MAX(static_cast<uint8_t>(bitLen));
      tgt[tgt_i] &= ~mask2;
      tgt[tgt_i] |= val & mask2;
    }


    template<typename SrcArrayT, typename TgtArrayT>
    inline void mvBitsRL_SameOffs
    (
     const SrcArrayT & src,
     const uint64_t srcBitPos,
     TgtArrayT & tgt,
     const uint64_t tgtBitPos,
     uint64_t bitLen
     ) {
      assert(srcBitPos % 64 == tgtBitPos % 64);

      uint64_t src_i = (srcBitPos + bitLen) / 64;
      uint64_t tgt_i = (tgtBitPos + bitLen) / 64;
      uint8_t offset = (srcBitPos + bitLen) % 64;
      if (offset == 0) {
        offset = 64;
        --src_i;
        --tgt_i;
      }
      const uint64_t mask1 = bits::UINTW_MAX(offset);
      uint64_t val = src[src_i] & mask1;
      val |= tgt[tgt_i] & ~mask1;
      for (bitLen += 64 - offset; bitLen > 64; bitLen -= 64) {
        tgt[tgt_i--] = val;
        val = src[--src_i];
      }
      assert(64 - bitLen <= 64);
      const uint64_t mask2 = bits::UINTW_MAX(static_cast<uint8_t>(64 - bitLen));
      tgt[tgt_i] &= mask2;
      tgt[tgt_i] |= val & ~mask2;
    }


    /*!
     * @brief Move bits (by scanning bits from left to right) from src-region to tgt-region (see also ::mvBits).
     * @attention
     *   This function should be used when bits are moved backwardly, i.e., "position of src-region" >= "position of tgt-region", or non-overlapping.
     */
    template<typename SrcArrayT, typename TgtArrayT>
    inline void mvBitsLR
    (
     const SrcArrayT & src,
     const uint64_t srcBitPos,
     TgtArrayT & tgt,
     const uint64_t tgtBitPos,
     uint64_t bitLen
     ) {
      if (srcBitPos % 64 != tgtBitPos % 64) {
        mvBitsLR_DiffOffs(src, srcBitPos, tgt, tgtBitPos, bitLen);
      } else {
        mvBitsLR_SameOffs(src, srcBitPos, tgt, tgtBitPos, bitLen);
      }
    }


    /*!
     * @brief Move bits (by scanning from right to left) from src-region to tgt-region (see also ::mvBits).
     * @attention This function should be used when bits are moved forwardly, i.e., "position of src-region" < "position of tgt-region".
     * @attention Src-region is designated by the bit-region *ending (excluded)* at 'srcOffset' bit of 'src' and of length 'bitLen'.
     *            Tgt-region is similarly defined.
     */
    template<typename SrcArrayT, typename TgtArrayT>
    inline void mvBitsRL
    (
     const SrcArrayT & src,
     const uint64_t srcBitPos,
     TgtArrayT & tgt,
     const uint64_t tgtBitPos,
     uint64_t bitLen
     ) {
      if (srcBitPos % 64 != tgtBitPos % 64) {
        mvBitsRL_DiffOffs(src, srcBitPos, tgt, tgtBitPos, bitLen);
      } else {
        mvBitsRL_SameOffs(src, srcBitPos, tgt, tgtBitPos, bitLen);
      }
    }


    /*!
     * @brief Move bits from src-region to tgt-region.
     * @attention When src-region and tgt-region overlap, the overlapped part of src-region is overwritten.
     *            The other part of src-region is not changed.
     * @note
     *   The process will be forwarded to one of the four functions based on the following two (orthogonal) criteria.
     *   - move bits forwardly or backwardly.
     *   - offsets are same or not.
     */
    template<typename SrcArrayT, typename TgtArrayT>
    inline void mvBits
    (
     const SrcArrayT & src,
     const uint64_t srcBitPos,
     TgtArrayT & tgt,
     const uint64_t tgtBitPos,
     uint64_t bitLen
     ) {
      if (srcBitPos >= tgtBitPos) {
        mvBitsLR(src, srcBitPos, tgt, tgtBitPos, bitLen);
      } else {
        mvBitsRL(src, srcBitPos, tgt, tgtBitPos, bitLen);
      }
    }


    template<typename SrcArrayT, typename TgtArrayT>
    inline void cpBits_DiffOffs
    (
     const SrcArrayT & src,
     const uint64_t srcBitPos,
     TgtArrayT & tgt,
     const uint64_t tgtBitPos,
     uint64_t bitLen
     ) {
      mvBitsLR_DiffOffs(src, srcBitPos, tgt, tgtBitPos, bitLen);
    }


    template<typename SrcArrayT, typename TgtArrayT>
    inline void cpBits_SameOffs
    (
     const SrcArrayT & src,
     const uint64_t srcBitPos,
     TgtArrayT & tgt,
     const uint64_t tgtBitPos,
     uint64_t bitLen
     ) {
      assert(srcBitPos % 64 == tgtBitPos % 64);

      mvBitsLR_SameOffs(src, srcBitPos, tgt, tgtBitPos, bitLen);
    }


    /*!
     * @brief Copy bits from src-region to tgt-region (redirect to ::mvBitsLR). 
     * @attention
     *   This function should be used when bits are moved backwardly, i.e., "position of src-region" >= "position of tgt-region", or non-overlapping.
     */
    template<typename SrcArrayT, typename TgtArrayT>
    inline void cpBits
    (
     const SrcArrayT & src,
     const uint64_t srcBitPos,
     TgtArrayT & tgt,
     const uint64_t tgtBitPos,
     uint64_t bitLen
     ) {
      mvBitsLR(src, srcBitPos, tgt, tgtBitPos, bitLen);
    }


    inline void mvBytes
    (
     const void * src,
     void * tgt,
     const uint64_t bytes
     ) {
      std::memmove(tgt, src, bytes); // tgt and src are in this order in memmove
    }


    inline void cpBytes
    (
     const void * src,
     void * tgt,
     const uint64_t bytes
     ) {
      std::memcpy(tgt, src, bytes); // tgt and src are in this order in memcpy
    }

  } // namespace bits
} // namespace itmmti


#endif
