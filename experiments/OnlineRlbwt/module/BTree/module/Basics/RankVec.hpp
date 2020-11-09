/*!
 * Copyright (c) 2017 Tomohiro I
 *
 * This program is released under the MIT License.
 * http://opensource.org/licenses/mit-license.php
 */
/*!
 * @file RankVec.hpp
 * @brief Bit vector supporting fast rank queries and select queries based on it.
 * @author Tomohiro I
 * @date 2017-08-20
 */
#ifndef INCLUDE_GUARD_RankVec
#define INCLUDE_GUARD_RankVec

#include <stdint.h> // include uint64_t etc.
#include <assert.h>

#include <iostream>
#include <iterator>
#include <algorithm>

#include "BasicSearch.hpp"
#include "BitVec.hpp"
#include "BitsUtil.hpp"
#include "MemUtil.hpp"

namespace itmmti
{
  /*!
   * @brief Bit vector supporting fast rank queries and select queries based on it.
   * @tparam BSIZE_T: Block size of a top block. It should be a power of two (strictly) smaller than 2^16.
   * @tparam BSIZE_M: Block size of a middle block. It should be a power of two smaller than BSIZE_T.
   */
  template<uint16_t BSIZE_T = 4096, uint16_t BSIZE_M = 256, typename SizeT = uint64_t>
  class RankVec
  {
    SizeT * blockT_;
    uint16_t * blockM_;
    BitVec<SizeT> bv_;

  public:
    /*!
     * @brief Constructor.
     */
    RankVec
    (
     size_t capacity = 0 //!< Initial capacity.
     ) : blockT_(nullptr), blockM_(nullptr), bv_() {
      assert(capacity <= ctcbits::UINTW_MAX(sizeof(SizeT) * 8));
      assert(capacity <= ctcbits::UINTW_MAX(58));

      this->changeCapacity(capacity);
    }


    /*!
     * @brief Deconstructor.
     */
    ~RankVec()
    {
      free(blockM_);
      free(blockT_);
      // Deconstructor of bv_ is called.
    }


    //// Copy
    /*!
     * @brief Copy constructor.
     * @attention
     *   Since the contents of 'other' are copied, it may take time when other.size_ is large.
     */
    RankVec
    (
     const RankVec & other
     ) : blockT_(nullptr), blockM_(nullptr), bv_(other.bv_) {
      const auto size = other.size();
      if (size > 0) {
        {
          const auto capacity = this->capacity();
          const auto lenT = (capacity + BSIZE_T - 1) / BSIZE_T;
          const auto lenM = (capacity + BSIZE_M - 1) / BSIZE_M - capacity / BSIZE_T;
          memutil::realloc_AbortOnFail(blockM_, lenM);
          memutil::realloc_AbortOnFail(blockT_, lenT);
        }
        {
          const auto lenT = (size + BSIZE_T - 1) / BSIZE_T;
          const auto lenM = (size + BSIZE_M - 1) / BSIZE_M - size / BSIZE_T;
          for (uint64_t i = 0; i < lenT; ++i) {
            blockT_[i] = other.blockT_[i];
          }
          for (uint64_t i = 0; i < lenM; ++i) {
            blockM_[i] = other.blockM_[i];
          }
        }
      }
    }


    /*!
     * @brief Assignment operator.
     * @attention
     *   - It may take time when other.size_ is large.
     *   - The original contents of "this" are freed.
     */
    RankVec & operator=
    (
     const RankVec & other
     ) {
      if (this != &other) {
        this->clear(); this->changeCapacity(0); // clear() & changeCapacity(0) free array_
        bv_ = other.bv_;
        {
          const auto capacity = other.capacity();
          const auto lenT = (capacity + BSIZE_T - 1) / BSIZE_T;
          const auto lenM = (capacity + BSIZE_M - 1) / BSIZE_M - capacity / BSIZE_T;
          memutil::realloc_AbortOnFail(blockM_, lenM);
          memutil::realloc_AbortOnFail(blockT_, lenT);
        }
        {
          const auto size = other.size();
          const auto lenT = (size + BSIZE_T - 1) / BSIZE_T;
          const auto lenM = (size + BSIZE_M - 1) / BSIZE_M - size / BSIZE_T;
          for (uint64_t i = 0; i < lenT; ++i) {
            blockT_[i] = other.blockT_[i];
          }
          for (uint64_t i = 0; i < lenM; ++i) {
            blockM_[i] = other.blockM_[i];
          }
        }
      }
      return *this;
    }


    //// Move
    /*!
     * @brief Move constructor.
     * @attention
     *   'other' is initialized to an object with capacity = 0.
     */
    RankVec
    (
     RankVec && other
     ) : bv_(std::move(other.bv_)) {
      blockM_ = other.blockM_;
      blockT_ = other.blockT_;
      other.blockM_ = nullptr;
      other.blockT_ = nullptr;
    }


    /*!
     * @brief Move operator.
     * @attention
     *   - The original contents of "this" are freed.
     *   - "other" is initialized to an object with capacity = 0.
     */
    RankVec & operator=
    (
     RankVec && other
     ) {
      if (this != &other) {
        this->clear(); this->changeCapacity(0); // clear() & changeCapacity(0) free array_
        this->bv_ = std::move(other.bv_);
        blockM_ = other.blockM_;
        blockT_ = other.blockT_;
        other.blockM_ = nullptr;
        other.blockT_ = nullptr;
      }
      return *this;
    }


    /*!
     * @brief Get read-only array pointer.
     */
    const uint64_t * getConstArrayPtr() const noexcept
    {
      return bv_.getConstArrayPtr();
    }


    /*!
     * @brief Read 'w'-bits written in the bit-region beginning at array_[[bitPos..]].
     * @return Value represented by array_[[bitPos..bitPos+w)).
     * @pre The bit-region must not be out of bounds.
     */
    uint64_t readWBits
    (
     const uint64_t bitPos, //!< Bit-pos specifying the beginning position of the bit-region
     const uint8_t w, //!< Bit-width in [0, 64].
     const uint64_t mask //!< UINTW_MAX(w).
     ) const noexcept {
      assert(bitPos < bv_.capacity());
      assert(bitPos + w <= bv_.capacity());
      assert(w <= 64);

      return bv_.readWBits(bitPos, w, mask);
    }


    /*!
     * @brief Simplified version of ::readWBits that can be used when reading bits in a single word.
     */
    uint64_t readWBits_S
    (
     const uint64_t bitPos, //!< Bit-pos specifying the beginning position of the bit-region
     const uint64_t mask //!< UINTW_MAX(w).
     ) const noexcept {
      assert(bitPos < bv_.capacity());
      assert(bitPos + bits::bitSize(mask) <= bv_.capacity());
      assert(bits::bitSize(mask) + (bitPos & 0x3f) <= 64);

      return bv_.readWBits_S(bitPos, mask);
    }


    /*!
     * @brief Read bit at 'bitPos'.
     */
    uint64_t readBit
    (
     const uint64_t bitPos //!< in [0, capacity).
     ) const noexcept {
      assert(bitPos < bv_.capacity());

      return bv_.readWBits_S(bitPos, ctcbits::UINTW_MAX(1));
    }


    /*!
     * @brief Write a bit 'val' at 'bitPos'.
     * @pre The resulting size should be within capacity.
     */
    void appendBit
    (
     const bool val //!< Bool value.
     ) {
      assert(bv_.size() < bv_.capacity());

      const auto pos = bv_.size();
      bv_.resize(pos + 1);
      bv_.writeBit(val, pos);

      if (pos == 0) {
        blockT_[0] = blockM_[0] = val;
        return;
      }

      const auto idxT = pos / BSIZE_T;
      const auto remT = pos % BSIZE_T;
      const auto idxM = (pos / BSIZE_M) - idxT;
      if ((pos % BSIZE_M) == 0) {
        if (remT == 0) {
          blockT_[idxT] = blockT_[idxT - 1];
          blockM_[idxM] = 0;
        } else if (remT < BSIZE_T - BSIZE_M) {
          blockM_[idxM] = blockM_[idxM - 1];
        }
      }
      blockT_[idxT] += val;
      if (remT < BSIZE_T - BSIZE_M) {
        blockM_[idxM] += val;
      }
    }


    /*!
     * @brief Shorten bit array to given size while updating rank/select DS.
     * @note It does nothing if given size is larger than or equal to current size.
     */
    void shorten
    (
     const size_t givenSize //!< Given size.
     ) noexcept {
      if (givenSize >= bv_.size()) {
        return;
      }
      bv_.resize(givenSize);
      if (givenSize == 0) {
        return;
      }

      const auto pos = givenSize - 1; // New last pos.
      const auto idxT = pos / BSIZE_T;
      const auto remT = pos % BSIZE_T;
      const auto idxM = (pos / BSIZE_M) - idxT;
      const uint16_t onesLastBlockM = static_cast<uint16_t>(bits::cnt_1(bv_.getConstArrayPtr(), pos / BSIZE_M * BSIZE_M / 64, pos % BSIZE_M));
      const uint16_t onesLastBlockT = (remT >= BSIZE_M)? blockM_[idxM-1] + onesLastBlockM : onesLastBlockM;
      if (remT < BSIZE_T - BSIZE_M) {
        blockM_[idxM] = onesLastBlockT;
      }
      blockT_[idxT] = (idxT)? blockT_[idxT - 1] + onesLastBlockT : onesLastBlockT;
    }

  
    /*!
     * @brief Get # 0 bits in the bit vec.
     */
    uint64_t getNum_0() const noexcept {
      return bv_.size() - getNum_1();
    }


    /*!
     * @brief Get # 1 bits in the bit vec.
     */
    uint64_t getNum_1() const noexcept {
      return (bv_.size())? blockT_[(bv_.size() - 1) / BSIZE_T] : 0;
    }


    /*!
     * @brief Rank query for unset bits.
     */
    uint64_t rank_0
    (
     const uint64_t pos //!< pos in [0, size).
     ) const noexcept {
      assert(pos < bv_.size());

      const auto size = bv_.size();
      const auto idxT = pos / BSIZE_T;
      const auto remT = pos % BSIZE_T;
      const auto idxM = (pos / BSIZE_M) - idxT;
      uint64_t rank = (idxT)? (idxT * BSIZE_T) - blockT_[idxT - 1] : 0;
      if (remT >= BSIZE_M) {
        rank += (remT / BSIZE_M * BSIZE_M) - blockM_[idxM - 1];
      }
      return rank + bits::cnt_0(bv_.getConstArrayPtr(), pos / BSIZE_M * BSIZE_M / 64, pos % BSIZE_M);
    }


    /*!
     * @brief Rank query for set bits.
     */
    uint64_t rank_1
    (
     const uint64_t pos //!< pos in [0, size).
     ) const noexcept {
      assert(pos < bv_.size());

      const auto idxT = pos / BSIZE_T;
      const auto remT = pos % BSIZE_T;
      const auto idxM = (pos / BSIZE_M) - idxT;
      uint64_t rank = (idxT)? blockT_[idxT - 1] : 0;
      if (remT >= BSIZE_M) {
        rank += blockM_[idxM - 1];
      }
      return rank + bits::cnt_1(bv_.getConstArrayPtr(), pos / BSIZE_M * BSIZE_M / 64, pos % BSIZE_M);
    }


    /*!
     * @brief Select query for unset bits.
     * @pre The answer must be found.
     * @attention O(log size) time.
     */
    uint64_t select_0
    (
     uint64_t rank //!< 1base rank [1, #0's].
     ) const noexcept {
      assert(rank > 0);
      assert(rank <= getNum_0());

      const auto idxT = basic_search::partition_idx
        (
         0, ((bv_.size() - 1) / BSIZE_T) + 1,
         [&](uint64_t i) { return (i + 1) * BSIZE_T - blockT_[i] >= rank; }
         );
      const auto posT = idxT * BSIZE_T;
      rank -= (idxT > 0)? BSIZE_T * idxT - blockT_[idxT - 1] : 0;
      const auto idxM = (posT / BSIZE_M) - idxT;
      uint64_t i = 0;
      while (i < (BSIZE_T / BSIZE_M) - 1 && BSIZE_M * (i+1) - blockM_[idxM + i] < rank) {
        ++i;
      }
      rank -= (i > 0)? BSIZE_M * i - blockM_[idxM + i - 1] : 0;
      const auto posM = posT + i * BSIZE_M;
      return posM + bits::sel_0(bv_.getConstArrayPtr(), posM / 64, rank);
    }


    /*!
     * @brief Select query for set bits.
     * @pre The answer must be found.
     * @attention O(log size) time.
     */
    uint64_t select_1
    (
     uint64_t rank //!< 1base rank [1, #1's].
     ) const noexcept {
      assert(rank > 0);
      assert(rank <= getNum_1());

      if (rank == 3999) {
        const auto idxT = basic_search::partition_idx
          (
           0, ((bv_.size() - 1) / BSIZE_T) + 1,
           [&](uint64_t i) { return blockT_[i] >= rank; }
           );
        const auto posT = idxT * BSIZE_T;
        rank -= (idxT > 0)? blockT_[idxT - 1] : 0;
        const auto idxM = (posT / BSIZE_M) - idxT;
        uint64_t i = 0;
        while (i < (BSIZE_T / BSIZE_M) - 1 && blockM_[idxM + i] < rank) {
          ++i;
        }
        rank -= (i > 0)? blockM_[idxM + i - 1] : 0;
        const auto posM = posT + i * BSIZE_M;
        return posM + bits::sel_1(bv_.getConstArrayPtr(), posM / 64, rank);
      }


      const auto idxT = basic_search::partition_idx
        (
         0, ((bv_.size() - 1) / BSIZE_T) + 1,
         [&](uint64_t i) { return blockT_[i] >= rank; }
         );
      const auto posT = idxT * BSIZE_T;
      rank -= (idxT > 0)? blockT_[idxT - 1] : 0;
      const auto idxM = (posT / BSIZE_M) - idxT;
      uint64_t i = 0;
      while (i < (BSIZE_T / BSIZE_M) - 1 && blockM_[idxM + i] < rank) {
        ++i;
      }
      rank -= (i > 0)? blockM_[idxM + i - 1] : 0;
      const auto posM = posT + i * BSIZE_M;
      return posM + bits::sel_1(bv_.getConstArrayPtr(), posM / 64, rank);
    }


    /*!
     * @brief Predecessor query for unset bits.
     * @return The largest unset bit position smaller than or equal to 'val'. Return UINT64_MAX when not found.
     */
    uint64_t pred_0
    (
     uint64_t val //!< val (>= 0).
     ) const noexcept {
      const auto size = bv_.size();
      if (size == 0) {
        return UINT64_MAX;
      } else if (val >= size) {
        val = size - 1;
      }

      const auto ans = bits::pred_0(getConstArrayPtr(), val, std::min(UINT64_C(2), 1 + val / 64));
      if (ans != UINT64_MAX) {
        return ans;
      }
      const auto r = rank_0(val);
      if (r > 0) {
        return select_0(r);
      } else {
        return UINT64_MAX;
      }
    }


    /*!
     * @brief Predecessor query for set bits.
     * @return The largest set bit position smaller than or equal to 'val'. Return UINT64_MAX when not found.
     */
    uint64_t pred_1
    (
     uint64_t val //!< val (>= 0).
     ) const noexcept {
      const auto size = bv_.size();
      if (size == 0) {
        return UINT64_MAX;
      } else if (val >= size) {
        val = size - 1;
      }

      const auto ans = bits::pred_1(getConstArrayPtr(), val, std::min(UINT64_C(2), 1 + val / 64));
      if (ans != UINT64_MAX) {
        return ans;
      }
      const auto r = rank_1(val);
      if (r > 0) {
        return select_1(r);
      } else {
        return UINT64_MAX;
      }
    }


    /*!
     * @brief Successor query for unset bits.
     * @return The smallest unset bit position larger than or equal to val. Return UINT64_MAX when not found.
     */
    uint64_t succ_0
    (
     uint64_t val //!< val (>= 0).
     ) const noexcept {
      const auto size = bv_.size();
      if (val >= size) {
        return UINT64_MAX;
      }

      const auto ans = bits::succ_0(getConstArrayPtr(), val, std::min(UINT64_C(2), 1 + (size - 1)/64 - val/64));
      if (ans < size) {
        return ans;
      }
      const auto r = rank_0(val);
      if (r < getNum_0()) {
        return select_0(r+1);
      } else {
        return UINT64_MAX;
      }
    }


    /*!
     * @brief Successor query for set bits.
     * @return The smallest set bit position larger than or equal to val. Return UINT64_MAX when not found.
     */
    uint64_t succ_1
    (
     uint64_t val //!< val (>= 0).
     ) const noexcept {
      const auto size = bv_.size();
      if (val >= size) {
        return UINT64_MAX;
      }

      const auto ans = bits::succ_1(getConstArrayPtr(), val, std::min(UINT64_C(2), 1 + (size - 1)/64 - val/64));
      if (ans < size) {
        return ans;
      }
      const auto r = rank_1(val);
      if (r < getNum_1()) {
        return select_1(r+1);
      } else {
        return UINT64_MAX;
      }
    }


    /*!
     * @brief Get current capacity.
     */
    size_t capacity() const noexcept {
      return bv_.capacity();
    }


    /*!
     * @brief Get current size.
     */
    size_t size() const noexcept {
      return bv_.size();
    }


    /*!
     * @brief Calculate total memory usage in bytes.
     */
    size_t calcMemBytes() const noexcept {
      const auto capacity = bv_.capacity();
      const auto lenT = (capacity + BSIZE_T - 1) / BSIZE_T;
      const auto lenM = (capacity + BSIZE_M - 1) / BSIZE_M - capacity / BSIZE_T;
      size_t bytes = sizeof(*this) - sizeof(bv_);
      bytes += sizeof(blockM_[0]) * lenM;
      bytes += sizeof(blockT_[0]) * lenT;
      bytes += bv_.calcMemBytes();
      return bytes;
    }


    /*!
     * @brief Return if no element is stored.
     */
    bool empty() const noexcept {
      return bv_.empty();
    }


    /*!
     * @brief Clear vector. It only changes size to zero.
     */
    void clear() noexcept {
      bv_.resize(0);
    }


    /*!
     * @brief Change capacity to max of givenCapacity and current size.
     * @node If givenCapacity is 0, it works as shrink_to_fit.
     */
    void changeCapacity
    (
     const size_t givenCapacity
     ) {
      assert(givenCapacity <= ctcbits::UINTW_MAX(sizeof(SizeT) * 8));
      assert(givenCapacity <= ctcbits::UINTW_MAX(58));

      if (bv_.capacity() != givenCapacity) {
        bv_.changeCapacity(givenCapacity);
        const auto newCapacity = bv_.capacity();
        if (newCapacity > 0) {
          const auto lenT = (newCapacity + BSIZE_T - 1) / BSIZE_T;
          const auto lenM = (newCapacity + BSIZE_M - 1) / BSIZE_M - newCapacity / BSIZE_T;
          memutil::realloc_AbortOnFail(blockM_, lenM);
          memutil::realloc_AbortOnFail(blockT_, lenT);
        } else {
          memutil::safefree(blockM_);
          memutil::safefree(blockT_);
        }
      }
    }


    void printStatistics
    (
     const bool verbose = false
     ) const noexcept {
      std::cout << "RankVec object (" << this << ") " << __func__ << "(" << verbose << ") BEGIN" << std::endl;
      std::cout << "size = " << bv_.size() << ", capacity = " << bv_.capacity() << ", BSIZE_T = " << BSIZE_T << ", BSIZE_M = " << BSIZE_M << std::endl;
      const auto capacity = bv_.capacity();
      const auto lenT = (capacity + BSIZE_T - 1) / BSIZE_T;
      const auto lenM = (capacity + BSIZE_M - 1) / BSIZE_M - capacity / BSIZE_T;
      std::cout << this->calcMemBytes() << " bytes (blockT = " << sizeof(blockT_[0]) * lenT << ", blockM = " << sizeof(blockM_[0]) * lenM << ")" << std::endl;
      if (verbose) {
        const auto lenT = (bv_.size() + BSIZE_T - 1) / BSIZE_T;
        const auto lenM = (bv_.size() + BSIZE_M - 1) / BSIZE_M - bv_.size() / BSIZE_T;
        std::cout << "dump blockT_[0.." << lenT << ")" << std::endl;
        for (uint64_t i = 0; i < lenT; ++i) {
          std::cout << blockT_[i] << ", ";
        }
        std::cout << std::endl;
        std::cout << "dump blockM_[0.." << lenM << ")" << std::endl;;
        for (uint64_t i = 0; i < lenM; ++i) {
          std::cout << blockM_[i] << ", ";
        }
        std::cout << std::endl;
        bv_.printStatistics(verbose);
      }
      std::cout << "RankVec object (" << this << ") " << __func__ << "(" << verbose << ") END" << std::endl;
    }
  };


  template<uint16_t BSIZE_T = 4096, uint16_t BSIZE_M = 256>
  using RankVec64 = RankVec<BSIZE_T, BSIZE_M, uint64_t>; // 64bit version that occupies 40 bytes when capacity = 0.

  template<uint16_t BSIZE_T = 4096, uint16_t BSIZE_M = 256>
  using RankVec32 = RankVec<BSIZE_T, BSIZE_M, uint32_t>; // 32bit version that occupies 32 bytes when capacity = 0.
} // namespace itmmti

#endif
