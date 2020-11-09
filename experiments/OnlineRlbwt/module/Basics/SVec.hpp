/*!
 * Copyright (c) 2017 Tomohiro I
 *
 * This program is released under the MIT License.
 * http://opensource.org/licenses/mit-license.php
 */
/*!
 * @file SVec.hpp
 * @brief Sarray implementation supporting semi-dynamic update.
 * @author Tomohiro I
 * @date 2017-08-22
 */
#ifndef INCLUDE_GUARD_SVec
#define INCLUDE_GUARD_SVec

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
   * @brief Implimentation of Sarray [1] supporting semi-dynamic update.
   * @tparam RSDicT: Data structure supporting fast rank/select queries.
   * @detail
   *   Vector stores a sequence of increasing uints in 1.92 m + m lg (n/m) + o(m) bits,
   *   where n is the max uint and m is the size of sequence.
   *   Each uint is partitioned into upper bits and lower bits (say the bit-width of lower bits is 'loW'):
   *   The higher part is stored by gapped-encoding in Rank/Select data structure, and
   *   the lower part is stored in 'loW'-bits packed array.
   *   To obtain the optimal space saving, 'n >> loW' should be around 1.44 m.
   * @par
   *   The underlying idea can date back to [2, 3]. So, the data structure is also called as 
   *   the Elias-Fano coding, Elias-Fano indexable dictionary, Elias-Fano representation for monotone sequence, etc.
   * @par References
   *   [1] Daisuke Okanohara and Kunihiko Sadakane. Practical Entropy-Compressed Rank/Select Dictionary, ALENEX, 2007
   *   [2] Peter Elias. Efficient storage and retrieval by content and address of static files, Journal of the ACM, 1974
   *   [3] Robert Fano. On the number of bits required to implement an associative memory. Memorandum 61. Computer Structures Group, Project MAC, MIT, 1971
   */
  template<class RSDicT>
  class SVec
  {
    RSDicT rsv_;
    WBitsVec wbv_;

  public:
    static uint8_t calcOptimalLoW
    (
     size_t max, //!< Max uint supposed to be stored.
     size_t size //!< #uints.
     ) noexcept {
      assert(size > 0);

      return bits::bitSize(static_cast<uint64_t>(max / (size * 1.44)));
    }


  public:
    SVec
    (
     uint8_t loW = 8, //!< Initial bit-width for lower part of uint stored in 'wbv_'.
     size_t capacity = 0 //!< Initial capacity.
     ) : rsv_(capacity), wbv_(loW, capacity) {
      assert(capacity <= ctcbits::UINTW_MAX(58));
      assert(0 < loW && loW <= 64);
    }


    /*!
     * @brief Deconstructor.
     */
    ~SVec()
    {
      // Deconstructors of rsv_ and wbv_ are called.
    }


    //// Copy
    /*!
     * @brief Copy constructor.
     * @attention
     *   Since the contents of "other" are copied, it may take time when other.size_ is large.
     */
    SVec
    (
     const SVec & other
     ) : rsv_(other.rsv_), wbv_(other.wbv_) {
    }


    /*!
     * @brief Assignment operator.
     * @attention
     *   - It may take time when other.size_ is large.
     *   - The original contents of "this" are freed.
     */
    SVec & operator=
    (
     const SVec & other
     ) {
      if (this != &other) {
        rsv_ = other.rsv_;
        wbv_ = other.wbv_;
      }
      return *this;
    }


    //// Move
    /*!
     * @brief Move constructor.
     * @attention
     *   "other" is initialized to an object with capacity = 0.
     */
    SVec
    (
     SVec && other
     ) : rsv_(std::move(other.rsv_)), wbv_(std::move(other.wbv_)) {
    }


    /*!
     * @brief Move operator.
     * @attention
     *   - The original contents of "this" are freed.
     *   - "other" is initialized to an object with capacity = 0.
     */
    SVec & operator=
    (
     SVec && other
     ) {
      if (this != &other) {
        rsv_ = std::move(other.rsv_);
        wbv_ = std::move(other.wbv_);
      }
      return *this;
    }


    /*!
     * @brief Get the current bit-width for least significant bits.
     */
    uint64_t getLoW() const noexcept {
      return wbv_.getW();
    }


    /*!
     * @brief Get #unset-bits in the bit vec.
     */
    uint64_t getNum_0() const noexcept {
      auto num1 = getNum_1();
      if (num1 > 0) {
        return this->getMax() - num1;
      } else {
        return 0;
      }
    }


    /*!
     * @brief Get #set-bits in the bit vec.
     */
    uint64_t getNum_1() const noexcept {
      return wbv_.size();
    }


    /*!
     * @brief Get max value (the last position of set-bit).
     */
    uint64_t getMax() const noexcept {
      assert(wbv_.size());

      return (rsv_.getNum_0() << wbv_.getW()) + wbv_.read(wbv_.size() - 1);
    }


    /*!
     * @brief Predecessor query for set-bits.
     * @return The largest set-bit position smaller than or equal to val. Return UINT64_MAX when not found.
     */
    uint64_t pred_1
    (
     uint64_t val //!< val (>= 0).
     ) const noexcept {
      const auto size = wbv_.size();
      if (size == 0) {
        return UINT64_MAX;
      } else if (val >= getMax()) {
        return size;
      }

      const auto r = rank_1(val);
      if (r) {
        return select_1(r);
      } else {
        return UINT64_MAX;
      }
    }


    /*!
     * @brief Successor query for set-bits.
     * @return The smallest set-bit position smaller than or equal to val. Return UINT64_MAX when not found.
     */
    uint64_t succ_1
    (
     uint64_t val //!< val (>= 0).
     ) const noexcept {
      const auto size = wbv_.size();
      if (size == 0 || val > getMax()) {
        return UINT64_MAX;
      }

      const auto r = rank_1(val);
      const auto s = select_1(r);
      if (s < val) {
        return select_1(r + 1);
      } else {
        return s;
      }
    }


    /*!
     * @brief Rank query.
     */
    uint64_t rank_1
    (
     const uint64_t pos //!< pos in [0, UINT64_MAX].
     ) const noexcept {
      const auto size = wbv_.size();
      if (size == 0) {
        return 0;
      }

      const auto loW = wbv_.getW();
      const auto hiBits = pos >> loW;
      const auto hiMax = rsv_.getNum_0();
      if (hiBits > hiMax) {
        return size;
      }
      const auto rvPos = (hiBits)? rsv_.select_0(hiBits) + 1 : 0;
      auto rank_lb = rsv_.rank_1(rvPos);
      if (rsv_.readBit(rvPos) == 0) {
        return rank_lb;
      }
      // Now, rank_lb >= 1.
      const auto rank_ub = (hiBits < hiMax)? rank_lb + rsv_.succ_0(rvPos) - rvPos : size + 1;
      const auto key = pos & bits::UINTW_MAX(loW);
      if (key < wbv_.read(rank_ub - 2)) {
        return basic_search::partition_idx(rank_lb - 1, rank_ub - 1, [&](uint64_t i) { return key < wbv_[i]; } );
      } else {
        return rank_ub - 1;
      }
    }


    /*!
     * @brief Select query.
     * @pre The answer must be found.
     * @note
     *   O(t) time, where t is the time for select_1 (depends on implementation of rsv_).
     */
    uint64_t select_1
    (
     uint64_t rank //!< 1base rank [1, #set-bits].
     ) const noexcept {
      assert(rank > 0);
      assert(rank <= getNum_1());

      return (rsv_.rank_0(rsv_.select_1(rank)) << wbv_.getW()) + wbv_.read(rank - 1);
    }


    /*!
     * @brief Rank query for unset-bits.
     */
    uint64_t rank_0
    (
     const uint64_t pos //!< pos in [0, UINT64_MAX].
     ) const noexcept {
      return pos + 1 - rank_1(pos);
    }


    /*!
     * @brief Select query for unset-bits.
     * @pre The answer must be found.
     * @note
     *   O(t lg size) time, where t is the time for select_1 (depends on implementation of rsv_).
     */
    uint64_t select_0
    (
     uint64_t rank //!< 1base rank [1, #unset-bits].
     ) const noexcept {
      assert(rank > 0);
      assert(rank <= getNum_0());

      if (rank < (rsv_.succ_1(0) << wbv_.getW()) + wbv_.read(0)) {
        return rank - 1;
      }
      const auto idx = basic_search::partition_idx
        (
         0, wbv_.size(),
         [&](uint64_t i) { return rank <= select_1(i+1) - i; }
         );

      return rank + idx - 1;
    }


    /*!
     * @brief Get current capacity.
     */
    size_t capacity() const noexcept {
      return wbv_.capacity();
    }


    /*!
     * @brief Get current size.
     */
    size_t size() const noexcept {
      return wbv_.size();
    }


    /*!
     * @brief Calculate total memory usage in bytes.
     */
    size_t calcMemBytes() const noexcept {
      size_t bytes = sizeof(*this);
      bytes += wbv_.calcMemBytes();
      bytes += rsv_.calcMemBytes();
      return bytes;
    }


    /*!
     * @brief Return if no element is stored.
     */
    bool empty() const noexcept {
      return wbv_.empty();
    }


    /*!
     * @brief Clear vector. It only changes size to zero.
     */
    void clear() noexcept {
      wbv_.clear(0);
      rsv_.clear(0);
    }


    /*!
     * @brief Change capacity to max of givenCapacity and current size.
     */
    void changeCapacity
    (
     const size_t givenCapacity
     ) {
      assert(givenCapacity <= ctcbits::UINTW_MAX(58));

      wbv_.changeCapacity(givenCapacity);
    }


    /*!
     * @brief Write uint "val" at the end.
     * @pre The resulting size should be within capacity.
     */
    void append
    (
     const uint64_t val, //!< Uint val.
     double marginFactor = 1.5 //!< For RankVec, "new_rsv_size * marginFactor" space is reserved (default = 1.5).
     ) {
      assert(wbv_.size() < wbv_.capacity());
      assert(size() == 0 || getMax() < val);

      const auto pos = wbv_.size();
      const auto loW = wbv_.getW();
      wbv_.resize(pos + 1);
      wbv_.write(val & bits::UINTW_MAX(loW), pos);

      const auto diff0 = (val >> loW) - rsv_.getNum_0();
      const auto rvSizeNew = rsv_.size() + diff0 + 1;
      if (rvSizeNew > rsv_.capacity()) {
        rsv_.changeCapacity(static_cast<size_t>(rvSizeNew * marginFactor));
      }
      for (uint64_t i = 0; i < diff0; ++i) {
        rsv_.appendBit(false);
      }
      rsv_.appendBit(true);
    }


    /*!
     * @brief Rebalance bit-widths of hi/lo bits (give new "loW').
     * @note
     *   - It takes linear time.
     *   - Give optional arguments (minCapacity, marginFactor, doShrink) to control how reallocation is done when it is required.
     */
    void convert
    (
     const uint8_t loW, //!< New bit-width.
     size_t minCapacity = 0, //!< Minimum support for capacity (default = 0, which turns into 'wbv_.size_').
     double marginFactor = 1.0, //!< For RankVec, 'new_rsv_size * marginFactor' space is reserved (default = 1.0).
     bool doShrink = false //!< If true, 'wbv_' is reallocated to fit max(size_, minCapacity).
     ) {
      assert(0 < loW && loW <= 64);
      assert(minCapacity <= ctcbits::UINTW_MAX(58));
      assert(marginFactor >= 1.0);

      const auto loW_old = wbv_.getW();
      if (loW == loW_old) {
        if (doShrink) {
          { // Update lower part.
            wbv_.convert(loW, minCapacity, doShrink);
          }
          { // Update higher part.
            const uint64_t sizeWithMargin = static_cast<uint64_t>(rsv_.size() * marginFactor);
            if (sizeWithMargin < rsv_.capacity()) {
              rsv_.changeCapacity(sizeWithMargin);
            }
          }
        }
        return;
      }

      const auto size = this->size();
      if (minCapacity < size) {
        minCapacity = size;
      }

      if (loW > loW_old) {
        this->shrink_to_fit();
        const auto diffW = loW - loW_old;
        { // Update lower part.
          WBitsVec wbv_new(loW, minCapacity);
          wbv_new.resize(size);
          for (uint64_t i = 0, rvPos = 0; i < size; ++i, ++rvPos) {
            rvPos = rsv_.succ_1(rvPos);
            const auto val = ((rvPos - i) << loW_old) + wbv_.read(i);
            wbv_new.write(val & bits::UINTW_MAX(loW), i);
          }
          wbv_ = std::move(wbv_new);
        }
        { // Update higher part.
          const size_t rsv_new_size = size + (rsv_.getNum_0() >> diffW);
          RSDicT rsv_new(static_cast<uint64_t>(rsv_new_size * marginFactor));
          for (uint64_t i = 0, rvPos = 0, cur = 0; i < size; ++i, ++rvPos) {
            rvPos = rsv_.succ_1(rvPos);
            const auto next = (rvPos - i) >> diffW;
            while (cur < next) {
              rsv_new.appendBit(false);
              ++cur;
            }
            rsv_new.appendBit(true);
          }
          rsv_ = std::move(rsv_new);
        }
      } else {
        // rsv_.changeCapacity();
        const auto diffW = loW_old - loW;
        { // Update higher part.
          const size_t rsv_new_size = size + (rsv_.getNum_0() << diffW) + (wbv_.read(size - 1) >> loW);
          RSDicT rsv_new(static_cast<uint64_t>(rsv_new_size * marginFactor));
          for (uint64_t i = 0, rvPos = 0, cur = 0; i < size; ++i, ++rvPos) {
            rvPos = rsv_.succ_1(rvPos);
            const auto next = ((rvPos - i) << diffW) + (wbv_.read(i) >> loW);
            while (cur < next) {
              rsv_new.appendBit(false);
              ++cur;
            }
            rsv_new.appendBit(true);
          }
          rsv_ = std::move(rsv_new);
        }
        { // Update lower part.
          wbv_.convert(loW, minCapacity, doShrink);
        }
      }
    }


    /*!
     * @brief Shrink vector to fit current bit-length in use.
     */
    void shrink_to_fit() {
      wbv_.changeCapacity(0);
      rsv_.changeCapacity(0);
    }


    void printStatistics
    (
     const bool verbose = false
     ) const noexcept {
      std::cout << "SVec object (" << this << ") " << __func__ << "(" << verbose << ") BEGIN" << std::endl;
      const auto size = this->size();
      std::cout << "size = " << size << ", capacity = " << this->capacity() << ", loW = " << (uint64_t)(wbv_.getW());
      if (size) {
        const auto max = this->getMax();
        const auto bits = bits::bitSize(max);
        std::cout << ", max = " << max << ", maxW = " << (uint64_t)(bits) << ", num_zeros = " << this->getNum_0() << ", num_ones = " << this->getNum_1();
      }
      std::cout << std::endl;
      std::cout << this->calcMemBytes() << " bytes (hi = " << rsv_.calcMemBytes() << ", lo = " << wbv_.calcMemBytes() << ")" << std::endl;
      if (verbose) {
        wbv_.printStatistics(verbose);
        rsv_.printStatistics(verbose);
      }
      std::cout << "SVec object (" << this << ") " << __func__ << "(" << verbose << ") END" << std::endl;
    }
  };
} // namespace itmmti

#endif
