/*!
 * Copyright (c) 2017 Tomohiro I
 *
 * This program is released under the MIT License.
 * http://opensource.org/licenses/mit-license.php
 */
/*!
 * @file WBitsVec.hpp
 * @brief W-bits packed vector and its iterator.
 * @author Tomohiro I
 * @date 2017-01-24
 */
#ifndef INCLUDE_GUARD_WBitsVec
#define INCLUDE_GUARD_WBitsVec

#include <stdint.h> // include uint64_t etc.
#include <assert.h>

#include <iostream>
#include <iterator>
#include <algorithm>

#include "BitsUtil.hpp"
#include "ArrayUtil.hpp"
#include "MemUtil.hpp"

namespace itmmti
{
  class WBitsVec;


  /*!
   * @brief Iterator for WBitsVec.
   * @attention
   *   An iterator is invalidated without notification when WBitsVec object is destroyed or reallocated.
   *   Using invalidated iterator would cause problems.
   * @note
   *   Unfortunately packed vector is only qualified to be input_iterator due to the current design of STL.
   *   (see e.g. http://www.boost.org/doc/libs/1_63_0/libs/iterator/doc/index.html).
   */
  class WBitsVecIterator :
    public std::iterator<std::input_iterator_tag, uint64_t>
  {
    uint64_t * array_;
    uint64_t pos_;
    uint8_t w_;


  public:
    /*!
     * @brief Constructor.
     */
    WBitsVecIterator(uint64_t * array, uint64_t pos, uint8_t w) noexcept
    : array_(array), pos_(pos), w_(w)
    {
      assert(0 < w && w <= 64);
    }
    //// copy
    WBitsVecIterator(const WBitsVecIterator & itr) noexcept = default;
    WBitsVecIterator& operator=(const WBitsVecIterator & itr) noexcept = default;
    //// move
    WBitsVecIterator(WBitsVecIterator && itr) noexcept = default;
    WBitsVecIterator& operator=(WBitsVecIterator && itr) noexcept = default;


    /*!
     * @brief Read value at the position pointed by iterator.
     */
    uint64_t read() const {
      return bits::readWBits(array_, pos_, w_, bits::UINTW_MAX(w_));
    }


    /*!
     * @brief Write 'val' at the position pointed by iterator.
     */
    void write
    (
     const uint64_t val
     ) {
      bits::writeWBits(val, array_, pos_, w_, bits::UINTW_MAX(w_));
    }


    ////// operator
    //// *itr
    uint64_t operator*() const {
      return this->read();
    }

    //// ++itr
    WBitsVecIterator & operator++() {
      pos_ += w_;
      return *this;
    }

    //// itr++
    WBitsVecIterator operator++(int) {
      WBitsVecIterator tmp(*this);
      ++(*this);
      return tmp;
    }

    //// --itr
    WBitsVecIterator & operator--() {
      pos_ -= w_;
      return *this;
    }

    //// itr--
    WBitsVecIterator operator--(int) {
      WBitsVecIterator tmp(*this);
      --(*this);
      return tmp;
    }

    //// *this != itr
    bool operator!=(const WBitsVecIterator & itr) {
      return (pos_ != itr.pos_ || array_ != itr.array_);
    }

    //// *this == itr
    bool operator==(const WBitsVecIterator & itr) {
      return !(*this != itr);
    }

    ////// add more operators
    //// itr += diff
    WBitsVecIterator & operator+=(const difference_type diff) noexcept {
      assert(static_cast<int64_t>(pos_) + static_cast<int64_t>(w_) * diff >= 0);

      pos_ = static_cast<uint64_t>(static_cast<int64_t>(pos_) + static_cast<int64_t>(w_) * diff);
      //    pos_ += static_cast<int64_t>(w_) * diff;
      return *this;
    }

    //// itr -= diff
    WBitsVecIterator & operator-=(const difference_type diff) noexcept {
      *this += (-1 * diff);
      return *this;
    }

    //// itr + diff
    friend WBitsVecIterator operator+(const WBitsVecIterator & itr, const difference_type diff) noexcept {
      const int64_t pos = static_cast<int64_t>(itr.pos_) + static_cast<int64_t>(itr.w_) * diff;
      assert(pos >= 0);
      return WBitsVecIterator(itr.array_, static_cast<uint64_t>(pos), itr.w_);
    }

    //// diff + itr
    friend WBitsVecIterator operator+(const difference_type diff, const WBitsVecIterator & itr) noexcept {
      return itr + diff;
    }

    //// itr - diff
    friend WBitsVecIterator operator-(const WBitsVecIterator & itr, const difference_type diff) noexcept {
      return itr + (-1 * diff);
    }

    //// lhs - rhs
    friend difference_type operator-(const WBitsVecIterator & lhs, const WBitsVecIterator & rhs) noexcept {
      assert(lhs.w_ == rhs.w_);
      return (static_cast<int64_t>(lhs.pos_) - static_cast<int64_t>(rhs.pos_)) / lhs.w_;
    }

    //// lhs < rhs
    friend bool operator<(const WBitsVecIterator & lhs, const WBitsVecIterator & rhs) noexcept {
      return (lhs.pos_ < rhs.pos_ || lhs.array_ < rhs.array_);
    }

    //// rhs < lhs
    friend bool operator>(const WBitsVecIterator & lhs, const WBitsVecIterator & rhs) noexcept {
      return (rhs < lhs);
    }

    //// lhs <= rhs
    friend bool operator<=(const WBitsVecIterator & lhs, const WBitsVecIterator & rhs) noexcept {
      return !(lhs > rhs);
    }

    //// lhs >= rhs
    friend bool operator>=(const WBitsVecIterator & lhs, const WBitsVecIterator & rhs) noexcept {
      return !(lhs < rhs);
    }


    /*!
     * @brief Move values from src-region to tgt-region with same bit-width.
     * @pre The bit-width of src and tgt should be same.
     */
    friend void mvWBA_SameW
    (
     WBitsVecIterator & src, //!< Iterator specifying the beginning position of src.
     WBitsVecIterator & tgt, //!< Iterator specifying the beginning position of tgt.
     const uint64_t num //!< Number of elements to move.
     ) {
      assert(src.w_ == tgt.w_);
      if (src.pos_ >= tgt.pos_) {
        bits::mvBitsLR(src.array_, src.pos_, tgt.array_, tgt.pos_, num * src.w_);
      } else {
        bits::mvBitsRL(src.array_, src.pos_, tgt.array_, tgt.pos_, num * src.w_);
      }
    }


    /*!
     * @brief Move values from src-region to tgt-region with different bit-width.
     */
    friend void mvWBA_DiffW
    (
     WBitsVecIterator & src, //!< Iterator specifying the beginning position of src.
     WBitsVecIterator & tgt, //!< Iterator specifying the beginning position of tgt.
     const uint64_t num //!< Number of elements to move.
     ) {
      for (uint64_t i = 0; i < num; ++i, ++src, ++tgt) {
        assert(src.read() <= bits::UINTW_MAX(tgt.w_));
        tgt.write(src.read());
      }
    }


    /*!
     * @brief Move values from src-region to tgt-region.
     * @attention When src-region and tgt-region overlap, the overlapped part of src-region is overwritten.
     *            The other part of src-region is not changed.
     * @note
     *   The bit-width of src and tgt can be different.
     */
    friend void mvWBA
    (
     WBitsVecIterator & src, //!< Iterator specifying the beginning position of src.
     WBitsVecIterator & tgt, //!< Iterator specifying the beginning position of tgt.
     const uint64_t num //!< Number of elements to move.
     ) {
      if (src.w_ == tgt.w_) {
        mvWBA_SameW(src, tgt, num);
      } else {
        mvWBA_DiffW(src, tgt, num);
      }
    }


    friend void mvWBA_SameW(WBitsVecIterator && src, WBitsVecIterator && tgt, const uint64_t num) {
      assert(src.w_ == tgt.w_);

      if (src.array_ + src.pos_ >= tgt.array_ + tgt.pos_) {
        bits::mvBitsLR(src.array_, src.pos_, tgt.array_, tgt.pos_, num * src.w_);
      } else {
        bits::mvBitsRL(src.array_, src.pos_, tgt.array_, tgt.pos_, num * src.w_);
      }
    }


    friend void mvWBA_DiffW(WBitsVecIterator && src, WBitsVecIterator && tgt, const uint64_t num) {
      for (uint64_t i = 0; i < num; ++i, ++src, ++tgt) {
        assert(src.read() <= bits::UINTW_MAX(tgt.w_));
        tgt.write(src.read());
      }
    }


    friend void mvWBA(WBitsVecIterator && src, WBitsVecIterator && tgt, const uint64_t num) {
      if (src.w_ == tgt.w_) {
        mvWBA_SameW(src, tgt, num);
      } else {
        mvWBA_DiffW(src, tgt, num);
      }
    }
  };




  /*!
   * @brief W-bits packed vector. Bit-width 'w' and capacity can be changed dynamically.
   * @attention
   *   For technical reason, capacity is limited to '2^58 - 1' so that 'capacity * w' does not overflow.
   */
  class WBitsVec
  {
  public:
    using iterator = WBitsVecIterator;


  private:
    uint64_t * array_; //!< Array to store values.
    size_t capacity_; //!< Current capacity (must be in [0, 2^58)). It is guaranteed that the reserved space can accomodate 'capacity_ * w_' bits.
    size_t size_; //!< Current size (must be in [0, capacity_]).
    uint8_t w_; //!< Current bit-width (must be in [1, 64]).


  public:
    WBitsVec
    (
     uint8_t w = 1, //!< Initial bit-width.
     size_t capacity = 0 //!< Initial capacity.
     ) : array_(nullptr), capacity_(0), size_(0), w_(w) {
      assert(capacity <= ctcbits::UINTW_MAX(58));
      assert(0 < w && w <= 64);

      changeCapacity(capacity);
    }


    ~WBitsVec()
    {
      free(array_);
    }


    //// Copy
    /*!
     * @brief Copy constructor.
     * @attention
     *   Since the contents of 'other' are copied, it may take time when other.size_ is large.
     */
    WBitsVec
    (
     const WBitsVec & other
     ) : array_(nullptr), capacity_(0), size_(other.size_), w_(other.w_) {
      if (size_ > 0) {
        changeCapacity(other.capacity());
        bits::cpBytes(other.array_, array_, (size_ * w_ + 7) / 8);
      }
    }


    /*!
     * @brief Assignment operator.
     * @attention
     *   - It may take time when other.size_ is large.
     *   - The original contents of "this" are freed.
     */
    WBitsVec & operator=
    (
     const WBitsVec & other
     ) {
      if (this != &other) {
        clear(); changeCapacity(0); // clear() & changeCapacity(0) free array_
        size_ = other.size_;
        w_ = other.w_;
        changeCapacity(other.capacity_);
        if (size_ > 0) {
          bits::cpBytes(other.array_, array_, (size() * w_ + 7) / 8);
        }
      }
      return *this;
    }


    //// Move
    /*!
     * @brief Move constructor.
     * @attention
     *   "other" is initialized to an object with capacity = 0.
     */
    WBitsVec
    (
     WBitsVec && other
     ) : array_(other.array_), capacity_(other.capacity_), size_(other.size_), w_(other.w_) {
      other.array_ = nullptr;
      other.size_ = other.capacity_ = 0;
    }


    /*!
     * @brief Move operator.
     * @attention
     *   - The original contents of "this" are freed.
     *   - "other" is initialized to an object with capacity = 0.
     */
    WBitsVec & operator=
    (
     WBitsVec && other
     ) {
      if (this != &other) {
        clear(); changeCapacity(0); // clear() & changeCapacity(0) free array_
        array_ = other.array_;
        capacity_ = other.capacity_;
        size_ = other.size_;
        w_ = other.w_;
        other.array_ = nullptr;
        other.size_ = other.capacity_ = 0;
      }
      return *this;
    }


    /*!
     * @brief Array subscript for writing
     */
    array_util::PackedArrayTypeValRef<WBitsVec> operator[]
    (
     const size_t idx
     ) {
      return array_util::PackedArrayTypeValRef<WBitsVec>(this, idx);
    }


    /*!
     * @brief Read only accessor.
     */
    uint64_t operator[]
    (
     const size_t idx //!< in [0, capacity_).
     ) const {
      assert(idx < capacity_);
    
      return this->read(idx);
    }


    /*!
     * @brief Get iterator at given idx.
     */
    WBitsVec::iterator getItrAt
    (
     size_t idx
     ) noexcept {
      return WBitsVecIterator(array_, idx * w_, w_);
    }


    /*!
     * @brief Get iterator points to the first element.
     */
    WBitsVec::iterator begin() noexcept {
      return getItrAt(0);
    }


    /*!
     * @brief Get iterator points to the end (right next to the last element).
     */
    WBitsVec::iterator end() noexcept {
      return getItrAt(size_);
    }


    /*!
     * @brief Read value at 'idx'.
     */
    uint64_t read
    (
     const size_t idx //!< in [0, capacity_).
     ) const {
      assert(idx < capacity_);
    
      return bits::readWBits(array_, idx * w_, w_, bits::UINTW_MAX(w_));
    }


    /*!
     * @brief Write 'val' at 'idx'.
     */
    void write
    (
     const uint64_t val, //!< in [0, 2^{w_}).
     const size_t idx //!< in [0, capacity_).
     ) {
      assert(idx < capacity_);
      assert(val <= bits::UINTW_MAX(w_));

      bits::writeWBits(val, array_, idx * w_, w_, bits::UINTW_MAX(w_));
    }


    /*!
     * @brief Get current bit-width of each element.
     */
    uint8_t getW() const noexcept {
      return w_;
    }


    /*!
     * @brief Get current capacity.
     */
    size_t capacity() const noexcept {
      return capacity_;
    }


    /*!
     * @brief Get current size.
     */
    size_t size() const noexcept {
      return size_;
    }


    /*!
     * @brief Calculate total memory usage in bytes.
     */
    size_t calcMemBytes() const noexcept {
      return sizeof(*this) + sizeof(uint64_t) * ((capacity_ * w_ + 63) / 64);
    }


    /*!
     * @brief Return if no element is stored.
     */
    bool empty() const noexcept {
      return (size_ == 0);
    }


    /*!
     * @brief Clear vector. It only changes size to zero.
     */
    void clear() noexcept {
      size_ = 0;
    }


    /*!
     * @brief Resize 'size_' to 'newSize'. It expands 'array_' if needed.
     * @note
     *   This function does not shrink 'capacity_'
     *   If 'newSize > capacity_', expand 'array_' by calling WBitsVec::reserve.
     */
    void resize
    (
     const size_t newSize
     ) {
      assert(newSize <= ctcbits::UINTW_MAX(58));

      if (newSize > capacity_) {
        changeCapacity(newSize);
      }
      size_ = newSize;
    }


    /*!
     * @brief Variant of WBitsVec::resize: If 'newSize > capacity_', it just returns false. Otherwise resize and return true.
     */
    bool resizeWithoutReserve
    (
     const size_t newSize
     ) noexcept {
      assert(newSize <= ctcbits::UINTW_MAX(58));

      if (newSize <= capacity_) {
        size_ = newSize;
        return true;
      }
      return false;
    }


    /*!
     * @brief Change capacity to max of givenCapacity and current size.
     * @node If givenCapacity is 0, it works as shrink_to_fit.
     */
    void changeCapacity
    (
     const size_t givenCapacity
     ) {
      if (capacity_ != givenCapacity) {
        const auto newCapacity = std::max(size_, givenCapacity);
        if (newCapacity > 0) {
          const size_t newLen = (newCapacity * w_ + 63) / 64; // +63 for roundup
          memutil::realloc_AbortOnFail<uint64_t>(array_, newLen);
          capacity_ = newLen * 64 / w_;
        } else {
          memutil::safefree(array_);
          capacity_ = 0;
        }
      }
    }


    /*!
     * @brief Change 'w_', and values accordingly. Give optional arguments ('minCapacity', 'doShrink') to control reallocation of array.
     * @note
     *   This function will do the following two tasks:
     *   1. Convert values of w_ bits to those of 'w' bits, and renew 'w_'.
     *      Note that If 'w < w_', the 'w_ - w' most significant bits for each value will be discarded.
     *   2. Realloc 'array_' if needed/required.
     *      If doShrink is true, 'array_' is reallocated to fit 'max(size_, minCapacity)' 'w'-bits uints.
     *      Otherwise, reallocation is executed only when array expands.
     * @note
     *   'size_' will not change.
     *
     * @par TIPS
     *   We do not provide a function to change only 'w_' (i.e., the setter for 'w_')
     *   because most of the case changing 'w_' would be accompanied by conversion (and reallocation if needed).
     *   If you want to change 'w_' but do not want to convert nor realloc, do the followings:
     *   @code
     *     WBAObj.clear();
     *     WBAObj.convert(w);
     *   @endcode
     */
    void convert
    (
     const uint8_t w, //!< New bit-width.
     size_t minCapacity = 0, //!< Minimum support for capacity (default = 0, which turns into 'size_').
     bool doShrink = false //!< If true, 'array_' is reallocated to fit max(size_, minCapacity) w-bits uints.
     ) {
      assert(0 < w && w <= 64);
      assert(minCapacity <= ctcbits::UINTW_MAX(58));

      // Convert the values in array for w < w_ where needed bits become smaller than before.
      if (w < w_) { // w < w_ (w_ - w most significant bits for each value will be discarded)
        for (uint64_t i = 0; i < this->size(); ++i) {
          bits::writeWBits(this->read(i) & bits::UINTW_MAX(w), array_, i * w, w, bits::UINTW_MAX(w));
        }
      }

      if (minCapacity < size_) {
        minCapacity = size_;
      }
      const size_t oldLen = (capacity_ * w_ + 63) / 64; // +63 for roundup
      const size_t minLen = (minCapacity * w + 63) / 64; // +63 for roundup
      if (minLen > oldLen || doShrink) {
        if (minLen) {
          memutil::realloc_AbortOnFail(array_, minLen);
        } else {
          memutil::safefree(array_);
        }
        capacity_ = (minLen * 64) / w;
      } else if (!doShrink) {
        capacity_ = (oldLen * 64) / w;
      }

      // Convert the values in array for w > w_ where needed bits become larger than before.
      if (w > w_) {
        for (uint64_t i = this->size() - 1; i != UINT64_MAX; --i) {
          bits::writeWBits(this->read(i), array_, i * w, w, bits::UINTW_MAX(w));
        }
      }
      w_ = w;
    }


    void printStatistics
    (
     const bool verbose = false
     ) const noexcept {
      std::cout << "WBitsVec object (" << this << ") " << __func__ << "(" << verbose << ") BEGIN" << std::endl;
      std::cout << "size = " << this->size() << ", capacity = " << this->capacity() << std::endl;
      if (verbose) {
        const auto size = this->size();
        std::cout << "dump stored values" << std::endl;
        for (uint64_t i = 0; i < size; ++i) {
          std::cout << this->read(i) << ", ";
        }
        std::cout << std::endl;
      }
      std::cout << "WBitsVec object (" << this << ") " << __func__ << "(" << verbose << ") END" << std::endl;
    }
  };
} // namespace itmmti

#endif
