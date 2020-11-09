/*!
 * Copyright (c) 2017 Tomohiro I
 *
 * This program is released under the MIT License.
 * http://opensource.org/licenses/mit-license.php
 */
/*!
 * @file WBitsBlockVec.hpp
 * @brief W-bits packed vector with BlockVec
 * @author Tomohiro I
 * @date 2018-01-23
 */
#ifndef INCLUDE_GUARD_WBitsBlockVec
#define INCLUDE_GUARD_WBitsBlockVec

#include <stdint.h> // include uint64_t etc.
#include <assert.h>

#include <iostream>
#include <iterator>
#include <algorithm>

#include "BlockVec.hpp"
#include "BitsUtil.hpp"
#include "ArrayUtil.hpp"
#include "MemUtil.hpp"

namespace itmmti
{
  /*!
   * @brief W-bits packed vector. Bit-width 'w' and capacity can be changed dynamically.
   * @tparam tparam_kBlockSize in 2^6, 2^7, 2^8, ..., 2^63
   * @attention
   *   For technical reason, capacity is limited to '2^58 - 1' so that 'capacity * w' does not overflow.
   */
  template<uint64_t tparam_kBlockSize>
  class WBitsBlockVec
  {
  public:
    // Public constant, alias etc.
    static constexpr uint64_t kBlockSize{tparam_kBlockSize};
    using WBitsBlockVecT = WBitsBlockVec<tparam_kBlockSize>;
    using BlockVecT = BlockVec<uint64_t, tparam_kBlockSize>;


  private:
    BlockVecT vec_; //!< Abstracted array to store values.
    size_t capacity_; //!< Current capacity (must be in [0, 2^58)). It is guaranteed that the reserved space can accomodate 'capacity_ * w_' bits.
    size_t size_; //!< Current size (must be in [0, capacity_]).
    uint8_t w_; //!< Current bit-width (must be in [1, 64]).


  public:
    WBitsBlockVec
    (
     uint8_t w = 1, //!< Initial bit-width.
     size_t capacity = 0 //!< Initial capacity.
     ) : vec_(), capacity_(0), size_(0), w_(w) {
      assert(capacity <= ctcbits::UINTW_MAX(58));
      assert(0 < w && w <= 64);

      reserve(capacity);
    }


    ~WBitsBlockVec()
    {
      clearAll();
    }


    void clearAll() {
      for (uint64_t i = 0; i < vec_.getNumBlocks(); ++i) {
        delete[] vec_.getBlockPtr(i);
      }
      vec_.freeBlocksPtrArray();
      capacity_ = size_ = 0;
    }


    //// Copy
    /*!
     * @brief Copy constructor.
     * @attention
     *   Since the contents of 'other' are copied, it may take time when other.size_ is large.
     */
    WBitsBlockVec
    (
     const WBitsBlockVecT & other
     ) : vec_(other.vec_), size_(other.size_), w_(other.w_) {
      capacity_ = (vec_.capacity() * 64) / other.w_;
    }


    /*!
     * @brief Assignment operator.
     * @attention
     *   - It may take time when other.size_ is large.
     *   - The original contents of "this" are freed.
     */
    WBitsBlockVecT & operator=
    (
     const WBitsBlockVecT & other
     ) {
      if (this != &other) {
        this->vec_ = other.vec_; // clear() & changeCapacity(0) free array_
        capacity_ = (vec_.capacity() * 64) / other.w_;
        size_ = other.size_;
        w_ = other.w_;
      }
      return *this;
    }


    //// Move
    /*!
     * @brief Move constructor.
     * @attention
     *   "other" is initialized to an object with capacity = 0.
     */
    WBitsBlockVec
    (
     WBitsBlockVecT && other
     ) : vec_(other.vec_), capacity_(other.capacity_), size_(other.size_), w_(other.w_) {
      other.size_ = other.capacity_ = 0;
    }


    /*!
     * @brief Move operator.
     * @attention
     *   - The original contents of "this" are freed.
     *   - "other" is initialized to an object with capacity = 0.
     */
    WBitsBlockVecT & operator=
    (
     WBitsBlockVecT && other
     ) {
      if (this != &other) {
        vec_ = other.vec_; // move
        capacity_ = (vec_.capacity() * 64) / other.w_;
        size_ = other.size_;
        w_ = other.w_;
        other.size_ = other.capacity_ = 0;
      }
      return *this;
    }


    /*!
     * @brief Array subscript for writing
     */
    array_util::PackedArrayTypeValRef<WBitsBlockVec> operator[]
    (
     const size_t idx
     ) {
      return array_util::PackedArrayTypeValRef<WBitsBlockVec>(this, idx);
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
     * @brief Read value at 'idx'.
     */
    uint64_t read
    (
     const size_t idx //!< in [0, capacity_).
     ) const {
      assert(idx < capacity_);

      const auto bitPos = idx * w_;
      return bits::readWBits(vec_, bitPos, w_, bits::UINTW_MAX(w_));
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

      const auto bitPos = idx * w_;
      bits::writeWBits(val, vec_, bitPos, w_, bits::UINTW_MAX(w_));
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
    size_t calcMemBytes
    (
     bool includeThis = true
     ) const noexcept {
      size_t size = sizeof(*this) * includeThis;
      return size + vec_.calcMemBytes(false) + (sizeof(uint64_t) * vec_.getNumBlocks() * kBlockSize);
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

      const size_t numUint = (newSize * w_ + 63) / 64; // +63 for roundup
      vec_.resize(numUint);
      size_ = newSize;
      capacity_ = (vec_.capacity() * 64) / w_;
    }


    /*!
     * @brief Reserve
     */
    void reserve
    (
     const size_t givenCapacity
     ) {
      assert(givenCapacity <= ctcbits::UINTW_MAX(58));

      const size_t numUint = (givenCapacity * w_ + 63) / 64; // +63 for roundup
      while (numUint > vec_.capacity()) {
        vec_.appendBlock();
      }
      capacity_ = (vec_.capacity() * 64) / w_;
    }


    /*!
     * @brief Append a block and expand "capacity_".
     */
    void appendBlock() {
      vec_.appendBlock();
      capacity_ = (vec_.capacity() * 64) / w_;
    }


    /*!
     * @brief Append a block and expand "capacity_" by kBlockSize.
     */
    void appendBlock
    (
     uint64_t * ptr
     ) {
      vec_.appendBlock(ptr);
      capacity_ = (vec_.capacity() * 64) / w_;
    }


    /*!
     * @brief Increase "w_" to "w", and change values accordingly.
     */
    void increaseW
    (
     const uint8_t w //!< New bit-width.
     ) {
      assert(0 < w && w <= 64);
      assert(w_ <= w);

      const size_t numUint = (size_ * w + 63) / 64; // +63 for roundup
      vec_.resize(numUint);

      // Convert the values in array for w > w_ where needed bits become larger than before.
      for (uint64_t i = this->size() - 1; i != UINT64_MAX; --i) {
        bits::writeWBits(this->read(i), vec_, i * w, w, bits::UINTW_MAX(w));
      }
      w_ = w;
      capacity_ = (vec_.capacity() * 64) / w_;
    }


    /*!
     * @brief Decrease "w_" to "w", and change values accordingly.
     */
    void decreaseW
    (
     const uint8_t w //!< New bit-width.
     ) {
      assert(0 < w && w <= 64);
      assert(w <= w_);

      // Convert the values in array for w < w_ where needed bits become smaller than before.
      // w < w_ (w_ - w most significant bits for each value will be discarded)
      for (uint64_t i = 0; i < this->size(); ++i) {
        bits::writeWBits(this->read(i) & bits::UINTW_MAX(w), vec_, i * w, w, bits::UINTW_MAX(w));
      }
      w_ = w;
      capacity_ = (vec_.capacity() * 64) / w_;
    }


    /*!
     * @brief Change capacity and discard unused blocks.
     */
    void shrink
    (
     const size_t givenCapacity
     ) {
      assert(size_ <= givenCapacity);
      assert(givenCapacity < capacity_);

      const size_t needNumBlocks = (((givenCapacity * w_ + 63) / 64) + (kBlockSize - 1) / kBlockSize);
      for (uint64_t i = needNumBlocks; i < vec_.getNumBlocks(); ++i) {
        free(vec_.getBlockPtr(i));
      }
      capacity_ = (vec_.capacity() * 64) / w_;
    }


    void printStatistics
    (
     const bool verbose = false
     ) const noexcept {
      std::cout << "WBitsBlockVec object (" << this << ") " << __func__ << "(" << verbose << ") BEGIN" << std::endl;
      std::cout << "size = " << this->size() << ", capacity = " << this->capacity() << std::endl;
      if (verbose) {
        const auto size = this->size();
        std::cout << "dump stored values" << std::endl;
        for (uint64_t i = 0; i < size; ++i) {
          std::cout << this->read(i) << ", ";
        }
        std::cout << std::endl;
      }
      std::cout << "WBitsBlockVec object (" << this << ") " << __func__ << "(" << verbose << ") END" << std::endl;
    }
  };
} // namespace itmmti

#endif
