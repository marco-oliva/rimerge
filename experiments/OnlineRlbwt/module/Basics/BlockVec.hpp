/*!
 * Copyright (c) 2018 Tomohiro I
 *
 * This program is released under the MIT License.
 * http://opensource.org/licenses/mit-license.php
 */
/*!
 * @file BlockVec.hpp
 * @brief Vector of blocks
 * @author Tomohiro I
 * @date 2018-01-22
 */
#ifndef INCLUDE_GUARD_BlockVec
#define INCLUDE_GUARD_BlockVec

#include <stdint.h> // include uint64_t etc.
#include <assert.h>

#include <iostream>
#include <iterator>
#include <algorithm>

#include "BitsUtil.hpp"
#include "MemUtil.hpp"

namespace itmmti
{
  /*!
   * @brief
   * @pram Specification
   *   - BlockVec is not responsible to freeing each block, i.e., blocks_[i]
   *   - BlockVec is not responsible to calculating memory for each block
   */
  template<typename ElemT, uint64_t tparam_kBlockSize>
  class BlockVec
  {
  public:
    // Public constant, alias etc.
    static constexpr uint64_t kBlockSize{tparam_kBlockSize};


  private:
    ElemT ** blocks_; //!< Pointers to blocks
    size_t capacity_; //!< Current capacity
    size_t size_; //!< Current size (must be in [0..capacity_]).
    size_t numBlocksCapacity_; //!< Current capacity of blocks_.
    size_t numBlocks_; //!< Current size of blocks_ (must be in [0..numBlocksCapacity_]).


  public:
    BlockVec
    (
     size_t givenCapacity = 0 //!< Initial capacity.
     ) : blocks_(nullptr), capacity_(0), size_(0), numBlocksCapacity_(0), numBlocks_(0) {
      while (givenCapacity > capacity_) {
        appendBlock();
      }
    }


    ~BlockVec()
    {
      freeBlocksPtrArray();
    }


    void freeBlocksPtrArray() {
      memutil::safefree(blocks_);
      capacity_ = size_ = numBlocksCapacity_ = numBlocks_ = 0;
    }


    //// Copy
    /*!
     * @brief Copy constructor.
     * @attention
     *   Since the contents of 'other' are copied, it may take time when other.size_ is large.
     */
    BlockVec
    (
     const BlockVec & other
     ) : blocks_(nullptr), capacity_(0), size_(0), numBlocksCapacity_(0), numBlocks_(0) {
      resize(other.size());
      std::cout << "size_ = " << size_ << " other.size() = " << other.size() << std::endl;
      for (size_t i = 0; i < size_; ++i) {
        (*this)[i] = other[i];
      }
    }


    /*!
     * @brief Assignment operator.
     * @attention
     *   - It may take time when other.size_ is large.
     *   - The original contents of "this" are freed.
     */
    BlockVec & operator=
    (
     const BlockVec & other
     ) {
      if (this != &other) {
        resize(other.size());
        for (size_t i = 0; i < size_; ++i) {
          (*this)[i] = other[i];
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
    BlockVec
    (
     BlockVec && other
     ) : blocks_(other.blocks_), capacity_(other.capacity_), size_(other.size_), numBlocksCapacity_(other.numBlocks_), numBlocks_(other.numBlocks_) {
      other.freeBlocksPtrArray();
    }


    /*!
     * @brief Move operator.
     * @attention
     *   - Potential memory leak: Be sure that it does not free each block in "*this".
     *   - "other" is initialized to an object with capacity = 0.
     */
    BlockVec & operator=
    (
     BlockVec && other
     ) {
      if (this != &other) {
        freeBlocksPtrArray(); // Potential memory leak: Be sure that it does not free each block in "*this"
        blocks_ = other.blocks_;
        capacity_ = other.capacity_;
        size_ = other.size_;
        other.freeBlocksPtrArray();
      }
      return *this;
    }


    /*!
     * @brief Read only accessor.
     */
    ElemT & operator[]
    (
     const size_t idx //!< in [0, capacity_).
     ) {
      assert(idx < capacity_);
    
      return blocks_[idx / kBlockSize][idx % kBlockSize];
    }


    /*!
     * @brief Read only accessor.
     */
    const ElemT & operator[]
    (
     const size_t idx //!< in [0, capacity_).
     ) const {
      assert(idx < capacity_);
    
      return blocks_[idx / kBlockSize][idx % kBlockSize];
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
     * @brief Get current capacity of "blocks_".
     */
    size_t getNumBlocksCapacity() const noexcept {
      return numBlocksCapacity_;
    }


    /*!
     * @brief Get current size of "blocks_".
     */
    size_t getNumBlocks() const noexcept {
      return numBlocks_;
    }

    
    /*!
     * @brief Get const block pointer
     */
    const ElemT * getConstBlockPtr
    (
     size_t i
     ) const noexcept {
      assert(i < numBlocksCapacity_);

      return blocks_[i];
    }


    /*!
     * @brief Get block pointer
     */
    ElemT * getBlockPtr
    (
     size_t i
     ) noexcept {
      assert(i < numBlocksCapacity_);

      return blocks_[i];
    }


    /*!
     * @brief Set block pointer at "i".
     */
    void setBlockPtr
    (
     size_t i,
     ElemT * ptr
     ) const noexcept {
      assert(i < numBlocksCapacity_);

      blocks_[i] = ptr;
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
      while (newSize > capacity_) {
        appendBlock();
      }
      size_ = newSize;
    }


    /*!
     * @brief Append a block and expand "capacity_" by kBlockSize.
     */
    void appendBlock() {
      if (++numBlocks_ > numBlocksCapacity_) {
        numBlocksCapacity_ = numBlocks_ * 2;
        memutil::realloc_AbortOnFail<ElemT *>(blocks_, numBlocksCapacity_);
      }
      blocks_[numBlocks_ - 1] = new ElemT[kBlockSize];
      capacity_ += kBlockSize;
    }


    /*!
     * @brief Append a block and expand "capacity_" by kBlockSize.
     */
    void appendBlock
    (
     ElemT * ptr
     ) {
      if (++numBlocks_ > numBlocksCapacity_) {
        numBlocksCapacity_ = numBlocks_ * 2;
        memutil::realloc_AbortOnFail<ElemT *>(blocks_, numBlocksCapacity_);
      }
      blocks_[numBlocks_ - 1] = ptr;
      capacity_ += kBlockSize;
    }


    /*!
     * @brief Append a block and expand "capacity_" by kBlockSize.
     * @attention It does not free blocks.
     */
    void reduceNumBlocks
    (
     const size_t givenNumBlocks
     ) {
      assert(givenNumBlocks <= numBlocks_);

      numBlocks_ = givenNumBlocks;
    }


    /*!
     * @brief Calculate total memory usage in bytes.
     * @attention It is only responsible to space for pointers in blocks_.
     */
    size_t calcMemBytes
    (
     bool includeThis = true
     ) const noexcept {
      size_t size = sizeof(*this) * includeThis;
      return size + sizeof(ElemT *) * numBlocksCapacity_;
    }


    void printStatistics
    (
     std::ostream & os,
     const bool verbose = false
     ) const noexcept {
      os << "BlockVec object (" << this << ") " << __func__ << "(" << verbose << ") BEGIN" << std::endl;
      os << "size = " << this->size() << ", capacity = " << this->capacity() << std::endl;
      os << "numBlocks = " << this->getNumBlocks() << ", numBlocksCapacity = " << this->getNumBlocksCapacity() << std::endl;
      // if (verbose) {
      //   const auto size = this->size();
      //   os << "dump stored values" << std::endl;
      //   for (uint64_t i = 0; i < size; ++i) {
      //     if (i % kBlockSize) {
      //       os << i % kBlockSize << ": (" << blocks_[i % kBlockSize] << ")" << std::endl;
      //     }
      //     os << (*this)[i] << ", ";
      //   }
      //   os << std::endl;
      // }
      os << "BlockVec object (" << this << ") " << __func__ << "(" << verbose << ") END" << std::endl;
    }
  };
} // namespace itmmti

#endif
