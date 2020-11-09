/*!
 * Copyright (c) 2018 Tomohiro I
 *
 * This program is released under the MIT License.
 * http://opensource.org/licenses/mit-license.php
 */
/*!
 * @file ArrayUtil.hpp
 * @brief Useful functions to work on array.
 * @author Tomohiro I
 * @date 2018-01-29
 */
#ifndef INCLUDE_GUARD_ArrayUtil
#define INCLUDE_GUARD_ArrayUtil

#include <stdint.h> // include uint64_t etc.
#include <assert.h>
#include <iostream>
#include <cstring>

namespace itmmti
{
  /*!
   * @namespace array_util
   * @brief Useful functions to work on array.
   */
  namespace array_util
  {
    /*!
     * @brief Proxiy class to call write function with [] operator for packed array data type.
     * @attention
     *   It becomes roughly 2-times slower than direct read/write if compiler optimization is disabled (recommend at least -O2).
     */
    template<typename ArrayT>
    class PackedArrayTypeValRef
    {
      friend ArrayT;


    private:
      ArrayT * const obj_;
      const uint64_t idx_;


      PackedArrayTypeValRef
      (
       ArrayT * obj,
       uint64_t idx
       ) :
        obj_(obj),
        idx_(idx)
      {}


    public:
      uint64_t operator=
      (
       uint64_t val
       ) {
        obj_->write(val, idx_);
        return val;
      }


      operator uint64_t() const {
        return obj_->read(idx_);
      }
    };




    /*!
     * @brief Move values (by scanning from left to right) from src to tgt.
     * @note Read function and write function can be given
     */
    template<typename ReadFunc, typename WriteFunc>
    inline void mvValsLR
    (
     const uint64_t srcIdx,
     const uint64_t tgtIdx,
     uint64_t num,
     ReadFunc readFunc,
     WriteFunc writeFunc
     ) {
      for (uint64_t i = 0; i < num; ++i) {
        writeFunc(readFunc(srcIdx + i), tgtIdx + i);
      }
    }


    /*!
     * @brief Move values (by scanning from right to left) from src to tgt.
     * @note Read function and write function can be given
     */
    template<typename ReadFunc, typename WriteFunc>
    inline void mvValsRL
    (
     const uint64_t srcIdx,
     const uint64_t tgtIdx,
     uint64_t num,
     ReadFunc readFunc,
     WriteFunc writeFunc
     ) {
      for (uint64_t i = num - 1; i != UINT64_MAX; --i) {
        writeFunc(readFunc(srcIdx + i), tgtIdx + i);
      }
    }


    /*!
     * @brief Move values from src to tgt.
     * @note Read function and write function can be given
     */
    template<typename ReadFunc, typename WriteFunc>
    inline void mvVals
    (
     const uint64_t srcIdx,
     const uint64_t tgtIdx,
     uint64_t num,
     ReadFunc readFunc,
     WriteFunc writeFunc
     ) {
      if (srcIdx >= tgtIdx) {
        mvValsLR(srcIdx, tgtIdx, num, readFunc, writeFunc);
      } else {
        mvValsRL(srcIdx, tgtIdx, num, readFunc, writeFunc);
      }
    }

    
    /*!
     * @brief Move values in array type from src to tgt while scanning left to right
     */
    template<typename SrcArrayT, typename TgtArrayT>
    inline void mvValsLR
    (
     const SrcArrayT & src,
     const uint64_t srcIdx,
     TgtArrayT & tgt,
     const uint64_t tgtIdx,
     uint64_t num
     ) {
      mvValsLR(srcIdx, tgtIdx, num,
               [src](uint64_t i) {
                 return src[i];
               },
               [tgt](decltype(tgt[0]) val, uint64_t i) {
                 tgt[i] = val;
               });
    }


    /*!
     * @brief Move values in array type from src to tgt while scanning right to left
     */
    template<typename SrcArrayT, typename TgtArrayT>
    inline void mvValsRL
    (
     const SrcArrayT & src,
     const uint64_t srcIdx,
     TgtArrayT & tgt,
     const uint64_t tgtIdx,
     uint64_t num
     ) {
      mvValsRL(srcIdx, tgtIdx, num,
               [src](uint64_t i) {
                 src[i];
               },
               [tgt](decltype(tgt[0]) val, uint64_t i) {
                 tgt[i] = val;
               });
    }


    /*!
     * @brief Move values in array type from src to tgt
     */
    template<typename SrcArrayT, typename TgtArrayT>
    inline void mvVals
    (
     const SrcArrayT & src,
     const uint64_t srcIdx,
     TgtArrayT & tgt,
     const uint64_t tgtIdx,
     uint64_t num
     ) {
      mvVals(srcIdx, tgtIdx, num,
             [src](uint64_t i) {
               src[i];
             },
             [tgt](decltype(tgt[0]) val, uint64_t i) {
               tgt[i] = val;
             });
    }
  }
} // namespace itmmti


#endif
