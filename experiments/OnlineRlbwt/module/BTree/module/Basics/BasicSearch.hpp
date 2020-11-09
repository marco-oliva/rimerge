/*!
 * Copyright (c) 2017 Tomohiro I
 *
 * This program is released under the MIT License.
 * http://opensource.org/licenses/mit-license.php
 */
/*!
 * @file BasicSearch.hpp
 * @brief Useful functions to work on bits.
 * @author Tomohiro I
 * @date 2017-08-20
 */
#ifndef INCLUDE_GUARD_BasicSearch
#define INCLUDE_GUARD_BasicSearch

#include <stdint.h> // include uint64_t etc.
#include <assert.h>


/*!
 * @namespace basic_search
 * @brief Basic search on array.
 */
namespace basic_search
{
  /*!
   * @brief My implimentation of lower_bound working on ranges with predicate function.
   * @tparam LS: Thereshold to switch to liniear search.
   * @tparam
   *   Predicate: Fuction that returns true/false for a given integer in [lb, ub).
   *   The function should partition [lb, ub) into two such that the first part is all false
   *   (first part is possibly empty) and the second part is all true.
   * @return The smallest idx in [lb, ub) such that func(idx) = true.
   */
  template<uint64_t LS = 0, class Predicate>
    uint64_t partition_idx
  (
   uint64_t lb, //!< Initial lower bound (inclusive).
   uint64_t ub, //!< Initial upper bound (exclusive).
   Predicate func
   ) noexcept {
    assert(lb < ub);
    assert(ub <= UINT64_MAX - LS);

    if (LS != UINT64_MAX) { // If LS == UINT64_MAX, it is just a linear search.
      while (lb + LS + 1 < ub) { // invariant: the answer is in [lb..ub).
        auto mid = (lb + ub) / 2;
        if (!func(mid - 1)) {
          lb = mid;
        } else {
          ub = mid;
        }
      }
    }
    if (LS) {
      while (!func(lb)) {
        ++lb;
      }
    }

    return lb;
  }
}

#endif
