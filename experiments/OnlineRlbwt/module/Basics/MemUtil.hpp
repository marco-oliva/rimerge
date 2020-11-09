#ifndef INCLUDE_GUARD_MemUtil
#define INCLUDE_GUARD_MemUtil

#include <stdint.h> // include uint64_t etc.
#include <assert.h>
#include <iostream>

namespace itmmti
{
  /*!
   * @namespace memutil
   * @brief Utilities for memory allocation.
   */
  namespace memutil
  {
    /*!
     * @brief malloc that aborts on fail.
     * @tparam T Type of elements to allocate.
     * @return Pointer to allocated object.
     * @par Example
     *   @code
     *   uint64_t * p = memutil::malloc_AbortOnFail<uint64_t>(10);
     *   @endcode
     *   The above code would reserve the space for 10 64bit uints, i.e., 10 * 8 = 80 bytes.
     */
    template <typename T>
    T * malloc_AbortOnFail
    (
     size_t n //!< Num of elements of 'T' to allocate.
     ) noexcept {
      assert(n > 0); //! @pre 'n' > 0.

      T * ptr = static_cast<T *>(malloc(n * sizeof(T)));
      if (ptr == NULL) {
        std::cerr << "ABORTED ON FAIL: " << __func__ << ", sizeof(T) = " << sizeof(T) << ", n = " << n << std::endl;
        exit(EXIT_FAILURE);
      }
      return ptr;
    }


    /*!
     * @brief realloc that aborts on fail.
     * @tparam T Type of elements to allocate.
     * @return Pointer to allocated object.
     * @par Example
     *   @code
     *   uint64_t * p = memutil::malloc_AbortOnFail<uint64_t>(10);
     *   // Do something on p.
     *   memutil::realloc_AbortOnFail(p, 20); // Expand p.
     *   @endcode
     */
    template <typename T>
    void realloc_AbortOnFail
    (
     T *& ptr,
     size_t n //!< Num of elements of 'T' to allocate.
     ) noexcept {
      assert(n > 0); //! @pre 'n' > 0.

      ptr = static_cast<T *>(realloc(ptr, n * sizeof(T)));
      if (ptr == nullptr) {
        std::cerr << "ABORTED ON FAIL: " << __func__ << ", sizeof(T) = " << sizeof(T) << ", n = " << n << std::endl;
        exit(EXIT_FAILURE);
      }
    }


    /*!
     * @brief Safe free that sets NULL after free.
     */
    template <typename T>
    void safefree(T *& ptr) {
      free(ptr);
      ptr = nullptr;
    }


    /*!
     * @brief Safe delete that sets NULL after delete.
     */
    template <typename T>
    void safedelete(T *& ptr) {
      delete ptr;
      ptr = nullptr;
    }
  } // namespace memutil
} // namespace itmmti

#endif
