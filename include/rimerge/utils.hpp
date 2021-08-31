/*
  Copyright (c) 2018 Jouni Siren

  Author: Jouni Siren <jouni.siren@iki.fi>

  Permission is hereby granted, free of charge, to any person obtaining a copy
  of this software and associated documentation files (the "Software"), to deal
  in the Software without restriction, including without limitation the rights
  to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
  copies of the Software, and to permit persons to whom the Software is
  furnished to do so, subject to the following conditions:

  The above copyright notice and this permission notice shall be included in all
  copies or substantial portions of the Software.

  THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
  IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
  FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
  AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
  LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
  OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
  SOFTWARE.
*/

#ifndef rimerge_utils_hpp
#define rimerge_utils_hpp

#include <algorithm>
#include <cstdint>
#include <fstream>
#include <iostream>
#include <string>
#include <vector>
#include <array>
#include <map>
#include <algorithm>
#include <mutex>
#include <atomic>
#include <limits>
#include <stdio.h>


#include <sdsl/bit_vectors.hpp>
#include <sdsl/wavelet_trees.hpp>
#include <sdsl/suffix_arrays.hpp>


#include <spdlog/spdlog.h>

#include <mio/mio.hpp>

#include <rimerge/malloc_count.h>

#include <omp.h>

// Parallel sorting is only available with libstdc++ parallel mode.
#ifdef __GLIBCXX__
#include <parallel/algorithm>
#endif

namespace rimerge
{

typedef unsigned long long size_type;
typedef unsigned long long rank_type;
typedef unsigned short_type;
typedef unsigned char  byte_type;

//------------------------------------------------------------------------------

constexpr size_type BYTE_BITS    = 8;
constexpr size_type S_WORD_BITS  = 32;
constexpr size_type WORD_BITS    = 64;

constexpr size_type LOW_MASK     = 0xFFFFFFFFULL;

constexpr size_type KILOBYTE     = 1024;
constexpr size_type MEGABYTE     = KILOBYTE * KILOBYTE;
constexpr size_type GIGABYTE     = KILOBYTE * MEGABYTE;

constexpr double KILOBYTE_DOUBLE = 1024.0;
constexpr double MILLION_DOUBLE  = 1000000.0;
constexpr double MEGABYTE_DOUBLE = KILOBYTE_DOUBLE * KILOBYTE_DOUBLE;
constexpr double GIGABYTE_DOUBLE = KILOBYTE_DOUBLE * MEGABYTE_DOUBLE;

constexpr size_type MILLION      = 1000000;
constexpr size_type BILLION      = 1000 * MILLION;

//------------------------------------------------------------------------------

constexpr byte_type IMPL_TERMINATOR   = 0;
constexpr byte_type DATA_TERMINATOR   = 1;
constexpr byte_type STRING_TERMINATOR = 3;

bool is_terminator(byte_type c);

//------------------------------------------------------------------------------

typedef sdsl::int_vector<0>        text_type;
typedef sdsl::int_vector_buffer<0> text_buffer_type;

//------------------------------------------------------------------------------

typedef sdsl::sd_vector<>       vector_type;
typedef std::vector<byte_type>  string_type;

string_type to_string_type(const char* s);

//------------------------------------------------------------------------------

typedef std::pair<byte_type, size_type> run_type;

//------------------------------------------------------------------------------


enum sample_genre : byte_type
{
    NOT =   0,
    START = 1,
    END =   2,
    START_END = 3 // START & AND
};

//------------------------------------------------------------------------------


constexpr size_type alphabet_max_size() { return std::numeric_limits<byte_type>::max() + 1; }
constexpr size_type invalid_value() { return std::numeric_limits<size_type>::max(); }

/*
  range_type stores a closed range [starts, second]. Empty ranges are indicated by
  starts > second. The emptiness check uses +1 to handle the common special case
  [0, -1].
*/

typedef std::pair<size_type, size_type> range_type;

struct Range
{
    static size_type length(range_type range)
    {
        return range.second + 1 - range.first;
    }
    
    static bool empty(range_type range)
    {
        return (range.first + 1 > range.second + 1);
    }
    
    static bool empty(size_type sp, size_type ep)
    {
        return (sp + 1 > ep + 1);
    }
    
    static size_type bound(size_type value, range_type bounds)
    {
        return bound(value, bounds.first, bounds.second);
    }
    
    static size_type bound(size_type value, size_type low, size_type high)
    {
        return std::max(std::min(value, high), low);
    }
    
    // FIXME Change to constexpr when moving to C++14.
    static range_type empty_range()
    {
        return range_type(1, 0);
    }
    
    /*
      Partition the range approximately evenly between the blocks. The actual number of
      blocks will not be greater than the length of the range.
    */
    static std::vector<range_type> partition(range_type range, size_type blocks);
    
    /*
      Range to which the value belong
    */
    static size_type bin(size_type value, const std::vector<range_type>& ranges);
};

template<class A, class B>
std::ostream& operator<<(std::ostream& out, const std::pair<A, B>& data)
{
    return out << "(" << data.first << ", " << data.second << ")";
}

template<class A>
std::ostream& operator<<(std::ostream& out, const std::vector<A>& data)
{
    out << "{ ";
    for(const A& element : data) { out << element << " "; }
    out << "}";
    return out;
}

//------------------------------------------------------------------------------

template<class IntegerType>
    size_type
    bit_length(IntegerType val)
    {
    return sdsl::bits::hi(val) + 1;
}

//------------------------------------------------------------------------------

/*
  Global verbosity setting for index construction. Used in conditions of type
  if(Verbosity::level >= Verbosity::THRESHOLD). While the level can be set directly,
  Verbosity::set() does a few sanity checks.

  SILENT    no status information
  BASIC     basic statistics on the input and the final index
  EXTENDED  intermediate statistics for each batch
  FULL      further details of each batch
*/

struct Verbosity
{
    static size_type level;
    
    static void set(size_type new_level);
    static std::string levelName();
    
    constexpr static size_type SILENT   = 0;
    constexpr static size_type BASIC    = 1;
    constexpr static size_type EXTENDED = 2;
    constexpr static size_type DEFAULT  = 3;
    constexpr static size_type FULL     = 3;
};

//------------------------------------------------------------------------------

/*
Thomas Wang's integer hash function. In many implementations, std::hash
is identity function for integers, which leads to performance issues.
*/

inline size_type
wang_hash_64(size_type key)
{
    key = (~key) + (key << 21); // key = (key << 21) - key - 1;
    key = key ^ (key >> 24);
    key = (key + (key << 3)) + (key << 8); // key * 265
    key = key ^ (key >> 14);
    key = (key + (key << 2)) + (key << 4); // key * 21
    key = key ^ (key >> 28);
    key = key + (key << 31);
    return key;
}

//------------------------------------------------------------------------------

inline double
inMegabytes(size_type bytes)
{
    return bytes / MEGABYTE_DOUBLE;
}

inline double
inGigabytes(size_type bytes)
{
    return bytes / GIGABYTE_DOUBLE;
}

inline double
inBPC(size_type bytes, size_type size)
{
    return (8.0 * bytes) / size;
}

inline double
inMicroseconds(double seconds)
{
    return seconds * MILLION_DOUBLE;
}

constexpr size_type DEFAULT_INDENT = 18;

void printHeader(const std::string& header, size_type indent = DEFAULT_INDENT);
void printTime(const std::string& header, size_type queries, double seconds, size_type indent = DEFAULT_INDENT);
void printTimeLength(const std::string& header, size_type queries, size_type total_length, double seconds, size_type indent = DEFAULT_INDENT);

//------------------------------------------------------------------------------

double readTimer();       // Seconds from an arbitrary time point.
size_type memoryUsage();  // Peak memory usage in bytes.

//------------------------------------------------------------------------------

/*
Temporary file names have the pattern "prefix_hostname_pid_counter", where
- prefix is given as an argument to getName();
- hostname is the name of the host;
- pid is the process id; and
- counter is a running counter starting from 0.

The generated names are stored until the file is deleted with remove(). All
remaining temporary files are deleted when the program exits (normally or
with std::exit()).

TempFile is not thread-safe!
*/

namespace TempFile
{
    extern const std::string DEFAULT_TEMP_DIR;
    extern std::string temp_dir;
    
    void setDirectory(const std::string& directory);
    std::string getName(const std::string& name_part);
    void remove(std::string& filename);  // Also clears the filename.
}

// Returns the total length of the rows, excluding line ends.
// If 'rows' is nonempty, appends to it.
size_type readRows(const std::string& filename, std::vector<std::string>& rows, bool skip_empty_rows);

size_type fileSize(std::ifstream& file);
size_type fileSize(std::ofstream& file);

//------------------------------------------------------------------------------

/*
  parallelQuickSort() uses less working space than parallelMergeSort(). Calling omp_set_nested(1)
  improves the speed of parallelQuickSort().
*/

template<class Iterator, class Comparator>
void
parallelQuickSort(Iterator first, Iterator last, const Comparator& comp)
{
    #ifdef __GLIBCXX__
        int nested = omp_get_nested();
    omp_set_nested(1);
    __gnu_parallel::sort(first, last, comp, __gnu_parallel::balanced_quicksort_tag());
    omp_set_nested(nested);
    #else
        std::sort(first, last, comp);
    #endif
}

template<class Iterator>
void
parallelQuickSort(Iterator first, Iterator last)
{
    #ifdef __GLIBCXX__
        int nested = omp_get_nested();
    omp_set_nested(1);
    __gnu_parallel::sort(first, last, __gnu_parallel::balanced_quicksort_tag());
    omp_set_nested(nested);
    #else
        std::sort(first, last);
    #endif
}

template<class Iterator, class Comparator>
void
parallelMergeSort(Iterator first, Iterator last, const Comparator& comp)
{
    #ifdef __GLIBCXX__
        __gnu_parallel::sort(first, last, comp, __gnu_parallel::multiway_mergesort_tag());
    #else
        std::sort(first, last, comp);
    #endif
}

template<class Iterator>
void
parallelMergeSort(Iterator first, Iterator last)
{
    #ifdef __GLIBCXX__
        __gnu_parallel::sort(first, last, __gnu_parallel::multiway_mergesort_tag());
    #else
        std::sort(first, last);
    #endif
}

template<class Iterator, class Comparator>
void
sequentialSort(Iterator first, Iterator last, const Comparator& comp)
{
    #ifdef __GLIBCXX__
        __gnu_parallel::sort(first, last, comp, __gnu_parallel::sequential_tag());
    #else
        std::sort(first, last, comp);
    #endif
}

template<class Iterator>
void
sequentialSort(Iterator first, Iterator last)
{
    #ifdef __GLIBCXX__
        __gnu_parallel::sort(first, last, __gnu_parallel::sequential_tag());
    #else
        std::sort(first, last);
    #endif
}

template<class Element>
void
removeDuplicates(std::vector<Element>& vec, bool parallel)
{
    if(parallel) { parallelQuickSort(vec.begin(), vec.end()); }
    else         { sequentialSort(vec.begin(), vec.end()); }
    vec.resize(std::unique(vec.begin(), vec.end()) - vec.begin());
}

//------------------------------------------------------------------------------

template<typename value_type>
void
concatenate_file(std::ofstream& dest, std::string& source)
{
    // dest.setf(std::ios::app);
    std::ifstream in(source, std::ios::binary);
    
    std::array<value_type, 4096> tmp;
    while (true)
    {
        in.read(reinterpret_cast<char*>(tmp.data()), 4096 * sizeof(value_type));
        std::size_t cnt = in.gcount() / sizeof(value_type);
        if (cnt == 0) break;
        dest.write(reinterpret_cast<char*>(tmp.data()), cnt * sizeof(value_type));
    }
    
    in.close();
}

} // namespace rimerge

#endif // rimerge_utils_hpp
