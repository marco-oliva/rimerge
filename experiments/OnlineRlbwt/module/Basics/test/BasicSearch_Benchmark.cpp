#define BMTEST_BasicSearch
#ifdef BMTEST_BasicSearch
#include <stdint.h>
#include <iostream>
#include <iomanip>
#include <algorithm>
#include <random>
#include <chrono>
#include "cmdline.h"
#include "BasicSearch.hpp"
#include "WBitsVec.hpp"

using namespace itmmti;

int main(int argc, char *argv[])
{
  cmdline::parser parser;
  parser.add<uint32_t>("seed", 's', "random seed", false, 8);
  parser.add<uint8_t>("dummy", 0, "dummy argument (do not input this)", false, 0);
  parser.add("help", 0, "print help");

  parser.parse_check(argc, argv);
  const uint32_t seed = parser.get<uint32_t>("seed");
  const uint8_t dummy = parser.get<uint8_t>("dummy");

  std::cout << "seed = " << seed << std::endl;

  const uint64_t LS = 0;

  {
    const uint64_t step = 1000;
    const uint64_t sizes[] = {128, 1024, 8192, 65536, 131072, 262144, 524288, 1048576};
    // const uint64_t sizes[] = {100, 127, 128, 129, 1000, 1023, 1024, 1025, 1000000, 1048576};
    // const uint64_t sizes[] = {1000000};
    std::cout << "Baseline on Array vs Mymethod on Array: ";
    std::cout << "LS = " << LS << ", step = " << step << std::endl << "array_sizes = {";
    for (uint32_t i = 0; i < sizeof(sizes) / sizeof(sizes[0]); ++i) {
      std::cout << sizes[i] << ", ";
    }
    std::cout << "}" << std::endl;

    for (uint32_t e = 0; e < sizeof(sizes) / sizeof(sizes[0]); ++e) {
      uint64_t checksum0 = 0;
      uint64_t checksum1 = 0;
      double time0 = 0;
      double time1 = 0;

      std::random_device rnd;
      std::mt19937 mt(seed); // 32bit

      const auto size = sizes[e];
      uint32_t * array = new uint32_t[size];
      array[0] = mt() % (step + 1);
      for (uint64_t i = 1; i < size; ++i) {
        array[i] = array[i - 1] + mt() % (step + 1);
      }

      const auto max = array[size - 1];
      const auto vstep = max/size;

      std::cout << "size = " << size << ", max = " << max << ", vstep = " << vstep << std::endl;
      // for (uint64_t i = 0; i < size; ++i) {
      //   std::cout << array[i] << ", ";
      // }
      // std::cout << std::endl;

      {
        uint32_t key = 0;
        auto t1 = std::chrono::high_resolution_clock::now();
        for (uint32_t k = 0; k < size; ++k, key += vstep) {
          checksum0 += (uint64_t)(std::lower_bound(array, array + size, key) - array);
        }
        auto t2 = std::chrono::high_resolution_clock::now();
        time0 += std::chrono::duration_cast<std::chrono::microseconds>(t2 - t1).count();
      }

      {
        uint32_t key = 0;
        auto t1 = std::chrono::high_resolution_clock::now();
        for (uint32_t k = 0; k < size; ++k, key += vstep) {
          checksum1 += basic_search::partition_idx<LS>
            (
             0, size,
             [=](uint64_t i) { return array[i] >= key; }
             );
        }
        auto t2 = std::chrono::high_resolution_clock::now();
        time1 += std::chrono::duration_cast<std::chrono::microseconds>(t2 - t1).count();
      }

      std::cout << "baseline: " << std::fixed << std::setprecision(8) << (time0 / size) << " micro sec per query." << std::endl;
      std::cout << "mymethod: " << std::fixed << std::setprecision(8) << (time1 / size) << " micro sec per query." << std::endl;
      if (checksum0 != checksum1) {
        std::cout << "ERROR: checksums are different." << std::endl;
      }

      delete[] array;
    }
  }







  { // compliment search
    const uint64_t step = 1000;
    const uint64_t sizes[] = {128, 1024, 8192, 65536, 131072, 262144, 524288, 1048576};
    // const uint64_t sizes[] = {100, 127, 128, 129, 1000, 1023, 1024, 1025, 1000000, 1048576};
    // const uint64_t sizes[] = {1000000};
    std::cout << "Baseline on Array vs Mymethod on Array (compliment search): ";
    std::cout << "LS = " << LS << ", step = " << step << std::endl << "array_sizes = {";
    for (uint32_t i = 0; i < sizeof(sizes) / sizeof(sizes[0]); ++i) {
      std::cout << sizes[i] << ", ";
    }
    std::cout << "}" << std::endl;

    for (uint32_t e = 0; e < sizeof(sizes) / sizeof(sizes[0]); ++e) {
      uint64_t checksum0 = 0;
      uint64_t checksum1 = 0;
      double time0 = 0;
      double time1 = 0;

      std::random_device rnd;
      std::mt19937 mt(seed); // 32bit

      const auto size = sizes[e];
      uint32_t * array = new uint32_t[size];
      array[0] = mt() % (step + 1);
      for (uint64_t i = 1; i < size; ++i) {
        array[i] = array[i - 1] + mt() % (step + 1);
      }

      const auto max = array[size - 1];
      const auto vstep = (size * step - max)/size;

      std::cout << "size = " << size << ", max = " << max << ", vstep = " << vstep << std::endl;
      // for (uint64_t i = 0; i < size; ++i) {
      //   std::cout << array[i] << ", ";
      // }
      // std::cout << std::endl;

      {
        uint32_t key = 0;
        auto t1 = std::chrono::high_resolution_clock::now();
        for (uint32_t k = 0; k < size; ++k, key += vstep) {
          checksum0 += (uint64_t)(std::partition_point(array, array + size,
                                                       [=](const uint32_t & a) { return static_cast<uint64_t>(&a - array + 1) * step - a < key; }) - array);
        }
        auto t2 = std::chrono::high_resolution_clock::now();
        time0 += std::chrono::duration_cast<std::chrono::microseconds>(t2 - t1).count();
      }

      {
        uint32_t key = 0;
        auto t1 = std::chrono::high_resolution_clock::now();
        for (uint32_t k = 0; k < size; ++k, key += vstep) {
          checksum1 += basic_search::partition_idx<LS>
            (
             0, size,
             [=](uint64_t i) { return (i + 1) * step - array[i] >= key; }
             );
        }
        auto t2 = std::chrono::high_resolution_clock::now();
        time1 += std::chrono::duration_cast<std::chrono::microseconds>(t2 - t1).count();
      }

      std::cout << "baseline: " << std::fixed << std::setprecision(8) << (time0 / size) << " micro sec per query." << std::endl;
      std::cout << "mymethod: " << std::fixed << std::setprecision(8) << (time1 / size) << " micro sec per query." << std::endl;
      if (checksum0 != checksum1) {
        std::cout << "ERROR: checksums are different." << std::endl;
      }

      delete[] array;
    }
  }





  { // Seach on WBitsVec
    const uint64_t step = 1000;
    const uint64_t sizes[] = {128, 1024, 8192, 65536, 131072, 262144, 524288, 1048576};
    // const uint64_t sizes[] = {100, 127, 128, 129, 1000, 1023, 1024, 1025, 1000000, 1048576};
    // const uint64_t sizes[] = {1000000};
    std::cout << "Baseline on Array vs Mymethod on WBitsArray: ";
    std::cout << "LS = " << LS << ", step = " << step << std::endl << "array_sizes = {";
    for (uint32_t i = 0; i < sizeof(sizes) / sizeof(sizes[0]); ++i) {
      std::cout << sizes[i] << ", ";
    }
    std::cout << "}" << std::endl;

    for (uint32_t e = 0; e < sizeof(sizes) / sizeof(sizes[0]); ++e) {
      uint64_t checksum0 = 0;
      uint64_t checksum1 = 0;
      double time0 = 0;
      double time1 = 0;

      std::random_device rnd;
      std::mt19937 mt(seed); // 32bit

      const auto size = sizes[e];
      uint32_t * array = new uint32_t[size];
      array[0] = mt() % (step + 1);
      for (uint64_t i = 1; i < size; ++i) {
        array[i] = array[i - 1] + mt() % (step + 1);
      }
      WBitsVec wbv(32, size);
      wbv.resize(size);
      for (uint64_t i = 0; i < size; ++i) {
        wbv.write(array[i], i);
      }

      const auto max = array[size - 1];
      const auto vstep = max/size;

      std::cout << "size = " << size << ", max = " << max << ", vstep = " << vstep << std::endl;
      // for (uint64_t i = 0; i < size; ++i) {
      //   std::cout << array[i] << ", ";
      // }
      // std::cout << std::endl;

      {
        uint32_t key = 0;
        auto t1 = std::chrono::high_resolution_clock::now();
        for (uint32_t k = 0; k < size; ++k, key += vstep) {
          checksum0 += (uint64_t)(std::lower_bound(array, array + size, key) - array);
        }
        auto t2 = std::chrono::high_resolution_clock::now();
        time0 += std::chrono::duration_cast<std::chrono::microseconds>(t2 - t1).count();
      }

      {
        uint32_t key = 0;
        auto t1 = std::chrono::high_resolution_clock::now();
        for (uint32_t k = 0; k < size; ++k, key += vstep) {
          checksum1 += basic_search::partition_idx<LS>
            (
             0, size,
             [&](uint64_t i) { return wbv[i] >= key; }
             );
        }
        auto t2 = std::chrono::high_resolution_clock::now();
        time1 += std::chrono::duration_cast<std::chrono::microseconds>(t2 - t1).count();
      }

      std::cout << "baseline: " << std::fixed << std::setprecision(8) << (time0 / size) << " micro sec per query." << std::endl;
      std::cout << "mymethod: " << std::fixed << std::setprecision(8) << (time1 / size) << " micro sec per query." << std::endl;
      if (checksum0 != checksum1) {
        std::cout << "ERROR: checksums are different." << std::endl;
      }

      delete[] array;
    }
  }
}
#endif
