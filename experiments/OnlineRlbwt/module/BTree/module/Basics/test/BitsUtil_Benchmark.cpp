#include <stdint.h>
#include <iostream>
#include <iomanip>
#include <chrono>
#include "cmdline.h"
#include "BitsUtil.hpp"
#include "MemUtil.hpp"

using namespace itmmti;

int main(int argc, char *argv[])
{
  cmdline::parser parser;
  parser.add<uint64_t>("num",'n', "num of unsigned integers", false, 8388608);
  parser.add<uint64_t>("rep",'r', "num of iteration to compute average time", false, 2);
  // parser.add<uint8_t>("width", 'w', "bit width for packed vector", false, 8, cmdline::range(1, 64));
  parser.add<uint64_t>("jump", 'j', "jump for random access", false, 38201); // 38201 is a prime number
  parser.add<uint64_t>("val", 'v', "value to write (should fit in w bits)", false, 1);
  parser.add<uint8_t>("dummy", 0, "dummy argument (do not input this)", false, 0);
  parser.add("help", 0, "print help");

  parser.parse_check(argc, argv);
  const uint64_t num = parser.get<uint64_t>("num");
  const uint64_t rep = parser.get<uint64_t>("rep");
  // const uint8_t w = parser.get<uint8_t>("width");
  const uint64_t jump = parser.get<uint64_t>("jump");
  uint64_t val = parser.get<uint64_t>("val");
  const uint8_t dummy = parser.get<uint8_t>("dummy");

  const uint64_t idxMask = bits::UINTW_MAX(bits::bitSize(num-1));
  // std::cout << "jump: " << jump << std::endl;

  assert(num > 2);

  std::cout << "num=" << num << ", rep=" << rep << ", jump=" << jump << ", val=" << val << std::endl;

  uint64_t * array0 = memutil::malloc_AbortOnFail<uint64_t>(num);
  uint64_t * array1 = memutil::malloc_AbortOnFail<uint64_t>(num);
  uint64_t checksum0 = 0;
  for (uint64_t i = 0; i < num; ++i) {
    array0[i] = i;
    checksum0 += i;
  }

  {
    std::cout << "cpBits: microseconds to copy " << num << " uint64_t values" << std::endl;
    double time = 0;
    for (uint64_t r = 0; r < rep; ++r) {
      auto t1 = std::chrono::high_resolution_clock::now();
      bits::cpBits_SameOffs(array0, 0, array1, 0, num * 64);
      auto t2 = std::chrono::high_resolution_clock::now();
      time += std::chrono::duration_cast<std::chrono::microseconds>(t2 - t1).count();
      uint64_t checksum = 0;
      for (uint64_t i = 0; i < num; ++i) {
        checksum += array1[i];
      }
      if (checksum0 != checksum) {
        std::cout << std::endl << "error: checksum = " << checksum << " should be " << checksum0 << std::endl;
      }
    }
    std::cout << std::fixed << std::setprecision(8) << (time / rep) << std::endl;
  }

  {
    std::cout << "cpBytes: microseconds to copy " << num << " uint64_t values" << std::endl;
    double time = 0;
    for (uint64_t r = 0; r < rep; ++r) {
      auto t1 = std::chrono::high_resolution_clock::now();
      bits::cpBytes(array0, array1, num * 8);
      auto t2 = std::chrono::high_resolution_clock::now();
      time += std::chrono::duration_cast<std::chrono::microseconds>(t2 - t1).count();
      uint64_t checksum = 0;
      for (uint64_t i = 0; i < num; ++i) {
        checksum += array1[i];
      }
      if (checksum0 != checksum) {
        std::cout << std::endl << "error: checksum = " << checksum << " should be " << checksum0 << std::endl;
      }
    }
    std::cout << std::fixed << std::setprecision(8) << (time / rep) << std::endl;
  }

  {
    std::cout << "mvBitsLR: microseconds to copy " << num << " uint64_t values" << std::endl;
    double time = 0;
    for (uint64_t r = 0; r < rep; ++r) {
      auto t1 = std::chrono::high_resolution_clock::now();
      bits::mvBitsLR(array0, 0, array1, 64, (num - 1) * 64);
      auto t2 = std::chrono::high_resolution_clock::now();
      time += std::chrono::duration_cast<std::chrono::microseconds>(t2 - t1).count();
      uint64_t checksum = 0;
      for (uint64_t i = 1; i < num; ++i) {
        checksum += array1[i];
      }
      checksum += array0[num - 1];
      if (checksum0 != checksum) {
        std::cout << std::endl << "error: checksum = " << checksum << " should be " << checksum0 << std::endl;
      }
    }
    std::cout << std::fixed << std::setprecision(8) << (time / rep) << std::endl;
  }

  {
    std::cout << "mvBitsRL: microseconds to copy " << num << " uint64_t values" << std::endl;
    double time = 0;
    for (uint64_t r = 0; r < rep; ++r) {
      bits::cpBytes(array0, array1, num * 8);
      auto t1 = std::chrono::high_resolution_clock::now();
      bits::mvBitsRL(array1, 0, array1, 64, (num - 1) * 64);
      auto t2 = std::chrono::high_resolution_clock::now();
      time += std::chrono::duration_cast<std::chrono::microseconds>(t2 - t1).count();
      uint64_t checksum = 0;
      for (uint64_t i = 1; i < num; ++i) {
        checksum += array1[i];
      }
      checksum += num - 1;
      if (checksum0 != checksum) {
        std::cout << std::endl << "error: checksum = " << checksum << " should be " << checksum0 << std::endl;
      }
    }
    std::cout << std::fixed << std::setprecision(8) << (time / rep) << std::endl;
  }

  {
    std::cout << "check mvBits: " << num << " uint64_t values" << std::endl;
    double time = 0;
    for (uint64_t r = 0; r < rep; ++r) {
      bits::cpBytes(array0, array1, num * 8);
      auto t1 = std::chrono::high_resolution_clock::now();
      bits::mvBits(array1, 64, array1, 63, (num - 1) * 64);
      bits::mvBits(array1, 63, array1, 61, (num - 1) * 64);
      bits::mvBits(array1, 61, array1, 58, (num - 1) * 64);
      bits::mvBits(array1, 58, array1, 54, (num - 1) * 64);
      bits::mvBits(array1, 54, array1, 49, (num - 1) * 64);
      bits::mvBits(array1, 49, array1, 43, (num - 1) * 64);
      bits::mvBits(array1, 43, array1, 36, (num - 1) * 64);
      bits::mvBits(array1, 36, array1, 64, (num - 1) * 64);
      auto t2 = std::chrono::high_resolution_clock::now();
      time += std::chrono::duration_cast<std::chrono::microseconds>(t2 - t1).count();
      uint64_t checksum = 0;
      for (uint64_t i = 1; i < num; ++i) {
        checksum += array1[i];
      }
      checksum += array0[0];
      if (checksum0 != checksum) {
        std::cout << std::endl << "error: checksum = " << checksum << " should be " << checksum0 << std::endl;
      }
    }
    std::cout << std::fixed << std::setprecision(8) << (time / rep) << std::endl;
  }

  memutil::safefree(array0);
  memutil::safefree(array1);
}
