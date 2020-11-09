#define BMTEST_BlockVec
#ifdef BMTEST_BlockVec
#include <stdint.h>
#include <iostream>
#include <iomanip>
#include "cmdline.h"
#include "BlockVec.hpp"
#include <chrono>

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

  std::cout << "num=" << num << ", rep=" << rep << ", jump=" << jump << ", val=" << val << ": microseconds per read/write" << std::endl;

  {
    BlockVec<uint64_t, 1024> vec(num);
    vec.resize(num);
    std::cout << "BlockVec: kBlockSize = " << vec.kBlockSize << std::endl;

    {
      std::cout << "sequential read: " << std::endl;
      uint64_t checksum0 = 0;
      for (uint64_t i = 0; i < num; ++i) {
        checksum0 += i + dummy;
      }
      for (uint64_t i = 0; i < num; ++i) { // write some values for checksum
        vec[i] = i + dummy + dummy;
      }

      double time = 0;
      for (uint64_t r = 0; r < rep; ++r) {
        uint64_t checksum = 0;
        auto t1 = std::chrono::high_resolution_clock::now();
        for (uint64_t i = 0; i < num; ++i) {
          checksum += vec[i];
        }
        auto t2 = std::chrono::high_resolution_clock::now();
        time += std::chrono::duration_cast<std::chrono::microseconds>(t2 - t1).count();
        if (checksum0 != checksum) {
          std::cout << std::endl << "error: checksum = " << checksum << " should be " << checksum0 << std::endl;
        }
      }
      std::cout << std::fixed << std::setprecision(8) << (time / (rep * num)) << ";" << std::endl;
    }

    {
      std::cout << "sequential write: " << std::endl;
      double time = 0;
      for (uint64_t r = 0; r < rep; ++r) {
        auto t1 = std::chrono::high_resolution_clock::now();
        for (uint64_t i = 0; i < num; ++i) {
          vec[i] = val;
        }
        auto t2 = std::chrono::high_resolution_clock::now();
        time += std::chrono::duration_cast<std::chrono::microseconds>(t2 - t1).count();
      }
      std::cout << std::fixed << std::setprecision(8) << (time / (rep * num)) << ";" << std::endl;
    }

    {
      std::cout << "random read: " << std::endl;
      uint64_t checksum0 = 0;
      for (uint64_t i = 0, pos = jump & idxMask; i < num; ++i, pos = (pos + jump) & idxMask) {
        checksum0 += pos;
      }
      for (uint64_t i = 0; i < num; ++i) { // write some values for checksum
        vec[i] = i + dummy + dummy;
      }

      double time = 0;
      for (uint64_t r = 0; r < rep; ++r) {
        uint64_t checksum = 0;
        auto t1 = std::chrono::high_resolution_clock::now();
        for (uint64_t i = 0, pos = jump & idxMask; i < num; ++i, pos = (pos + jump) & idxMask) {
          checksum += vec[pos];
        }
        auto t2 = std::chrono::high_resolution_clock::now();
        time += std::chrono::duration_cast<std::chrono::microseconds>(t2 - t1).count();
        if (checksum0 != checksum) {
          std::cout << std::endl << "error?: checksum = " << checksum << std::endl;
        }
      }
      std::cout << std::fixed << std::setprecision(8) << (time / (rep * num)) << ";" << std::endl;
    }

    {
      std::cout << "random write: " << std::endl;
      double time = 0;
      for (uint64_t r = 0; r < rep; ++r) {
        auto t1 = std::chrono::high_resolution_clock::now();
        for (uint64_t i = 0, pos = jump & idxMask; i < num; ++i, pos = (pos + jump) & idxMask) {
          vec[pos] = val;
        }
        auto t2 = std::chrono::high_resolution_clock::now();
        time += std::chrono::duration_cast<std::chrono::microseconds>(t2 - t1).count();
      }
      std::cout << std::fixed << std::setprecision(8) << (time / (rep * num)) << ";" << std::endl;
    }
  }






  {
    BlockVec<uint64_t, 2048> vec(num);
    vec.resize(num);
    std::cout << "BlockVec: kBlockSize = " << vec.kBlockSize << std::endl;

    {
      std::cout << "sequential read: " << std::endl;
      uint64_t checksum0 = 0;
      for (uint64_t i = 0; i < num; ++i) {
        checksum0 += i + dummy;
      }
      for (uint64_t i = 0; i < num; ++i) { // write some values for checksum
        vec[i] = i + dummy + dummy;
      }

      double time = 0;
      for (uint64_t r = 0; r < rep; ++r) {
        uint64_t checksum = 0;
        auto t1 = std::chrono::high_resolution_clock::now();
        for (uint64_t i = 0; i < num; ++i) {
          checksum += vec[i];
        }
        auto t2 = std::chrono::high_resolution_clock::now();
        time += std::chrono::duration_cast<std::chrono::microseconds>(t2 - t1).count();
        if (checksum0 != checksum) {
          std::cout << std::endl << "error: checksum = " << checksum << " should be " << checksum0 << std::endl;
        }
      }
      std::cout << std::fixed << std::setprecision(8) << (time / (rep * num)) << ";" << std::endl;
    }

    {
      std::cout << "sequential write: " << std::endl;
      double time = 0;
      for (uint64_t r = 0; r < rep; ++r) {
        auto t1 = std::chrono::high_resolution_clock::now();
        for (uint64_t i = 0; i < num; ++i) {
          vec[i] = val;
        }
        auto t2 = std::chrono::high_resolution_clock::now();
        time += std::chrono::duration_cast<std::chrono::microseconds>(t2 - t1).count();
      }
      std::cout << std::fixed << std::setprecision(8) << (time / (rep * num)) << ";" << std::endl;
    }

    {
      std::cout << "random read: " << std::endl;
      uint64_t checksum0 = 0;
      for (uint64_t i = 0, pos = jump & idxMask; i < num; ++i, pos = (pos + jump) & idxMask) {
        checksum0 += pos;
      }
      for (uint64_t i = 0; i < num; ++i) { // write some values for checksum
        vec[i] = i + dummy + dummy;
      }

      double time = 0;
      for (uint64_t r = 0; r < rep; ++r) {
        uint64_t checksum = 0;
        auto t1 = std::chrono::high_resolution_clock::now();
        for (uint64_t i = 0, pos = jump & idxMask; i < num; ++i, pos = (pos + jump) & idxMask) {
          checksum += vec[pos];
        }
        auto t2 = std::chrono::high_resolution_clock::now();
        time += std::chrono::duration_cast<std::chrono::microseconds>(t2 - t1).count();
        if (checksum0 != checksum) {
          std::cout << std::endl << "error?: checksum = " << checksum << std::endl;
        }
      }
      std::cout << std::fixed << std::setprecision(8) << (time / (rep * num)) << ";" << std::endl;
    }

    {
      std::cout << "random write: " << std::endl;
      double time = 0;
      for (uint64_t r = 0; r < rep; ++r) {
        auto t1 = std::chrono::high_resolution_clock::now();
        for (uint64_t i = 0, pos = jump & idxMask; i < num; ++i, pos = (pos + jump) & idxMask) {
          vec[pos] = val;
        }
        auto t2 = std::chrono::high_resolution_clock::now();
        time += std::chrono::duration_cast<std::chrono::microseconds>(t2 - t1).count();
      }
      std::cout << std::fixed << std::setprecision(8) << (time / (rep * num)) << ";" << std::endl;
    }
  }






  {
    std::cout << "std::vector:" << std::endl;
    std::vector<uint64_t> vec(num);
    vec.resize(num);

    {
      std::cout << "sequential read: " << std::endl;
      uint64_t checksum0 = 0;
      for (uint64_t i = 0; i < num; ++i) {
        checksum0 += i + dummy;
      }
      for (uint64_t i = 0; i < num; ++i) { // write some values for checksum
        vec[i] = i + dummy + dummy;
      }

      double time = 0;
      for (uint64_t r = 0; r < rep; ++r) {
        uint64_t checksum = 0;
        auto t1 = std::chrono::high_resolution_clock::now();
        for (uint64_t i = 0; i < num; ++i) {
          checksum += vec[i];
        }
        auto t2 = std::chrono::high_resolution_clock::now();
        time += std::chrono::duration_cast<std::chrono::microseconds>(t2 - t1).count();
        if (checksum0 != checksum) {
          std::cout << std::endl << "error: checksum = " << checksum << " should be " << checksum0 << std::endl;
        }
      }
      std::cout << std::fixed << std::setprecision(8) << (time / (rep * num)) << ";" << std::endl;
    }

    {
      std::cout << "sequential write: " << std::endl;
      double time = 0;
      for (uint64_t r = 0; r < rep; ++r) {
        auto t1 = std::chrono::high_resolution_clock::now();
        for (uint64_t i = 0; i < num; ++i) {
          vec[i] = val;
        }
        auto t2 = std::chrono::high_resolution_clock::now();
        time += std::chrono::duration_cast<std::chrono::microseconds>(t2 - t1).count();
      }
      std::cout << std::fixed << std::setprecision(8) << (time / (rep * num)) << ";" << std::endl;
    }

    {
      std::cout << "random read: " << std::endl;
      uint64_t checksum0 = 0;
      for (uint64_t i = 0, pos = jump & idxMask; i < num; ++i, pos = (pos + jump) & idxMask) {
        checksum0 += pos;
      }
      for (uint64_t i = 0; i < num; ++i) { // write some values for checksum
        vec[i] = i + dummy + dummy;
      }

      double time = 0;
      for (uint64_t r = 0; r < rep; ++r) {
        uint64_t checksum = 0;
        auto t1 = std::chrono::high_resolution_clock::now();
        for (uint64_t i = 0, pos = jump & idxMask; i < num; ++i, pos = (pos + jump) & idxMask) {
          checksum += vec[pos];
        }
        auto t2 = std::chrono::high_resolution_clock::now();
        time += std::chrono::duration_cast<std::chrono::microseconds>(t2 - t1).count();
        if (checksum0 != checksum) {
          std::cout << std::endl << "error?: checksum = " << checksum << std::endl;
        }
      }
      std::cout << std::fixed << std::setprecision(8) << (time / (rep * num)) << ";" << std::endl;
    }

    {
      std::cout << "random write: " << std::endl;
      double time = 0;
      for (uint64_t r = 0; r < rep; ++r) {
        auto t1 = std::chrono::high_resolution_clock::now();
        for (uint64_t i = 0, pos = jump & idxMask; i < num; ++i, pos = (pos + jump) & idxMask) {
          vec[pos] = val;
        }
        auto t2 = std::chrono::high_resolution_clock::now();
        time += std::chrono::duration_cast<std::chrono::microseconds>(t2 - t1).count();
      }
      std::cout << std::fixed << std::setprecision(8) << (time / (rep * num)) << ";" << std::endl;
    }
  }






  {
    {
      double time = 0;
      for (uint64_t r = 0; r < rep; ++r) {
        BlockVec<uint64_t, 1024> vec;
        std::cout << "appendBlock: kBlockSize = " << vec.kBlockSize << ", r = " << r << std::endl;
        std::cout << "size  = " << vec.size() << ", capacity = " << vec.capacity() << std::endl;

        auto t1 = std::chrono::high_resolution_clock::now();
        for (uint64_t i = 0; i < num; ++i) {
          vec.resize(i + 1);
          vec[i] = val;
        }
        auto t2 = std::chrono::high_resolution_clock::now();
        time += std::chrono::duration_cast<std::chrono::microseconds>(t2 - t1).count();

        std::cout << "size  = " << vec.size() << ", capacity = " << vec.capacity() << std::endl;
      }
      std::cout << std::fixed << std::setprecision(8) << (time / (rep * num)) << ";" << std::endl;
    }

    {
      double time = 0;
      for (uint64_t r = 0; r < rep; ++r) {
        std::vector<uint64_t> vec;
        std::cout << "vector.push_back: r = " << r << std::endl;
        std::cout << "size  = " << vec.size() << ", capacity = " << vec.capacity() << std::endl;

        auto t1 = std::chrono::high_resolution_clock::now();
        for (uint64_t i = 0; i < num; ++i) {
          vec.push_back(val);
        }
        auto t2 = std::chrono::high_resolution_clock::now();
        time += std::chrono::duration_cast<std::chrono::microseconds>(t2 - t1).count();

        std::cout << "size  = " << vec.size() << ", capacity = " << vec.capacity() << std::endl;
      }
      std::cout << std::fixed << std::setprecision(8) << (time / (rep * num)) << ";" << std::endl;
    }
  }
}
#endif
