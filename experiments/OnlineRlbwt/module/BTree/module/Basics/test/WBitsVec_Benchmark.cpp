#include "BitsUtil.hpp"
#include <stdint.h>


#define BMTEST_BitsUtil
#ifdef BMTEST_BitsUtil
#include <iostream>
#include <iomanip>
#include "cmdline.h"
#include "WBitsVec.hpp"
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

  WBitsVec wbVec(64, num);
  wbVec.resize(num);

  std::cout << "num=" << num << ", rep=" << rep << ", jump=" << jump << ", val=" << val << ": microseconds per read/write for each [w]" << std::endl;

  std::cout << "sequential read (itr): " << std::endl;
  for (uint8_t w = 6; w <= 58; w += 7) {
    w += dummy;
    wbVec.resize(0);
    wbVec.convert(w);
    const uint64_t valMask = bits::UINTW_MAX(w);
    uint64_t checksum0 = 0;
    for (uint64_t i = 0; i < num; ++i) {
      checksum0 += i & valMask;
    }
    for (uint64_t i = 0; i < num; ++i) { // write some values for checksum
      wbVec.write(i & valMask, i);
    }

    double time = 0;
    for (uint64_t r = 0; r < rep; ++r) {
      uint64_t checksum = 0;
      auto t1 = std::chrono::high_resolution_clock::now();
      auto itr = wbVec.begin();
      for (uint64_t i = 0; i < num; ++i, ++itr) {
        checksum += itr.read();
      }
      auto t2 = std::chrono::high_resolution_clock::now();
      time += std::chrono::duration_cast<std::chrono::microseconds>(t2 - t1).count();
      if (checksum0 != checksum) {
        std::cout << std::endl << "error: checksum = " << checksum << " should be " << checksum0 << std::endl;
      }
    }
    std::cout << "[" << static_cast<uint64_t>(w) << "] " << std::fixed << std::setprecision(8) << (time / (rep * num)) << "; ";
  }
  std::cout << std::endl;

  std::cout << "sequential read: " << std::endl;
  for (uint8_t w = 6; w <= 58; w += 7) {
    w += dummy;
    wbVec.resize(0);
    wbVec.convert(w);
    const uint64_t valMask = bits::UINTW_MAX(w);
    uint64_t checksum0 = 0;
    for (uint64_t i = 0; i < num; ++i) {
      checksum0 += i & valMask;
    }
    for (uint64_t i = 0; i < num; ++i) { // write some values for checksum
      wbVec.write(i & valMask, i);
    }

    double time = 0;
    for (uint64_t r = 0; r < rep; ++r) {
      uint64_t checksum = 0;
      auto t1 = std::chrono::high_resolution_clock::now();
      for (uint64_t i = 0; i < num; ++i) {
        checksum += wbVec.read(i);
      }
      auto t2 = std::chrono::high_resolution_clock::now();
      time += std::chrono::duration_cast<std::chrono::microseconds>(t2 - t1).count();
      if (checksum0 != checksum) {
        std::cout << std::endl << "error: checksum = " << checksum << " should be " << checksum0 << std::endl;
      }
    }
    std::cout << "[" << static_cast<uint64_t>(w) << "] " << std::fixed << std::setprecision(8) << (time / (rep * num)) << "; ";
  }
  std::cout << std::endl;

  std::cout << "sequential read with [] operator: " << std::endl;
  for (uint8_t w = 6; w <= 58; w += 7) {
    w += dummy;
    wbVec.resize(0);
    wbVec.convert(w);
    const uint64_t valMask = bits::UINTW_MAX(w);
    uint64_t checksum0 = 0;
    for (uint64_t i = 0; i < num; ++i) {
      checksum0 += i & valMask;
    }
    for (uint64_t i = 0; i < num; ++i) { // write some values for checksum
      wbVec[i] = i & valMask;
    }

    double time = 0;
    for (uint64_t r = 0; r < rep; ++r) {
      uint64_t checksum = 0;
      auto t1 = std::chrono::high_resolution_clock::now();
      for (uint64_t i = 0; i < num; ++i) {
        checksum += wbVec[i];
      }
      auto t2 = std::chrono::high_resolution_clock::now();
      time += std::chrono::duration_cast<std::chrono::microseconds>(t2 - t1).count();
      if (checksum0 != checksum) {
        std::cout << std::endl << "error: checksum = " << checksum << " should be " << checksum0 << std::endl;
      }
    }
    std::cout << "[" << static_cast<uint64_t>(w) << "] " << std::fixed << std::setprecision(8) << (time / (rep * num)) << "; ";
  }
  std::cout << std::endl;




  std::cout << "sequential write (itr): " << std::endl;
  for (uint8_t w = 6; w <= 58; w += 7) {
    w += dummy;
    wbVec.resize(0);
    wbVec.convert(w);

    double time = 0;
    for (uint64_t r = 0; r < rep; ++r) {
      auto t1 = std::chrono::high_resolution_clock::now();
      auto itr = wbVec.begin();
      for (uint64_t i = 0; i < num; ++i, ++itr) {
        itr.write(val);
      }
      auto t2 = std::chrono::high_resolution_clock::now();
      time += std::chrono::duration_cast<std::chrono::microseconds>(t2 - t1).count();
    }
    std::cout << "[" << static_cast<uint64_t>(w) << "] " << std::fixed << std::setprecision(8) << (time / (rep * num)) << "; ";
  }
  std::cout << std::endl;

  std::cout << "sequential write: " << std::endl;
  for (uint8_t w = 6; w <= 58; w += 7) {
    w += dummy;
    wbVec.resize(0);
    wbVec.convert(w);

    double time = 0;
    for (uint64_t r = 0; r < rep; ++r) {
      auto t1 = std::chrono::high_resolution_clock::now();
      for (uint64_t i = 0; i < num; ++i) {
        wbVec.write(val, i);
      }
      auto t2 = std::chrono::high_resolution_clock::now();
      time += std::chrono::duration_cast<std::chrono::microseconds>(t2 - t1).count();
    }
    std::cout << "[" << static_cast<uint64_t>(w) << "] " << std::fixed << std::setprecision(8) << (time / (rep * num)) << "; ";
  }
  std::cout << std::endl;

  std::cout << "sequential write with [] operator: " << std::endl;
  for (uint8_t w = 6; w <= 58; w += 7) {
    w += dummy;
    wbVec.resize(0);
    wbVec.convert(w);

    double time = 0;
    for (uint64_t r = 0; r < rep; ++r) {
      auto t1 = std::chrono::high_resolution_clock::now();
      for (uint64_t i = 0; i < num; ++i) {
        wbVec[i] = val;
      }
      auto t2 = std::chrono::high_resolution_clock::now();
      time += std::chrono::duration_cast<std::chrono::microseconds>(t2 - t1).count();
    }
    std::cout << "[" << static_cast<uint64_t>(w) << "] " << std::fixed << std::setprecision(8) << (time / (rep * num)) << "; ";
  }
  std::cout << std::endl;




  std::cout << "random read: " << std::endl;
  for (uint8_t w = 6; w <= 58; w += 7) {
    w += dummy;
    wbVec.resize(0);
    wbVec.convert(w);
    const uint64_t valMask = bits::UINTW_MAX(w);
    for (uint64_t i = 0; i < num; ++i) { // write some values for checksum
      wbVec.write(i & valMask, i);
    }

    double time = 0;
    for (uint64_t r = 0; r < rep; ++r) {
      uint64_t checksum = 0;
      auto t1 = std::chrono::high_resolution_clock::now();
      for (uint64_t i = 0, pos = jump & idxMask; i < num; ++i, pos = (pos + jump) & idxMask) {
        checksum += wbVec.read(pos);
      }
      auto t2 = std::chrono::high_resolution_clock::now();
      time += std::chrono::duration_cast<std::chrono::microseconds>(t2 - t1).count();
      if (checksum == 0) {
        std::cout << std::endl << "error?: checksum = " << checksum << std::endl;
      }
    }
    std::cout << "[" << (uint64_t)w << "] " << std::fixed << std::setprecision(8) << (time / (rep * num)) << "; ";
  }
  std::cout << std::endl;

  std::cout << "random read with [] operator: " << std::endl;
  for (uint8_t w = 6; w <= 58; w += 7) {
    w += dummy;
    wbVec.resize(0);
    wbVec.convert(w);
    const uint64_t valMask = bits::UINTW_MAX(w);
    for (uint64_t i = 0; i < num; ++i) { // write some values for checksum
      wbVec[i] = i & valMask;
    }

    double time = 0;
    for (uint64_t r = 0; r < rep; ++r) {
      uint64_t checksum = 0;
      auto t1 = std::chrono::high_resolution_clock::now();
      for (uint64_t i = 0, pos = jump & idxMask; i < num; ++i, pos = (pos + jump) & idxMask) {
        checksum += wbVec[pos];
      }
      auto t2 = std::chrono::high_resolution_clock::now();
      time += std::chrono::duration_cast<std::chrono::microseconds>(t2 - t1).count();
      if (checksum == 0) {
        std::cout << std::endl << "error?: checksum = " << checksum << std::endl;
      }
    }
    std::cout << "[" << (uint64_t)w << "] " << std::fixed << std::setprecision(8) << (time / (rep * num)) << "; ";
  }
  std::cout << std::endl;




  std::cout << "random write: " << std::endl;
  for (uint8_t w = 6; w <= 58; w += 7) {
    w += dummy;
    wbVec.resize(0);
    wbVec.convert(w);

    double time = 0;
    for (uint64_t r = 0; r < rep; ++r) {
      auto t1 = std::chrono::high_resolution_clock::now();
      for (uint64_t i = 0, pos = jump & idxMask; i < num; ++i, pos = (pos + jump) & idxMask) {
        wbVec.write(val, pos);
      }
      auto t2 = std::chrono::high_resolution_clock::now();
      time += std::chrono::duration_cast<std::chrono::microseconds>(t2 - t1).count();
    }
    std::cout << "[" << (uint64_t)w << "] " << std::fixed << std::setprecision(8) << (time / (rep * num)) << "; ";
  }
  std::cout << std::endl;

  std::cout << "random write with [] operator: " << std::endl;
  for (uint8_t w = 6; w <= 58; w += 7) {
    w += dummy;
    wbVec.resize(0);
    wbVec.convert(w);

    double time = 0;
    for (uint64_t r = 0; r < rep; ++r) {
      auto t1 = std::chrono::high_resolution_clock::now();
      for (uint64_t i = 0, pos = jump & idxMask; i < num; ++i, pos = (pos + jump) & idxMask) {
        wbVec[pos] = val;
      }
      auto t2 = std::chrono::high_resolution_clock::now();
      time += std::chrono::duration_cast<std::chrono::microseconds>(t2 - t1).count();
    }
    std::cout << "[" << (uint64_t)w << "] " << std::fixed << std::setprecision(8) << (time / (rep * num)) << "; ";
  }
  std::cout << std::endl;
}
#endif
