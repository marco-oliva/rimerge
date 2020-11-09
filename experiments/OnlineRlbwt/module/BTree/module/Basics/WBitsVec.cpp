#include "WBitsVec.hpp"


#define TEST_CORR_SIMPLE_
#ifdef TEST_CORR_SIMPLE_
#include <time.h>
#include <iostream>
#include <iomanip>
#include "cmdline.h"

//
// $ g++ -std=c++14 -march=native -O3 -DNDEBUG -W -Wall -Wno-deprecated -c BitsUtil.cpp
// $ g++ -std=c++14 -march=native -O3 -DNDEBUG -W -Wall -Wno-deprecated -o WBitsVec.out WBitsVec.cpp BitsUtil.o
// $ ./WBitsVec.out -n 100
//
int main(int argc, char *argv[])
{
  cmdline::parser parser;
  parser.add<uint64_t>("num",'n', "num of unsigned integers", true);
  // parser.add<uint8_t>("width", 'w', "bit width for packed vector", false, 8, cmdline::range(1, 64));
  parser.add<uint64_t>("jump", 'j', "jump for random access", false, 38201); // 38201 is a prime number
  parser.add<uint64_t>("val", 'v', "value to write (should fit in w bits)", false, 1);
  parser.add<uint8_t>("dummy", 0, "dummy argument (do not input this)", false, 0);
  parser.add("help", 0, "pring help");

  parser.parse_check(argc, argv);
  const uint64_t num = parser.get<uint64_t>("num");
  // const uint8_t w = parser.get<uint8_t>("width");
  // const uint64_t jump = parser.get<uint64_t>("jump");
  // uint64_t val = parser.get<uint64_t>("val");
  // const uint8_t dummy = parser.get<uint8_t>("dummy");

  // const uint64_t idxMask = bits::UINTW_MAX(bits::bitSize(num));
  // std::cout << "jump: " << jump << std::endl;

  assert(num > 2);

  WBitsVec wba[10];

  for (uint8_t i = 0, w = 6; w <= 58; ++i, w += 7) {
    wba[i].convert(w, num);
    wba[i].resize(num);
    std::cout << "calcSpace [w = " << static_cast<uint32_t>(wba[i].getW()) << "]: " << wba[i].calcMemBytes();
    std::cout << ", size / capacity_: " << wba[i].size() << " / " << wba[i].capacity() << std::endl;
  }
  // read
  std::cout << "read/write" << std::endl;
  for (uint8_t i = 0, w = 6; w <= 58; ++i, w += 7) {
    std::cout << static_cast<uint32_t>(i) << ": w = " << static_cast<uint32_t>(wba[i].getW()) << std::endl;
    for (uint64_t j = 0; j < num; ++j) {
      wba[i].write(j & bits::UINTW_MAX(w), j);
    }
    for (uint64_t j = 0; j < num; ++j) {
      std::cout << wba[i].read(j) << ", ";
    }
    std::cout << std::endl;
  }
  // right shift mv
  std::cout << "right 4-shift mv" << std::endl;
  for (uint8_t i = 0, w = 6; w <= 58; ++i, w += 7) {
    std::cout << static_cast<uint32_t>(i) << ": w = " << static_cast<uint32_t>(wba[i].getW()) << std::endl;
    mvWBA_SameW(wba[i].getItrAt(0), wba[i].getItrAt(4), num - 4);
    for (uint64_t j = 0; j < num; ++j) {
      std::cout << wba[i].read(j) << ", ";
    }
    std::cout << std::endl;
  }
  // right shift mv
  std::cout << "left 8-shift mv" << std::endl;
  for (uint8_t i = 0, w = 6; w <= 58; ++i, w += 7) {
    std::cout << static_cast<uint32_t>(i) << ": w = " << static_cast<uint32_t>(wba[i].getW()) << std::endl;
    auto itr1 = wba[i].getItrAt(8);
    auto itr2 = wba[i].getItrAt(0);
    mvWBA_SameW(std::move(itr1), std::move(itr2), num - 8);
    for (uint64_t j = 0; j < num; ++j) {
      std::cout << wba[i].read(j) << ", ";
    }
    std::cout << std::endl;
  }
  // convert decrease 4 bits
  std::cout << "convert (decrease 4 bits)" << std::endl;
  for (uint8_t i = 0, w = 6; w <= 58; ++i, w += 7) {
    std::cout << static_cast<uint32_t>(i) << ": w from " << static_cast<uint32_t>(wba[i].getW());
    wba[i].convert(w-4);
    std::cout << " to " << static_cast<uint32_t>(wba[i].getW()) << std::endl;
    for (uint64_t j = 0; j < num; ++j) {
      std::cout << wba[i].read(j) << ", ";
    }
    std::cout << std::endl;
  }
  // convert decrease 4 bits
  std::cout << "convert (increase 7 bits)" << std::endl;
  for (uint8_t i = 0, w = 6; w <= 58; ++i, w += 7) {
    std::cout << static_cast<uint32_t>(i) << ": w from " << static_cast<uint32_t>(wba[i].getW());
    wba[i].convert(wba[i].getW() + 7);
    std::cout << " to " << static_cast<uint32_t>(wba[i].getW()) << std::endl;
    for (uint64_t j = 0; j < num; ++j) {
      std::cout << wba[i].read(j) << ", ";
    }
    std::cout << std::endl;
  }
  for (uint8_t i = 0, w = 6; w <= 58; ++i, w += 7) {
    std::cout << "calcSpace [w = " << static_cast<uint32_t>(wba[i].getW()) << "]: " << wba[i].calcMemBytes();
    std::cout << ", size / capacity_: " << wba[i].size() << " / " << wba[i].capacity() << std::endl;
  }
  //
  std::cout << "copy between different wbVecs" << std::endl;
  for (uint8_t i = 0, w = 6; w <= 58; ++i, w += 7) {
    if (i == 0) continue;
    std::cout << static_cast<uint32_t>(i) << ": from " << (int)(wba[i-1].getW()) << " to " << (int)(wba[i].getW()) << std::endl;
    mvWBA(wba[i-1].getItrAt(i-1), wba[i].getItrAt(i), 10);
    for (uint64_t j = 0; j < num; ++j) {
      std::cout << wba[i-1].read(j) << ", ";
    }
    std::cout << std::endl;
    for (uint64_t j = 0; j < num; ++j) {
      std::cout << wba[i].read(j) << ", ";
    }
    std::cout << std::endl;
  }
  //
  std::cout << "copy between different wbVecs" << std::endl;
  wba[1].write(bits::UINTW_MAX(wba[1].getW()), num-1);
  for (uint8_t i = 0, w = 6; w <= 58; ++i, w += 7) {
    if (i == 0) continue;
    std::cout << static_cast<uint32_t>(i) << ": from " << (int)(wba[i].getW()) << " to " << (int)(wba[i-1].getW()) << std::endl;
    mvWBA(wba[i].getItrAt(num-10), wba[i-1].getItrAt(num-10), 10);
    for (uint64_t j = 0; j < num; ++j) {
      std::cout << wba[i-1].read(j) << ", ";
    }
    std::cout << std::endl;
    for (uint64_t j = 0; j < num; ++j) {
      std::cout << wba[i].read(j) << ", ";
    }
    std::cout << std::endl;
  }
}
#endif


// #define TEST_PERFORMANCE_
#ifdef TEST_PERFORMANCE_
#include <time.h>
#include <iostream>
#include <iomanip>
#include "cmdline.h"

int main(int argc, char *argv[])
{
  cmdline::parser parser;
  parser.add<uint64_t>("num",'n', "num of unsigned integers", true);
  // parser.add<uint8_t>("width", 'w', "bit width for packed vector", false, 8, cmdline::range(1, 64));
  parser.add<uint64_t>("jump", 'j', "jump for random access", false, 38201); // 38201 is a prime number
  parser.add<uint64_t>("val", 'v', "value to write (should fit in w bits)", false, 1);
  parser.add<uint8_t>("dummy", 0, "dummy argument (do not input this)", false, 0);
  parser.add("help", 0, "print help");

  parser.parse_check(argc, argv);
  const uint64_t num = parser.get<uint64_t>("num");
  // const uint8_t w = parser.get<uint8_t>("width");
  const uint64_t jump = parser.get<uint64_t>("jump");
  uint64_t val = parser.get<uint64_t>("val");
  const uint8_t dummy = parser.get<uint8_t>("dummy");

  const uint64_t idxMask = bits::UINTW_MAX(bits::bitSize(num));
  // std::cout << "jump: " << jump << std::endl;

  assert(num > 2);

  WBitsVec wbVec(64, num);
  wbVec.resize(num);

  clock_t start, end;
  double sec;  

  std::cout << __FILE__ << " sequential read (itr): " << std::endl;
  for (uint8_t w = 6; w <= 58; w += 7) {
    w += dummy;
    wbVec.resize(0);
    wbVec.convert(w);
    const uint64_t valMask = (1ULL << w) - 1;
    uint64_t checksum0 = 0;
    for (uint64_t i = 0; i < num; ++i) {
      checksum0 += i & valMask;
    }
    for (uint64_t i = 0; i < num; ++i) { // write some values for checksum
      wbVec.write(i & valMask, i);
    }
    uint64_t checksum = 0;
    start = clock();
    auto itr = wbVec.begin();
    for (uint64_t i = 0; i < num; ++i, ++itr) {
      checksum += itr.read();
    }
    end = clock();
    sec = (double)(end - start)/CLOCKS_PER_SEC;
    std::cout << "[" << static_cast<uint64_t>(w) << "] " << std::setw(8) << sec << "; ";
    if (checksum0 != checksum) {
      std::cout << std::endl << "error: checksum = " << checksum << " should be " << checksum0 << std::endl;
    }
  }
  std::cout << std::endl;

  std::cout << __FILE__ << " sequential write (itr): " << std::endl;
  for (uint8_t w = 6; w <= 58; w += 7) {
    w += dummy;
    wbVec.resize(0);
    wbVec.convert(w);
    start = clock();
    auto itr = wbVec.begin();
    for (uint64_t i = 0; i < num; ++i, ++itr) {
      itr.write(val);
    }
    end = clock();
    sec = (double)(end - start)/CLOCKS_PER_SEC;
    std::cout << "[" << static_cast<uint64_t>(w) << "] " << std::setw(8) << sec << "; ";
  }
  std::cout << std::endl;


  std::cout << __FILE__ << " sequential read: " << std::endl;
  for (uint8_t w = 6; w <= 58; w += 7) {
    w += dummy;
    wbVec.resize(0);
    wbVec.convert(w);
    const uint64_t valMask = (1ULL << w) - 1;
    uint64_t checksum0 = 0;
    for (uint64_t i = 0; i < num; ++i) {
      checksum0 += i & valMask;
    }
    for (uint64_t i = 0; i < num; ++i) { // write some values for checksum
      wbVec.write(i & valMask, i);
    }
    uint64_t checksum = 0;
    start = clock();
    for (uint64_t i = 0; i < num; ++i) {
      checksum += wbVec.read(i);
    }
    end = clock();
    sec = (double)(end - start)/CLOCKS_PER_SEC;
    std::cout << "[" << static_cast<uint64_t>(w) << "] " << std::setw(8) << sec << "; ";
    if (checksum0 != checksum) {
      std::cout << std::endl << "error: checksum = " << checksum << " should be " << checksum0 << std::endl;
    }
  }
  std::cout << std::endl;

  std::cout << __FILE__ << " sequential write: " << std::endl;
  for (uint8_t w = 6; w <= 58; w += 7) {
    w += dummy;
    wbVec.resize(0);
    wbVec.convert(w);
    start = clock();
    for (uint64_t i = 0; i < num; ++i) {
      wbVec.write(val, i);
    }
    end = clock();
    sec = (double)(end - start)/CLOCKS_PER_SEC;
    std::cout << "[" << static_cast<uint64_t>(w) << "] " << std::setw(8) << sec << "; ";
  }
  std::cout << std::endl;



  std::cout << __FILE__ << " random read: " << std::endl;
  for (uint8_t w = 6; w <= 58; w += 7) {
    w += dummy;
    wbVec.resize(0);
    wbVec.convert(w);
    const uint64_t valMask = (1ULL << w) - 1;
    uint64_t checksum0 = 0;
    for (uint64_t i = 0; i < num; ++i) { // write some values for checksum
      wbVec.write(i & valMask, i);
    }
    for (uint64_t i = 0, pos = jump & idxMask; i < num; ++i, pos = (pos + jump) & idxMask) {
      checksum0 += wbVec.read(pos);
    }
    for (uint64_t i = 0; i < num; ++i) { // write some values for checksum
      wbVec.write(i & valMask, i);
    }
    uint64_t checksum = 0;
    start = clock();
    for (uint64_t i = 0, pos = jump & idxMask; i < num; ++i, pos = (pos + jump) & idxMask) {
      checksum += wbVec.read(pos);
    }
    end = clock();
    sec = (double)(end - start)/CLOCKS_PER_SEC;
    std::cout << "[" << (uint64_t)w << "] " << std::setw(8) << sec << "; ";
    if (checksum0 != checksum) {
      std::cout << std::endl << "error: checksum = " << checksum << " should be " << checksum0 << std::endl;
    }
  }
  std::cout << std::endl;

  std::cout << __FILE__ << " random write: " << std::endl;
  for (uint8_t w = 6; w <= 58; w += 7) {
    w += dummy;
    wbVec.resize(0);
    wbVec.convert(w);
    start = clock();
    for (uint64_t i = 0, pos = jump & idxMask; i < num; ++i, pos = (pos + jump) & idxMask) {
      wbVec.write(val, pos);
    }
    end = clock();
    sec = (double)(end - start)/CLOCKS_PER_SEC;
    std::cout << "[" << (uint64_t)w << "] " << std::setw(8) << sec << "; ";
  }
  std::cout << std::endl;
}
#endif

