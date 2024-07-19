#include <cstdint>
#include <iostream>
#include <cassert>

#include <satutils/binary_ops.hpp>

void BitUnset(uint32_t& num, const uint8_t pos)
{
  assert(!(pos > 31));
  num &= ~(1 << pos);
}

void BitToggle(uint32_t& num, const uint8_t pos)
{
  assert(!(pos > 31));
  num ^= 1 << pos;
}

void PrintBinary(const uint8_t num)
{
  for (uint8_t i = 0; i < 8; i++) {
    std::cout << BitVal(num, 7 - i);
  }
  std::cout << std::endl;
}

void PrintBinary(const uint16_t num)
{
  for (uint8_t i = 0; i < 16; i++) {
    std::cout << BitVal(num, 15 - i);
  }
  std::cout << std::endl;
}

void PrintBinary(const uint32_t num)
{
  for (uint8_t i = 0; i < 32; i++) {
    std::cout << BitVal(num, 31 - i);
  }
  std::cout << std::endl;
}


