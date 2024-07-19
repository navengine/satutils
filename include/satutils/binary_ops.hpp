#ifndef NAVENGINE_INCLUDE_SATUTILS_BINARY_OPS_HPP
#define NAVENGINE_INCLUDE_SATUTILS_BINARY_OPS_HPP

#include <cstdint>
#include <array>
#include <cassert>

template<bool LsbIsZero = true>
bool BitVal(const uint32_t& num, const uint8_t pos)
{
  assert(!(pos > 31));

  if constexpr (LsbIsZero) {
    return (num >> pos) & 1;
  } else {
    return num & (0x80000000 >> pos);
  }
}

template<bool LsbIsZero = true>
void BitSet(uint32_t& num, const uint8_t pos)
{
  assert(!(pos > 31));
  if constexpr (LsbIsZero) {
    num |= 1 << pos;
  } else {
    num |= (0x80000000 >> pos);
  }
}

void BitUnset(uint32_t& num, const uint8_t pos);
void BitToggle(uint32_t& num, const uint8_t pos);

template<bool LsbIsZero = true>
void BitEqu(uint16_t& num, const uint8_t pos, bool val)
{
  assert(!(pos > 15));
  if constexpr (LsbIsZero) {
    num ^= (-(uint16_t)val ^ num) & (1 << pos);
  } else {
    num ^= (-(uint16_t)val ^ num) & (0x8000 >> pos);
  }
}

template<bool LsbIsZero = true>
void BitEqu(uint32_t& num, const uint8_t pos, bool val)
{
  assert(!(pos > 31));
  if constexpr (LsbIsZero) {
    num ^= (-(uint32_t)val ^ num) & (1 << pos);
  } else {
    num ^= (-(uint32_t)val ^ num) & (0x80000000 >> pos);
  }
}

// prints MSB to LSB
void PrintBinary(const uint8_t num);
void PrintBinary(const uint16_t num);
void PrintBinary(const uint32_t num);

template<int Size, bool IsLsbZero = true>
bool MultiXor(const uint32_t& num, const std::array<uint8_t,Size> positions)
{
  bool result = BitVal<IsLsbZero>(num, positions[0]);
  for (uint8_t i = 1; i < Size; i++) {
    result ^= BitVal<IsLsbZero>(num, positions[i]);
  }
  return result;
}

template<int Size, bool IsLsbZero = true>
bool MultiXor(const uint32_t& num, const uint8_t positions[])
{
  bool result = BitVal<IsLsbZero>(num, positions[0]);
  for (uint8_t i = 1; i < Size; i++) {
    result ^= BitVal<IsLsbZero>(num, positions[i]);
  }
  return result;
}

#endif
