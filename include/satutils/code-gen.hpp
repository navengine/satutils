/**
 * *code-gen.hpp*
 *
 * =======  ========================================================================================
 * @file    satutils/code-gen.hpp
 * @brief   Constellation gold code generators.
 * @date    January 2025
 * @author  Daniel Sturdivant <Auburn University GAVLAB>
 * @author  Blake Baker <Auburn University GAVLAB>
 * @ref     1. "IS-GPS-200N", 2022
 *          2. "Global Positioning System: Signals, Measurements, and Performance", 2nd Edition,
 *              2006 - Misra & Enge
 * =======  ========================================================================================
 */

#ifndef SATUTILS_CODE_GEN_HPP
#define SATUTILS_CODE_GEN_HPP

#include <array>
#include <cassert>
#include <cstdint>
#include <navtools/binary-ops.hpp>

namespace satutils {

inline void CodeGenCA(bool sequence[1023], uint8_t prn)
{
  assert(!((prn < 1) || (prn > 32)));
  prn -= 1;

  static constexpr uint8_t g2_out_taps[32][2] = {
      /*01 {2, 6} */ {1, 5},
      /*02 {3, 7} */ {2, 6},
      /*03 {4, 8} */ {3, 7},
      /*04 {5, 9} */ {4, 8},
      /*05 {1, 9} */ {0, 8},
      /*06 {2, 10}*/ {1, 9},
      /*07 {1, 8} */ {0, 7},
      /*08 {2, 9} */ {1, 8},
      /*09 {3, 10}*/ {2, 9},
      /*10 {2, 3} */ {1, 2},
      /*11 {3, 4} */ {2, 3},
      /*12 {5, 6} */ {4, 5},
      /*13 {6, 7} */ {5, 6},
      /*14 {7, 8} */ {6, 7},
      /*15 {8, 9} */ {7, 8},
      /*16 {9, 10}*/ {8, 9},
      /*17 {1, 4} */ {0, 3},
      /*18 {2, 5} */ {1, 4},
      /*19 {3, 6} */ {2, 5},
      /*20 {4, 7} */ {3, 6},
      /*21 {5, 8} */ {4, 7},
      /*22 {6, 9} */ {5, 8},
      /*23 {1, 3} */ {0, 2},
      /*24 {4, 6} */ {3, 5},
      /*25 {5, 7} */ {4, 6},
      /*26 {6, 8} */ {5, 7},
      /*27 {7, 9} */ {6, 8},
      /*28 {8, 10}*/ {7, 9},
      /*29 {1, 6} */ {0, 5},
      /*30 {2, 7} */ {1, 6},
      /*31 {3, 8} */ {2, 7},
      /*32 {4, 9} */ {3, 8}};

  // Linear-feedback shift registers
  uint32_t G1 = 0xFFFFFFFF;
  uint32_t G2 = 0xFFFFFFFF;

  uint8_t taps1[2] = {2, 9};              // 3,10
  uint8_t taps2[6] = {1, 2, 5, 7, 8, 9};  // 2,3,6,8,9,10

  for (std::size_t i = 0; i < 1023; i++) {
    // set value in sequence
    sequence[i] = navtools::GetBit<true>(G1, 9u) ^
                  navtools::GetBit<true>(G2, g2_out_taps[prn][0]) ^
                  navtools::GetBit<true>(G2, g2_out_taps[prn][1]);

    // shift the registers and set first bits
    bool feedback1 = navtools::MultiXor<2, true>(G1, taps1);
    bool feedback2 = navtools::MultiXor<6, true>(G2, taps2);
    G1 <<= 1;
    G2 <<= 1;
    navtools::SetBitTo<true>(G1, 0, feedback1);
    navtools::SetBitTo<true>(G2, 0, feedback2);
  }
};

inline void CodeGenCA(std::array<bool, 1023>& sequence, const uint8_t prn)
{
  CodeGenCA(sequence.data(), prn);
}


inline void CodeGenL5IQ(bool l5i[10230], bool l5q[10230], uint8_t prn)
{
  assert(!((prn < 1) || (prn > 63)));
  prn -= 1;

  static constexpr uint16_t xb_initial_states[63][2] = {
    /*PRN 01*/ {0b0101011100100000, 0b1001011001100000},
    /*PRN 02*/ {0b1100000110101000, 0b0100011110110000},
    /*PRN 03*/ {0b0100000001000000, 0b1111000100011000},
    /*PRN 04*/ {0b1011000100110000, 0b0011101101010000},
    /*PRN 05*/ {0b1110111010111000, 0b0011110110010000},
    /*PRN 06*/ {0b0110011111010000, 0b0101010101001000},
    /*PRN 07*/ {0b1010010011111000, 0b1111110000001000},
    /*PRN 08*/ {0b1011110100100000, 0b0110101101000000},
    /*PRN 09*/ {0b1111100101011000, 0b1011101000011000},
    /*PRN 10*/ {0b0111111011110000, 0b0010010000110000},
    /*PRN 11*/ {0b0000100111010000, 0b0001000000101000},
    /*PRN 12*/ {0b1110011111001000, 0b0101011000101000},
    /*PRN 13*/ {0b0001110011100000, 0b0100110100101000},
    /*PRN 14*/ {0b0100000100111000, 0b1010000111111000},
    /*PRN 15*/ {0b0110101011010000, 0b1011110001111000},
    /*PRN 16*/ {0b0001111001001000, 0b1101001011111000},
    /*PRN 17*/ {0b0100110001111000, 0b1110011001000000},
    /*PRN 18*/ {0b1111000011110000, 0b1011011100100000},
    /*PRN 19*/ {0b1100100011111000, 0b0011001011011000},
    /*PRN 20*/ {0b0110101101101000, 0b1100001110001000}, 
    /*PRN 21*/ {0b0010000001000000, 0b0110110010000000}, 
    /*PRN 22*/ {0b1110111101111000, 0b0010110001110000}, 
    /*PRN 23*/ {0b1000011111110000, 0b1000101111101000}, 
    /*PRN 24*/ {0b1100010110100000, 0b0110111110011000}, 
    /*PRN 25*/ {0b1101001101101000, 0b0100010011011000}, 
    /*PRN 26*/ {0b1010110010110000, 0b0101010111100000}, 
    /*PRN 27*/ {0b0101011011110000, 0b1000011111010000}, 
    /*PRN 28*/ {0b0111101010110000, 0b1111101000010000}, 
    /*PRN 29*/ {0b0101111100001000, 0b0101000100100000}, 
    /*PRN 30*/ {0b1000010110111000, 0b1000001111001000}, 
    /*PRN 31*/ {0b0001010011110000, 0b0101111100101000}, 
    /*PRN 32*/ {0b0000010111001000, 0b1001000101010000}, 
    /*PRN 33*/ {0b1101010000001000, 0b1011001000100000}, 
    /*PRN 34*/ {0b1101111111001000, 0b1111001000100000}, 
    /*PRN 35*/ {0b1111011011100000, 0b0110010110011000}, 
    /*PRN 36*/ {0b1001011001000000, 0b0011110101111000}, 
    /*PRN 37*/ {0b0011010010000000, 0b0010011010001000}, 
    /*PRN 38*/ {0b0101100000110000, 0b1111110011101000},
    /*PRN 39*/ {0b1001001100101000, 0b0101010011111000},
    /*PRN 40*/ {0b1100111001010000, 0b1000110101010000},
    /*PRN 41*/ {0b0111011011001000, 0b0010111100100000},
    /*PRN 42*/ {0b0011101101100000, 0b1011000100000000},
    /*PRN 43*/ {0b0011011111010000, 0b0011001011001000},
    /*PRN 44*/ {0b1001011010001000, 0b1000100101000000},
    /*PRN 45*/ {0b1001010111111000, 0b0000001111110000},
    /*PRN 46*/ {0b0111000111101000, 0b0000000010011000},
    /*PRN 47*/ {0b0000001000100000, 0b0101110011110000},
    /*PRN 48*/ {0b1000101010001000, 0b0001001000111000},
    /*PRN 49*/ {0b0011010001001000, 0b0011110000100000},
    /*PRN 50*/ {0b1000111110001000, 0b0100101011100000},
    /*PRN 51*/ {0b1011100101001000, 0b0010100011111000},
    /*PRN 52*/ {0b0100101011010000, 0b1101110011001000},
    /*PRN 53*/ {0b0000001000010000, 0b0011111101111000},
    /*PRN 54*/ {0b0110001101110000, 0b1100100110111000},
    /*PRN 55*/ {0b0000011001110000, 0b1001001100110000},
    /*PRN 56*/ {0b1110111011110000, 0b0100010011001000},
    /*PRN 57*/ {0b0001000010011000, 0b0000000001011000},
    /*PRN 58*/ {0b0000010100001000, 0b0000001101111000},
    /*PRN 59*/ {0b0100001100001000, 0b0101101101111000},
    /*PRN 60*/ {0b0100101001001000, 0b0100100001101000},
    /*PRN 61*/ {0b0011110011110000, 0b1101100101011000},
    /*PRN 62*/ {0b1011000110001000, 0b1010111000100000},
    /*PRN 63*/ {0b0101111001011000, 0b0010001101001000}
  };

  static constexpr uint8_t xa_taps [4] = {8,9,11,12};
  static constexpr uint8_t xb_taps [8] = {0,2,3,5,6,7,11,12};

  uint16_t XBI = xb_initial_states[prn][0];
  uint16_t XBQ = xb_initial_states[prn][1];
  uint16_t XA = 0xFFFF;

  for (std::size_t j = 0; j < 8190; j++) {
    l5i[j] = navtools::GetBit<false>(XA, 12) ^ navtools::GetBit<false>(XBI, 12);
    l5q[j] = navtools::GetBit<false>(XA, 12) ^ navtools::GetBit<false>(XBQ, 12);

    bool feedback_a = navtools::MultiXor<4,false>(XA,xa_taps);
    bool feedback_bi = navtools::MultiXor<8,false>(XBI,xb_taps);
    bool feedback_bq = navtools::MultiXor<8,false>(XBQ,xb_taps);
    
    XA >>= 1;
    XBI >>= 1;
    XBQ >>= 1;
    navtools::SetBitTo<false>(XA,0,feedback_a);
    navtools::SetBitTo<false>(XBI,0,feedback_bi);
    navtools::SetBitTo<false>(XBQ,0,feedback_bq);
  }
  XA = 0xFFFF;
  for (std::size_t j = 8190; j < 10230; j++) {
    l5i[j] = navtools::GetBit<false>(XA, 12) ^ navtools::GetBit<false>(XBI, 12);
    l5q[j] = navtools::GetBit<false>(XA, 12) ^ navtools::GetBit<false>(XBQ, 12);

    bool feedback_a = navtools::MultiXor<4,false>(XA,xa_taps);
    bool feedback_bi = navtools::MultiXor<8,false>(XBI,xb_taps);
    bool feedback_bq = navtools::MultiXor<8,false>(XBQ,xb_taps);

    XA >>= 1;
    XBI >>= 1;
    XBQ >>= 1;
    navtools::SetBitTo<false>(XA,0,feedback_a);
    navtools::SetBitTo<false>(XBI,0,feedback_bi);
    navtools::SetBitTo<false>(XBQ,0,feedback_bq);
  }
}

inline void CodeGenL5IQ(std::array<bool,10230> l5i, std::array<bool,10230> l5q,
                        uint8_t prn)
{
  CodeGenL5IQ(l5i.data(), l5q.data(), prn);
}


};  // namespace satutils
#endif
