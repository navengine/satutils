/**
 * *gps-lnav.hpp*
 *
 * =======  ========================================================================================
 * @file    satutils/gps-lnav.hpp
 * @brief   Implementation of GPS L1 C/A navigation message utils.
 * @date    January 2025
 * @author  Daniel Sturdivant <Auburn University GAVLAB>
 * @author  Blake Baker <Auburn University GAVLAB>
 * @ref     1. "Understanding GPS/GNSS Principles and Applications", 3rd Edition, 2017
 *            - Kaplan & Hegarty
 *          2. "Global Positioning System: Signals, Measurements, and Performance", 2nd Edition,
 *              2006 - Misra & Enge
 *          3. "A Software-Defined GPS and Galileo Receiver: A Single-Frequency Approach", 2007
 *            - Borre, Akos, Bertelsen, Rinder, Jensen
 * =======  ========================================================================================
 */

#ifndef SATUTILS_GPS_LNAV_HPP
#define SATUTILS_GPS_LNAV_HPP

#include <exception>
#include <iostream>
#include <navtools/binary-ops.hpp>

#include "satutils/atmosphere.hpp"
#include "satutils/ephemeris.hpp"
#include "satutils/gnss-constants.hpp"

namespace satutils {

template <typename T = double>
class GpsLnav : public KeplerElements<T>, KlobucharElements<T> {
 public:
  GpsLnav<T>() = default;
  GpsLnav<T>(const KeplerElements<T> &eph) {
    SetEphem(eph);
  };
  GpsLnav<T>(const KlobucharElements<T> &klob) {
    SetKlobuchar(klob);
  };
  GpsLnav<T>(const KeplerElements<T> &eph, const KlobucharElements<T> &klob) {
    SetEphem(eph);
    SetKlobuchar(klob);
  };

  /**
   * *=== SetEphem ===*
   * @brief set the ephemeris elements
   */
  void SetEphem(const KeplerElements<T> &eph) {
    this->iode = eph.iode;
    this->iodc = eph.iodc;
    this->toe = eph.toe;
    this->toc = eph.toc;
    this->tgd = eph.tgd;
    this->af2 = eph.af2;
    this->af1 = eph.af1;
    this->af0 = eph.af0;
    this->e = eph.e;
    this->sqrtA = eph.sqrtA;
    this->deltan = eph.deltan;
    this->m0 = eph.m0;
    this->omega0 = eph.omega0;
    this->omega = eph.omega;
    this->omegaDot = eph.omegaDot;
    this->i0 = eph.i0;
    this->iDot = eph.iDot;
    this->cuc = eph.cuc;
    this->cus = eph.cus;
    this->cic = eph.cic;
    this->cis = eph.cis;
    this->crc = eph.crc;
    this->crs = eph.crs;
    this->ura = eph.ura;
    this->health = eph.health;
  };

  /**
   * *=== SetKlobuchar ===*
   * @brief Set the Klobuchar elements
   */
  void SetKlobuchar(const KlobucharElements<T> &klob) {
    this->a0 = klob.a0;
    this->a1 = klob.a1;
    this->a2 = klob.a2;
    this->a3 = klob.a3;
    this->b0 = klob.b0;
    this->b1 = klob.b1;
    this->b2 = klob.b2;
    this->b3 = klob.b3;
  }

  /**
   * *=== SetNextBit ===*
   * @brief Read the next navigation data bit
   * @param bit Navigation data bit parsed by tracking loop
   * @return True|False based on if a subframe was successfully parsed
   */
  bool SetNextBit(const bool &bit) {
    // data bits ordered [-2 -1 0 ... 29]
    try {
      bool subframe_parsed = false;

      // shift prev 32 to the left and insert next bit
      prev_32_bits_ <<= 1;
      navtools::SetBitTo<false>(prev_32_bits_, 31, bit);
      bit_cnt_ += 1;
      // std::bitset<32> tmp(prev_32_bits_);

      // Step 1: find preamble
      if (!preamble_sync_) {
        // Check most recent 8 bits for a preamble
        uint8_t test = static_cast<uint8_t>(prev_32_bits_ & 0x000000FF);
        if ((test == LNAV_PREAMBLE_BITS) || (test == LNAV_INV_PREAMBLE_BITS)) {
          if (bits_since_preamble_ == 301) {
            // If here, an initial preamble has been detected
            bit_cnt_ = 8;
            word_cnt_ = 0;
            bits_since_preamble_ = 0;
          } else if (bits_since_preamble_ == 300) {
            // If here, a second preamble has been successfully detected 300 bits apart!
            preamble_sync_ = true;
            subframe_parsed = ParseSubframe();
            bit_cnt_ = 8;
            word_cnt_ = 0;
            bits_since_preamble_ = 0;
          } else {
            // There is another preamble detection here, log it
            preamble_idx_.push_back(bits_since_preamble_);
          }
        }

        // Check if preamble failed to sync
        if (bits_since_preamble_ == 300) {
          // Sync failed, reset to next detection
          uint16_t bits_to_shift = (preamble_idx_[0] % 30) + 1;
          uint16_t words_to_shift = preamble_idx_[0] / 30;

          // shift down full words
          if (words_to_shift > 0) {
            for (uint16_t i = 0; i < (10 - words_to_shift); i++) {
              subframe[i] = subframe[i + words_to_shift];
            }
          }

          // shift down remaining bits
          if (bits_to_shift > 0) {
            uint32_t curr_bits, next_bits;
            for (uint16_t i = 0; i < (9 - words_to_shift); i++) {
              curr_bits = ((subframe[i] & 0xFFFFFFFC) << (bits_to_shift - 1));
              next_bits = (subframe[i + 1] >> (31 - bits_to_shift));
              subframe[i] = curr_bits | next_bits;
            }
          }

          // update counters
          bits_since_preamble_ -= preamble_idx_[0];
          bit_cnt_ = 31 - bits_to_shift + 8;  // +8 from preamble
          word_cnt_ = 9 - words_to_shift;

          // git rid of used detection and update detection counter
          for (uint16_t &x : preamble_idx_) {
            x -= preamble_idx_[0];
          }
          preamble_idx_.erase(preamble_idx_.begin());
        }

        // increment bit counter
        if (bits_since_preamble_ < 300) bits_since_preamble_++;
      }

      // Step 2: Save words
      if (bit_cnt_ == 30) {
        subframe[word_cnt_] = prev_32_bits_;
        bit_cnt_ = 0;
        word_cnt_ += 1;
      }

      // Step 3: Parse subframe
      if (word_cnt_ == 10) {
        word_cnt_ = 0;
        if (preamble_sync_) {
          subframe_parsed = ParseSubframe();
        }
      }

      return subframe_parsed;
    } catch (std::exception &e) {
      std::cerr << "Error! " << e.what() << "\n";
      return false;
    }
  };

  /**
   * *=== ParseSubframe ===*
   * @brief Attempts to parse subframe
   * @return True|False based on if a subframe was successfully parsed
   */
  bool ParseSubframe() {
    // data bits ordered [-2 -1 0 ... 29]
    try {
      // Step 1: Validate received data bits
      bool D29star, D30star;
      for (uint16_t i = 0; i < 10; i++) {
        D29star = navtools::CheckBit<false>(subframe[i], 0);
        D30star = navtools::CheckBit<false>(subframe[i], 1);

        // check bit polarity
        if (D30star) subframe[i] ^= 0x3FFFFFC0;

        // check parity
        if (!ParityCheck(subframe[i], D29star, D30star)) {
          throw std::runtime_error("Invalid parity check!");
        }
      }

      // Step 2: Get subframe id
      uint8_t sub_id = static_cast<uint8_t>((subframe[1] & 0x00000700) >> 8);

      // Step 3: Parse subframe
      bool subframe_parsed = false;
      switch (sub_id) {
        case 1:
          LoadSubframe1();
          sub1_parsed_ = true;
          subframe_parsed = true;
          // std::cout << "Subframe 1 parsed!\n";
          break;
        case 2:
          LoadSubframe2();
          sub2_parsed_ = true;
          subframe_parsed = true;
          // std::cout << "Subframe 2 parsed!\n";
          break;
        case 3:
          LoadSubframe3();
          sub3_parsed_ = true;
          subframe_parsed = true;
          // std::cout << "Subframe 3 parsed!\n";
          ;
          break;
        case 4:
          subframe_parsed = true;
          // std::cout << "Subframe 4 parsed!\n";
          break;
        case 5:
          subframe_parsed = true;
          // std::cout << "Subframe 5 parsed!\n";
          break;
        default:
          throw std::range_error("GpsLnav::ParseSubframe Invalid subframe ID!");
          break;
      }

      return subframe_parsed;
    } catch (std::exception &e) {
      std::cerr << "Error! " << e.what() << "\n";
      return false;
    }
  };

  /**
   * *=== ParityCheck ===*
   * @brief IS-GPS-200N pg. 139
   * @param gpsword GPS word to evaluate
   * @param D29star
   * @param D30star
   * @return Parity success or failure
   */
  bool ParityCheck(const uint32_t &gpsword, const bool &D29star, const bool &D30star) {
    // Calculate the parity
    uint32_t parity = 0;
    navtools::SetBitTo<false>(
        parity, 26, D29star ^ navtools::MultiXor<14, false>(gpsword, GPS_D25));
    navtools::SetBitTo<false>(
        parity, 27, D30star ^ navtools::MultiXor<14, false>(gpsword, GPS_D26));
    navtools::SetBitTo<false>(
        parity, 28, D29star ^ navtools::MultiXor<14, false>(gpsword, GPS_D27));
    navtools::SetBitTo<false>(
        parity, 29, D30star ^ navtools::MultiXor<14, false>(gpsword, GPS_D28));
    navtools::SetBitTo<false>(
        parity, 30, D30star ^ navtools::MultiXor<15, false>(gpsword, GPS_D29));
    navtools::SetBitTo<false>(
        parity, 31, D29star ^ navtools::MultiXor<13, false>(gpsword, GPS_D30));

    // compare the parity
    if (parity == (gpsword & 0x0000003F)) {
      return true;
    }
    return false;
  };

  /**
   * *=== LoadPreamble ===*
   * @brief Reads words 1 and 2 of each subframe (IS-GPS-200N pg. 92)
   */
  void LoadPreamble() {
    // word 1
    // tlm_message = static_cast<uint16_t>((subframe[0] & 0x003FFF00) >> 8);        // bits 9-22
    // integrity_status_flag = static_cast<bool>((subframe[0] & 0x00000080) >> 7);  // bit 23

    // word 2
    ToW_ = 6.0 * static_cast<T>((subframe[1] & 0x3FFFE000) >> 13);  // bits 1-17
    // alert_flag = static_cast<bool>((subframe[1] & 0x00001000) >> 12);       // bit 18
    // anti_spoof_flag = static_cast<bool>((subframe[1] & 0x00000800) >> 11);  // bit 19
  };

  /**
   * *=== LoadSubframe1 ===*
   * @brief Reads GPS LNAV subframe 1 (IS-GPS-200N pg. 80 & 97)
   */
  void LoadSubframe1() {
    uint32_t tmp1, tmp2;

    // Word 1-2
    LoadPreamble();

    // Word 3
    week_ = static_cast<uint16_t>((subframe[2] & 0x3FF00000) >> 20);  // bits 1-10
    // l2_flag = static_cast<uint8_t>((subframe[2] & 0x000C0000) >> 18);     // bits 11-12
    this->ura = static_cast<uint8_t>((subframe[2] & 0x0003C000) >> 14);    // bits 13-16
    this->health = static_cast<uint8_t>((subframe[2] & 0x00003F00) >> 8);  // bits 17-22

    // Word 7
    tmp1 = (subframe[6] & 0x00003FC0) >> 6;
    this->tgd = navtools::TwosComp(tmp1, 8) * PowerOfTwo<-31, T>();  // bits 17-24

    // Word 8
    tmp1 = (subframe[2] & 0x000000C0) >> 6;           // bits 23-24 (word 3)
    tmp2 = (subframe[7] & 0x3FC00000) >> 22;          // bits 1-8
    this->iodc = static_cast<T>((tmp1 << 8) | tmp2);  //
    this->toc = static_cast<T>((subframe[7] & 0x003FFFC0) >> 6) * PowerOfTwo<4, T>();  // bits 9-24

    // Word 9
    tmp1 = (subframe[8] & 0x000000C0) >> 22;
    tmp2 = (subframe[8] & 0x003FFFC0) >> 6;
    this->af2 = navtools::TwosComp(tmp1, 8) * PowerOfTwo<-55, T>();   // bits 1-8
    this->af1 = navtools::TwosComp(tmp2, 16) * PowerOfTwo<-43, T>();  // bits 9-24

    // word 10
    tmp1 = (subframe[9] & 0x3FFFFF00) >> 8;
    this->af0 = navtools::TwosComp(tmp1, 22) * PowerOfTwo<-31, T>();  // bits 1-22
  };

  /**
   * *=== LoadSubframe2 ===*
   * @brief Reads GPS LNAV subframe 2 (IS-GPS-200N pg. 81 & 105)
   */
  void LoadSubframe2() {
    uint32_t tmp1, tmp2, tmp3;

    // Word 1-2
    LoadPreamble();

    // Word 3
    tmp3 = (subframe[2] & 0x003FFFC0) >> 6;
    this->iode = static_cast<T>((subframe[2] & 0x3FC00000) >> 22);   // bits 1-8
    this->crs = navtools::TwosComp(tmp3, 16) * PowerOfTwo<-5, T>();  // bits 9-24

    // Word 4-5
    tmp3 = (subframe[3] & 0x3FFFC000) >> 14;
    this->deltan = navtools::TwosComp(tmp3, 16) * GPS_PI<T> * PowerOfTwo<-43, T>();  // bits 1-16
    tmp1 = (subframe[3] & 0x00003FC0) >> 6;                                          // bits 17-24
    tmp2 = (subframe[4] & 0x3FFFFFC0) >> 6;                                          // bits 1-24
    tmp3 = (tmp1 << 24) | tmp2;
    this->m0 = navtools::TwosComp(tmp3, 32) * GPS_PI<T> * PowerOfTwo<-31, T>();

    // Word 6 and 7
    tmp3 = (subframe[5] & 0x3FFFC000) >> 14;
    this->cuc = navtools::TwosComp(tmp3, 16) * PowerOfTwo<-29, T>();  // bits 1-16
    tmp1 = (subframe[5] & 0x00003FC0) >> 6;                           // bits 17-24
    tmp2 = (subframe[6] & 0x3FFFFFC0) >> 6;                           // bits 1-24
    this->e = static_cast<T>((tmp1 << 24) | tmp2) * PowerOfTwo<-33, T>();

    // Word 8 and 9
    tmp3 = (subframe[7] & 0x3FFFC000) >> 14;
    this->cus = navtools::TwosComp(tmp3, 16) * PowerOfTwo<-29, T>();  // bits 1-16
    tmp1 = (subframe[7] & 0x00003FC0) >> 6;                           // bits 17-24
    tmp2 = (subframe[8] & 0x3FFFFFC0) >> 6;                           // bits 1-24
    this->sqrtA = static_cast<T>((tmp1 << 24) | tmp2) * PowerOfTwo<-19, T>();

    // Word 10
    this->toe = static_cast<T>((subframe[9] & 0x3FFFC000) >> 14) * PowerOfTwo<4, T>();  // bits 1-16
    // fit_interval_alert_flag = bool((subframe[9] & 0x00002000) >> 13);           // bit 17
  };

  /**
   * *=== LoadSubframe3 ===*
   * @brief Reads GPS LNAV subframe 3 (IS-GPS-200N pg. 82 & 105)
   */
  void LoadSubframe3() {
    uint32_t tmp1, tmp2, tmp3;

    // Word 1-2
    LoadPreamble();

    // word 3 and 4
    tmp3 = (subframe[2] & 0x3FFFC000) >> 14;
    this->cic = navtools::TwosComp(tmp3, 16) * PowerOfTwo<-29, T>();  // bits 1-16
    tmp1 = (subframe[2] & 0x00003FC0) >> 6;                           // bits 17-24
    tmp2 = (subframe[3] & 0x3FFFFFC0) >> 6;                           // bits 1-24
    tmp3 = (tmp1 << 24) | tmp2;
    this->omega0 = navtools::TwosComp(tmp3, 32) * GPS_PI<T> * PowerOfTwo<-31, T>();

    // Word 5 and 6
    tmp3 = (subframe[4] & 0x3FFFC000) >> 14;
    this->cis = navtools::TwosComp(tmp3, 16) * PowerOfTwo<-29, T>();  // bits 1-16
    tmp1 = (subframe[4] & 0x00003FC0) >> 6;                           // bits 17-24
    tmp2 = (subframe[5] & 0x3FFFFFC0) >> 6;                           // bits 1-24
    tmp3 = (tmp1 << 24) | tmp2;
    this->i0 = navtools::TwosComp(tmp3, 32) * GPS_PI<T> * PowerOfTwo<-31, T>();

    // word 7 and 8
    tmp3 = (subframe[6] & 0x3FFFC000) >> 14;
    this->crc = navtools::TwosComp(tmp3, 16) * PowerOfTwo<-5, T>();  // bits 1-16
    tmp1 = (subframe[6] & 0x00003FC0) >> 6;                          // bits 17-24
    tmp2 = (subframe[7] & 0x3FFFFFC0) >> 6;                          // bits 1-24
    tmp3 = (tmp1 << 24) | tmp2;
    this->omega = navtools::TwosComp(tmp3, 32) * GPS_PI<T> * PowerOfTwo<-31, T>();

    // Word 9
    tmp3 = (subframe[8] & 0x3FFFFFC0) >> 6;
    this->omegaDot = navtools::TwosComp(tmp3, 24) * GPS_PI<T> * PowerOfTwo<-43, T>();  // bits 1-24

    // Word 10
    this->iode = static_cast<T>((subframe[9] & 0x3FC00000) >> 22);  // bits 1-8
    tmp3 = (subframe[9] & 0x003FFF00) >> 8;
    this->iDot = navtools::TwosComp(tmp3, 14) * GPS_PI<T> * PowerOfTwo<-43, T>();  // bits 9-22
  };

  /**
   * *=== LoadSubframe4 ===*
   * @brief Reads GPS LNAV subframe 4 (IS-GPS-200N pg. )
   */
  void LoadSubframe4();

  /**
   * *=== LoadSubframe5 ===*
   * @brief Reads GPS LNAV subframe 5 (IS-GPS-200N pg. )
   */
  void LoadSubframe5();

  /**
   * *=== GetEphemerides ===*
   * @returns Current set of ephemerides
   */
  KeplerElements<T> GetEphemerides() {
    return KeplerElements<T>{this->iode,   this->iodc, this->toe,    this->toc,   this->tgd,
                             this->af2,    this->af1,  this->af0,    this->e,     this->sqrtA,
                             this->deltan, this->m0,   this->omega0, this->omega, this->omegaDot,
                             this->i0,     this->iDot, this->cuc,    this->cus,   this->cic,
                             this->cis,    this->crc,  this->crs,    this->ura,   this->health};
  };

  /**
   * *=== GetKlobuchar ===*
   * @returns Current set of ephemerides
   */
  KlobucharElements<T> GetKlobuchar() {
    return KlobucharElements<T>{
        this->a0, this->a1, this->a2, this->a3, this->b0, this->b1, this->b2, this->b3};
  };

  /**
   * *=== GetWeekNumber ===*
   * @returns Week number
   */
  uint16_t GetWeekNumber() {
    return week_;
  };

  /**
   * *=== GetTimeOfWeek ===*
   * @returns GPS time of week [gps seconds]
   */
  T GetTimeOfWeek() {
    return ToW_;
  };

  /**
   * *=== AreEphemeridesParsed ===*
   * @returns True|False based on if subframe 1,2 and 3 have been parsed
   */
  bool AreEphemeridesParsed() {
    return sub1_parsed_ & sub2_parsed_ & sub3_parsed_;
  };

 private:
  bool preamble_sync_{false};
  bool sub1_parsed_{false};
  bool sub2_parsed_{false};
  bool sub3_parsed_{false};
  uint32_t subframe[10];
  uint32_t prev_32_bits_{0};
  uint16_t word_cnt_{0};
  uint16_t bit_cnt_{0};
  uint16_t bits_since_preamble_{301};
  std::vector<uint16_t> preamble_idx_;
  uint16_t week_;
  T ToW_;
  // clang-format off
  inline static constexpr uint8_t GPS_D25[14] = {2,3,4,6,7,11,12,13,14,15,18,19,21,24};   // [1,2,3,5,6,10,11,12,13,14,17,18,20,23]
  inline static constexpr uint8_t GPS_D26[14] = {3,4,5,7,8,12,13,14,15,16,19,20,22,25};   // [2,3,4,6,7,11,12,13,14,15,18,19,21,24]
  inline static constexpr uint8_t GPS_D27[14] = {2,4,5,6,8,9,13,14,15,16,17,20,21,23};    // [1,3,4,5,7,8,12,13,14,15,16,19,20,22]
  inline static constexpr uint8_t GPS_D28[14] = {3,5,6,7,9,10,14,15,16,17,18,21,22,24};   // [2,4,5,6,8,9,13,14,15,16,17,20,21,23]
  inline static constexpr uint8_t GPS_D29[15] = {2,4,6,7,8,10,11,15,16,17,18,19,22,23,25};// [1,3,5,6,7,9,10,14,15,16,17,18,21,22,24]
  inline static constexpr uint8_t GPS_D30[13] = {4,6,7,9,10,11,12,14,16,20,23,24,25};     // [3,5,6,8,9,10,11,13,15,19,22,23,24]
  // clang-format on
};

}  // namespace satutils

#endif