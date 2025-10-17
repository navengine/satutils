/**
 * *gnss-constants.hpp*
 *
 * =======  ========================================================================================
 * @file    satutils/gnss-constants.hpp
 * @brief   Satellite constellation constants.
 * @date    January 2025
 * @author  Daniel Sturdivant <Auburn University GAVLAB>
 * @author  Blake Baker <Auburn University GAVLAB>
 * @ref     1. "IS-GPS-200N", 2022
 *          2. "A Software-Defined GPS and Galileo Receiver: A Single-Frequency Approach", 2007
 *              - Borre, Akos, Bertelsen, Rinder, Jensen
 * =======  ========================================================================================
 */

#ifndef SATUTILS_GNSS_CONSTANTS_HPP
#define SATUTILS_GNSS_CONSTANTS_HPP

#include <iostream>

#include <cmath>
#include <navtools/core/macros.hpp>
#include <navtools/core/constants.hpp>

namespace satutils {

//* ===== GPS Constants ======================================================================== *//

DEFINE_FP_CONSTANT(GPS_PI, 3.1415926535898);                 //! PI defined by IS-GPS-200N
DEFINE_FP_CONSTANT(TWO_GPS_PI, 2.0 * GPS_PI<T>);             //! 2 * GPS_PI
DEFINE_FP_CONSTANT(GPS_L1_FREQUENCY, 1.57542e9);             //! GPS L1 frequency [Hz]
DEFINE_FP_CONSTANT(GPS_L2_FREQUENCY, 1.22760e9);             //! GPS L2 frequency [Hz]
DEFINE_FP_CONSTANT(GPS_L5_FREQUENCY, 1.17645e9);             //! GPS L5 frequency [Hz]
DEFINE_FP_CONSTANT(GPS_CA_CODE_RATE, 1.023e6);               //! GPS CA chipping rate [chips/s]
DEFINE_FP_CONSTANT(GPS_L2_CODE_RATE, 511.5e3);               //! GPS L2 chipping rate [chips/s]
DEFINE_FP_CONSTANT(GPS_L5_CODE_RATE, 10.23e6);               //! GPS L5 chipping rate [chips/s]
DEFINE_FP_CONSTANT(GPS_DATA_BIT_RATE, 50.0);                 //! NAV message bit rate [bits/s]
inline constexpr std::size_t GPS_CA_CODE_LENGTH = 1023;      //! GPS L1CA Code length [chips]
inline constexpr std::size_t GPS_L2CM_CODE_LENGTH = 10230;   //!
inline constexpr std::size_t GPS_L2CL_CODE_LENGTH = 767250;  //!
inline constexpr std::size_t GPS_L5_CODE_LENGTH = 10230;     //! GPS L5CA Code length [chips]
inline constexpr std::size_t LNAV_SUBFRAME_SIZE = 300;
inline constexpr std::size_t LNAV_WORD_SIZE = 30;
inline constexpr uint8_t LNAV_PREAMBLE_BITS = 0b10001011;
inline constexpr uint8_t LNAV_INV_PREAMBLE_BITS = 0b01110100;

//* ===== Galileo Constants ==================================================================== *//

DEFINE_FP_CONSTANT(GALILEO_E5_FREQUENCY, 1191.795e6);         //! Galileo E5 center frequency [Hz]
DEFINE_FP_CONSTANT(GALILEO_E5A_FREQUENCY, 1176.45e6);         //! Galileo E5a frequency [Hz]
DEFINE_FP_CONSTANT(GALILEO_E5B_FREQUENCY, 1207.14e6);         //! Galileo E5b frequency [Hz]
DEFINE_FP_CONSTANT(GALILEO_E6_FREQUENCY, 1278.75e6);          //! Galileo E6 frequency [Hz]
DEFINE_FP_CONSTANT(GALILEO_E6_CODE_RATE, 5.115e6);            //! Galileo E6 chipping rate [chips/s]
DEFINE_FP_CONSTANT(GALILEO_E1_DATA_RATE, 250);                //! E1 NAV message bit rate [bits/s]
DEFINE_FP_CONSTANT(GALILEO_E5A_DATA_RATE, 50);                //! E5a NAV message bit rate [bits/s]
DEFINE_FP_CONSTANT(GALILEO_E5B_DATA_RATE, 250);               //! E5b NAV message bit rate [bits/s]
DEFINE_FP_CONSTANT(GALILEO_E6_DATA_RATE, 1000);               //! E6 NAV message bit rate [bits/s]
inline constexpr std::size_t GALILEO_E1_CODE_LENGTH = 4092;   //! Galileo E1 Code length [chips]
inline constexpr std::size_t GALILEO_E5_CODE_LENGTH = 10230;  //! Galileo E5 Code length [chips]
inline constexpr std::size_t GALILEO_E6_CODE_LENGTH = 5115;   //! Galileo E6 Code length [chips]

//* ===== SGP Models =========================================================================== *//

DEFINE_FP_CONSTANT(SGP_AE, 1.0);                      // distance per earth radii
DEFINE_FP_CONSTANT(SGP_XKMPER, 6378.135);             // km per earth radii
DEFINE_FP_CONSTANT(SGP_S, 1.01222928);                // s
DEFINE_FP_CONSTANT(SGP_QOMS2T, 1.88027916e-9);        // (q0 - s)^4
DEFINE_FP_CONSTANT(SGP_XKE, 7.43669161331734132e-2);  // sqrt(G*M)
DEFINE_FP_CONSTANT(SGP_CK2, 1.0 / 2.0 * nt::WGS84_J2<T> * SGP_AE<T> * SGP_AE<T>);
DEFINE_FP_CONSTANT(
    SGP_CK4, -3.0 / 8.0 * nt::WGS84_J4<T> * SGP_AE<T> * SGP_AE<T> * SGP_AE<T> * SGP_AE<T>);
DEFINE_FP_CONSTANT(SGP_A3OVK2, -nt::WGS84_J3<T> / SGP_CK2<T> * SGP_AE<T> * SGP_AE<T> * SGP_AE<T>);

//* ===== Ephemeris Constants ================================================================== *//

// Used for generating upper and lower limits on different parameters
template <int Power, typename Float = double>
constexpr Float PowerOfTwo() {
  return static_cast<Float>(std::pow(2.0, Power));
}

DEFINE_FP_CONSTANT(HALF_WEEK, 302400.0);  //! half GPS week [s]
DEFINE_FP_CONSTANT(WEEK, 604800.0);       //! GPS week [s]
DEFINE_FP_CONSTANT(S_PER_DAY, 86400.0);   //! seconds in a day [s]
// TODO: add scale factors for ephemeris other than GPS


//* ===== Constellation Enums ================================================================== *//

enum ConstellationId { GPS, Galileo };

enum CodeId {
    GPSCA,
    GPSL1C,
    GPSL2CM,
    GPSL2CL,
    GPSL5I,
    GPSL5Q,
    GalileoE1OS,
    GalileoE5A,
    GalileoE5B,
    GalileoE6CS
};

template<CodeId Code, typename Float = double>
constexpr Float CodeRate()
{
  if constexpr ((Code == GPSCA) || (Code == GPSL1C)) {
    return GPS_CA_CODE_RATE<Float>;
  }
  else if constexpr ((Code == GPSL2CM) || (Code == GPSL2CL)) {
    return GPS_L2_CODE_RATE<Float>;
  }
  else if constexpr ((Code == GPSL5I) || (Code == GPSL5Q)) {
    return GPS_L5_CODE_RATE<Float>;
  }
  else if constexpr (Code == GalileoE1OS) {
    return GPS_CA_CODE_RATE<Float>;
  }
  else if constexpr ((Code == GalileoE5A) || (Code == GalileoE5B)) {
    return GPS_L5_CODE_RATE<Float>;
  }
  else if constexpr (Code == GalileoE6CS) {
    return GALILEO_E6_CODE_RATE<Float>;
  }
  else {
    std::cerr << "invalid CodeId used for CodeRate<CodeId,Float>()\n";
    return 0;
  }
}

template<typename Float = double>
constexpr Float CodeRate(CodeId code)
{
  switch(code) {
    case GPSCA:
    case GPSL1C:
    case GalileoE1OS:
      return GPS_CA_CODE_RATE<Float>;
    case GPSL2CM:
    case GPSL2CL:
      return GPS_L2_CODE_RATE<Float>;
    case GPSL5I:
    case GPSL5Q:
    case GalileoE5A:
    case GalileoE5B:
      return GPS_L5_CODE_RATE<Float>;
    case GalileoE6CS:
      return GALILEO_E6_CODE_RATE<Float>;
    default:
      std::cerr << "invalid CodeId used for CodeRate<Float>(CodeId)\n";
      return 0;
  }
}

template<CodeId Code>
constexpr std::size_t CodeLength()
{
  if constexpr ((Code == GPSCA) || (Code == GPSL1C)) {
    return GPS_CA_CODE_LENGTH;
  }
  else if constexpr (Code == GPSL2CM) {
    return GPS_L2CM_CODE_LENGTH;
  }
  else if constexpr (Code == GPSL2CL) {
    return GPS_L2CL_CODE_LENGTH;
  }
  else if constexpr ((Code == GPSL5I) || (Code == GPSL5Q)) {
    return GPS_L5_CODE_LENGTH;
  }
  else if constexpr (Code == GalileoE1OS) {
    return GALILEO_E1_CODE_LENGTH;
  }
  else if constexpr ((Code == GalileoE5A) || (Code == GalileoE5B)) {
    return GALILEO_E5_CODE_LENGTH;
  }
  else if constexpr (Code == GalileoE6CS) {
    return GALILEO_E6_CODE_LENGTH;
  }
  else {
    std::cerr << "invalid CodeId used for CodeLength<CodeId>()\n";
    return 0;
  }
}

constexpr std::size_t CodeLength(CodeId code)
{
  switch(code) {
    case GPSCA:
    case GPSL1C:
      return GPS_CA_CODE_LENGTH;
    case GPSL2CM:
      return GPS_L2CM_CODE_LENGTH;
    case GPSL2CL:
      return GPS_L2CL_CODE_LENGTH;
    case GPSL5I:
    case GPSL5Q:
    case GalileoE1OS:
    case GalileoE5A:
    case GalileoE5B:
    case GalileoE6CS:
    default:
      std::cerr << "invalid CodeId used for CodeLength(CodeId)\n";
      return 0;
  }
}

template<CodeId Code, typename Float = double>
constexpr Float CarrierFrequency()
{
  if constexpr ((Code == GPSCA) || (Code == GPSL1C)) {
    return GPS_L1_FREQUENCY<Float>;
  }
  else if constexpr ((Code == GPSL2CM) || (Code == GPSL2CL)) {
    return GPS_L2_FREQUENCY<Float>;
  }
  else if constexpr ((Code == GPSL5I) || (Code == GPSL5Q)) {
    return GPS_L5_FREQUENCY<Float>;
  }
  else if constexpr (Code == GalileoE1OS) {
    return GPS_L1_FREQUENCY<Float>;
  }
  else if constexpr ((Code == GalileoE5A) || (Code == GalileoE5B)) {
    return GALILEO_E5_FREQUENCY<Float>;
  }
  else if constexpr (Code == GalileoE6CS) {
    return GALILEO_E6_FREQUENCY<Float>;
  }
  else {
    std::cerr << "invalid CodeId used for CarrierFrequency<CodeId,Float>()\n";
    return 0;
  }
}

template <typename Float = double>
constexpr Float CarrierFrequency(CodeId code) {
  switch (code) {
    case GPSCA:
      return GPS_L1_FREQUENCY<Float>;
    case GPSL1C:
      return GPS_L1_FREQUENCY<Float>;
    case GPSL2CM:
      return GPS_L2_FREQUENCY<Float>;
    case GPSL2CL:
      return GPS_L2_FREQUENCY<Float>;
    case GPSL5I:
      return GPS_L5_FREQUENCY<Float>;
    case GPSL5Q:
      return GPS_L5_FREQUENCY<Float>;
    case GalileoE1OS:
      return GPS_L1_FREQUENCY<Float>;
    case GalileoE5A:
      return GALILEO_E5_FREQUENCY<Float>;
    case GalileoE5B:
      return GALILEO_E5_FREQUENCY<Float>;
    case GalileoE6CS:
      return GALILEO_E6_FREQUENCY<Float>;
    default:
      std::cerr << "invalid CodeId used for CarrierFrequency<Float>(CodeId)\n";
      return 0;
  }
}

template<CodeId Code, typename Float = double>
constexpr Float AngularFrequency()
{
  return CarrierFrequency<Code,Float>()
        * std::numbers::pi_v<Float> * static_cast<Float>(2);
}

template<typename Float = double>
constexpr Float AngularFrequency(CodeId code)
{
  return CarrierFrequency<Float>(code)
        * std::numbers::pi_v<Float> * static_cast<Float>(2);
}

template<CodeId Code, typename Float = double>
constexpr Float DataRate()
{
  if constexpr (
    (Code == GPSCA) || (Code == GPSL1C) || (Code == GPSL2CL) ||
    (Code == GPSL2CM) || (Code == GPSL5I) || (Code == GPSL5Q)
  ) {
    return GPS_DATA_BIT_RATE<Float>;
  }
  else if constexpr (Code == GalileoE1OS) {
    return GALILEO_E1_DATA_RATE<Float>;
  }
  else if constexpr (Code == GalileoE5A) {
    return GALILEO_E5A_DATA_RATE<Float>;
  }
  else if constexpr (Code == GalileoE5B) {
    return GALILEO_E5B_DATA_RATE<Float>;
  }
  else if constexpr (Code == GalileoE6CS) {
    return GALILEO_E6_DATA_RATE<Float>;
  }
  else {
    std::cerr << "invalid CodeId used for DataRate<CodeId,Float>()\n";
    return 0;
  }
}

template<typename Float = double>
constexpr Float DataRate(CodeId code)
{
  switch(code) {
    case GPSCA:
    case GPSL1C:
    case GPSL2CM:
    case GPSL2CL:
    case GPSL5I:
    case GPSL5Q:
      return GPS_DATA_BIT_RATE<Float>;
    case GalileoE1OS:
      return GALILEO_E1_DATA_RATE<Float>;
    case GalileoE5A:
      return GALILEO_E5A_DATA_RATE<Float>;
    case GalileoE5B:
      return GALILEO_E5B_DATA_RATE<Float>;
    case GalileoE6CS:
      return GALILEO_E6_DATA_RATE<Float>;
    default:
      std::cerr << "invalid CodeId used for DataRate<Float>(CodeId)\n";
      return 0;
  }
}

}  // namespace satutils

#endif
