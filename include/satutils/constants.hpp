/**
|======================================== constants.hpp ===========================================|
|                                                                                                  |
|   @file     include/satutils/constants.hpp                                                       |
|   @brief    Satellite constellation constants.                                                   |
|   @date     July 2024                                                                            |
|                                                                                                  |
|==================================================================================================|
*/

#ifndef SATUTILS_CONSTANTS_HPP
#define SATUTILS_CONSTANTS_HPP

#include <cmath>
#include <navtools/types.hpp>

namespace satutils {

//* ===== GPS Constants ======================================================================== *//

DEFINE_FP_CONSTANT(GPS_PI, 3.1415926535898);                 //! PI defined by IS-GPS-200N
template <typename Float>                                    //
inline constexpr Float GPS_TWOPI = 2.0 * GPS_PI<Float>;      //! 2 * GPS_PI
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

//* ===== Ephemeris Constants ================================================================== *//

DEFINE_FP_CONSTANT(TWO_P4, std::pow(2.0, 4));     //! 2^4
DEFINE_FP_CONSTANT(TWO_N5, std::pow(2.0, -5));    //! 2^-5
DEFINE_FP_CONSTANT(TWO_N19, std::pow(2.0, -19));  //! 2^-19
DEFINE_FP_CONSTANT(TWO_N29, std::pow(2.0, -29));  //! 2^-29
DEFINE_FP_CONSTANT(TWO_N31, std::pow(2.0, -31));  //! 2^-31
DEFINE_FP_CONSTANT(TWO_N33, std::pow(2.0, -33));  //! 2^-33
DEFINE_FP_CONSTANT(TWO_N43, std::pow(2.0, -43));  //! 2^-43
DEFINE_FP_CONSTANT(TWO_N55, std::pow(2.0, -55));  //! 2^-55
DEFINE_FP_CONSTANT(HALF_WEEK, 302400.0);          //! half GPS week [s]
DEFINE_FP_CONSTANT(WEEK, 604800.0);               //! GPS week [s]
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
    assert(false);
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
      assert(false);
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
    assert(false);
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
      assert(false);
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
    assert(false);
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
            assert(false);
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
    assert(false);
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
      assert(false);
      return 0;
  }
}

}  // namespace satutils

#endif
