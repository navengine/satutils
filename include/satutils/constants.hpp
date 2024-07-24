#ifndef NAVENGINE_INCLUDE_SATUTILS_CONSTANTS_HPP 
#define NAVENGINE_INCLUDE_SATUTILS_CONSTANTS_HPP  

#include <cassert>
#include <numbers>

#include <navtools/constants.hpp>

namespace satutils {

using navtools::TWO_PI;
using navtools::COMPLEX_I;
using navtools::LIGHT_SPEED;

// These are the exact precision as specified in the IS-GPS-200 documentation
DEFINE_FP_CONSTANT(GPS_PI, 3.1415926535898);
DEFINE_FP_CONSTANT(GPS_2PI, 2.0 * GPS_PI<T>);

DEFINE_FP_CONSTANT(GPS_L1_FREQUENCY, 1.57542e9);
DEFINE_FP_CONSTANT(GPS_L2_FREQUENCY, 1.22760e9);
DEFINE_FP_CONSTANT(GPS_L5_FREQUENCY, 1.17645e9);
DEFINE_FP_CONSTANT(GPS_CA_CODE_RATE, 1.023e6);
DEFINE_FP_CONSTANT(GPS_L2_CODE_RATE, 511.5e3);
DEFINE_FP_CONSTANT(GPS_L5_CODE_RATE, 10.23e6);
DEFINE_FP_CONSTANT(GPS_DATA_BIT_RATE, 50.0);

inline constexpr std::size_t GPS_CA_CODE_LENGTH = 1023;
inline constexpr std::size_t GPS_L2CM_CODE_LENGTH = 10230;
inline constexpr std::size_t GPS_L2CL_CODE_LENGTH = 767250;
inline constexpr std::size_t GPS_L5_CODE_LENGTH = 10230;

DEFINE_FP_CONSTANT(GALILEO_E5_FREQUENCY, 1191.795e6);
DEFINE_FP_CONSTANT(GALILEO_E5A_FREQUENCY, 1176.45e6);
DEFINE_FP_CONSTANT(GALILEO_E5B_FREQUENCY, 1207.14e6);
DEFINE_FP_CONSTANT(GALILEO_E6_FREQUENCY, 1278.75e6);
DEFINE_FP_CONSTANT(GALILEO_E6_CODE_RATE, 5.115e6);
DEFINE_FP_CONSTANT(GALILEO_E1_DATA_RATE, 250);
DEFINE_FP_CONSTANT(GALILEO_E5A_DATA_RATE, 50);
DEFINE_FP_CONSTANT(GALILEO_E5B_DATA_RATE, 250);
DEFINE_FP_CONSTANT(GALILEO_E6_DATA_RATE, 1000);

inline constexpr std::size_t GALILEO_E1_CODE_LENGTH = 4092;
inline constexpr std::size_t GALILEO_E5_CODE_LENGTH = 10230;
inline constexpr std::size_t GALILEO_E6_CODE_LENGTH = 5115;


enum ConstellationId {
  GPS,
  Galileo
};

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

template<CodeId Code, typename Scalar = double>
constexpr Scalar CodeRate()
{
  if constexpr ((Code == GPSCA) || (Code == GPSL1C)) {
    return GPS_CA_CODE_RATE<Scalar>;
  }
  else if constexpr ((Code == GPSL2CM) || (Code == GPSL2CL)) {
    return GPS_L2_CODE_RATE<Scalar>;
  }
  else if constexpr ((Code == GPSL5I) || (Code == GPSL5Q)) {
    return GPS_L5_CODE_RATE<Scalar>;
  }
  else if constexpr (Code == GalileoE1OS) {
    return GPS_CA_CODE_RATE<Scalar>;
  }
  else if constexpr ((Code == GalileoE5A) || (Code == GalileoE5B)) {
    return GPS_L5_CODE_RATE<Scalar>;
  }
  else if constexpr (Code == GalileoE6CS) {
    return GALILEO_E6_CODE_RATE<Scalar>;
  }
  else {
    assert(false);
    return 0;
  }
}

// Runtime version
template<typename Scalar = double>
Scalar CodeRate(CodeId code)
{
  switch(code) {
    case GPSCA:
    case GPSL1C:
    case GalileoE1OS:
      return GPS_CA_CODE_RATE<Scalar>;
    case GPSL2CM:
    case GPSL2CL:
      return GPS_L2_CODE_RATE<Scalar>;
    case GPSL5I:
    case GPSL5Q:
    case GalileoE5A:
    case GalileoE5B:
      return GPS_L5_CODE_RATE<Scalar>;
    case GalileoE6CS:
      return GALILEO_E6_CODE_RATE<Scalar>;
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

std::size_t CodeLength(CodeId code);

template<CodeId Code, typename Scalar = double>
constexpr Scalar CarrierFrequency()
{
  if constexpr ((Code == GPSCA) || (Code == GPSL1C)) {
    return GPS_L1_FREQUENCY<Scalar>;
  }
  else if constexpr ((Code == GPSL2CM) || (Code == GPSL2CL)) {
    return GPS_L2_FREQUENCY<Scalar>;
  }
  else if constexpr ((Code == GPSL5I) || (Code == GPSL5Q)) {
    return GPS_L5_FREQUENCY<Scalar>;
  }
  else if constexpr (Code == GalileoE1OS) {
    return GPS_L1_FREQUENCY<Scalar>;
  }
  else if constexpr ((Code == GalileoE5A) || (Code == GalileoE5B)) {
    return GALILEO_E5_FREQUENCY<Scalar>;
  }
  else if constexpr (Code == GalileoE6CS) {
    return GALILEO_E6_FREQUENCY<Scalar>;
  }
  else {
    assert(false);
    return 0;
  }
}

// Runtime version
template<typename Scalar = double>
Scalar CarrierFrequency(CodeId code)
{
  switch(code) {
    case GPSCA:
    case GPSL1C:
    case GalileoE1OS:
      return GPS_L1_FREQUENCY<Scalar>;
    case GPSL2CM:
    case GPSL2CL:
      return GPS_L2_FREQUENCY<Scalar>;
    case GPSL5I:
    case GPSL5Q:
      return GPS_L5_FREQUENCY<Scalar>;
    case GalileoE5A:
    case GalileoE5B:
      return GALILEO_E5_FREQUENCY<Scalar>;
    case GalileoE6CS:
      return GALILEO_E6_FREQUENCY<Scalar>;
    default:
      assert(false);
      return 0;
  }
}

template<CodeId Code, typename Scalar = double>
constexpr Scalar AngularFrequency()
{
  return CarrierFrequency<Code,Scalar>()
        * std::numbers::pi_v<Scalar> * static_cast<Scalar>(2);
}

template<typename Scalar = double>
Scalar AngularFrequency(CodeId code)
{
  return CarrierFrequency<Scalar>(code)
        * std::numbers::pi_v<Scalar> * static_cast<Scalar>(2);
}

template<CodeId Code, typename Scalar = double>
constexpr Scalar DataRate()
{
  if constexpr (
    (Code == GPSCA) || (Code == GPSL1C) || (Code == GPSL2CL) ||
    (Code == GPSL2CM) || (Code == GPSL5I) || (Code == GPSL5Q)
  ) {
    return GPS_DATA_BIT_RATE<Scalar>;
  }
  else if constexpr (Code == GalileoE1OS) {
    return GALILEO_E1_DATA_RATE<Scalar>;
  }
  else if constexpr (Code == GalileoE5A) {
    return GALILEO_E5A_DATA_RATE<Scalar>;
  }
  else if constexpr (Code == GalileoE5B) {
    return GALILEO_E5B_DATA_RATE<Scalar>;
  }
  else if constexpr (Code == GalileoE6CS) {
    return GALILEO_E6_DATA_RATE<Scalar>;
  }
  else {
    assert(false);
    return 0;
  }
}

template<typename Scalar = double>
Scalar DataRate(CodeId code)
{
  switch(code) {
    case GPSCA:
    case GPSL1C:
    case GPSL2CM:
    case GPSL2CL:
    case GPSL5I:
    case GPSL5Q:
      return GPS_DATA_BIT_RATE<Scalar>;
    case GalileoE1OS:
      return GALILEO_E1_DATA_RATE<Scalar>;
    case GalileoE5A:
      return GALILEO_E5A_DATA_RATE<Scalar>;
    case GalileoE5B:
      return GALILEO_E5B_DATA_RATE<Scalar>;
    case GalileoE6CS:
      return GALILEO_E6_DATA_RATE<Scalar>;
    default:
      assert(false);
      return 0;
  }
}


} // namespace satutils
#endif
