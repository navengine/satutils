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

#include <cmath>
#include <navtools/constants.hpp>

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
DEFINE_FP_CONSTANT(SGP_CK2, 1.0 / 2.0 * navtools::J2<T> * SGP_AE<T> * SGP_AE<T>);
DEFINE_FP_CONSTANT(
    SGP_CK4, -3.0 / 8.0 * navtools::J4<T> * SGP_AE<T> * SGP_AE<T> * SGP_AE<T> * SGP_AE<T>);
DEFINE_FP_CONSTANT(SGP_A3OVK2, -navtools::J3<T> / SGP_CK2<T> * SGP_AE<T> * SGP_AE<T> * SGP_AE<T>);

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

}  // namespace satutils

#endif