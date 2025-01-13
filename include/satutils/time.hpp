/**
 * *time.hpp*
 *
 * =======  ========================================================================================
 * @file    satutils/time.hpp
 * @brief   GNSS time conversions.
 * @date    January 2025
 * @author  Daniel Sturdivant <Auburn University GAVLAB>
 * @author  Blake Baker <Auburn University GAVLAB>
 * @ref     1. "IS-GPS-200N", 2022
 *          2. "A Software-Defined GPS and Galileo Receiver: A Single-Frequency Approach", 2007
 *              - Borre, Akos, Bertelsen, Rinder, Jensen
 * =======  ========================================================================================
 */

#ifndef SATUTILS_TIME_HPP
#define SATUTILS_TIME_HPP

#include <chrono>
#include <cmath>

#include "satutils/gnss-constants.hpp"

namespace satutils {

constexpr std::chrono::system_clock::time_point GPS_REF_EPOCH =
    std::chrono::sys_days{std::chrono::January / 6 / 1980};

/**
 * *=== Date2GpsTime ===*
 * @brief Convert a date-time into a GPS week and second of week
 * @param dt    std::chrono system clock time point
 * @param gpsT  GPS time
 * @return struct containing GPS week and second of week
 */
template <typename I, typename F>
void Date2GpsTime(I& week, F& second, const std::chrono::system_clock::time_point& dt) {
  // count seconds from Jan. 6, 1980 00:00:00 (gps_time)
  //  - counts nanoseconds and converts to seconds.frac_seconds
  F gps_time_s = std::chrono::duration<F, std::nano>(dt - GPS_REF_EPOCH).count() * 1e-9;

  // GPS week and time of week (second)
  week = static_cast<int>(gps_time_s / 604800.0);
  second = std::fmod(gps_time_s, 604800.0);
};

/**
 * *=== CheckGpsSecond ===*
 * @brief Ensure that the time is within a GPS week
 * @param t Time to check [gps seconds]
 */
template <typename T>
T CheckGpsSecond(T t) {
  if (t > HALF_WEEK<T>) {
    t -= WEEK<T>;
  } else if (t < -HALF_WEEK<T>) {
    t += WEEK<T>;
  }
  return t;
};

}  // namespace satutils

#endif