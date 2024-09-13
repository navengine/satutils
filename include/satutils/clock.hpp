/**
|=========================================== clock.hpp ============================================|
|                                                                                                  |
|   @file     include/satutils/clock.hpp                                                           |
|   @brief    GNSS time constants and conversions.                                                 |
|   @date     July 2024                                                                            |
|                                                                                                  |
|==================================================================================================|
*/

#ifndef SATUTILS_CLOCK_HPP
#define SATUTILS_CLOCK_HPP

#include <chrono>
#include <cmath>
#include <iostream>

#include "satutils/constants.hpp"

namespace satutils {

//* ===== Time Points ========================================================================== *//

using time_point = std::chrono::system_clock::time_point;
using namespace std::chrono_literals;

constexpr time_point gps_reference_epoch =
        std::chrono::sys_days{std::chrono::January / 6 / 1980} + 0h + 0min + 0s;

struct GpsTime {
    int week;       // week number
    double second;  // time of week
};

//* ===== Conversions ========================================================================== *//

//! === DATE2GPSTIME ===
/// @brief      Convert a date-time into a GPS week and second of week
/// @param dt       std::chrono system clock time point
/// @param gpsT     GPS time
/// @returns    struct containing GPS week and second of week
void date2gpsTime(GpsTime& gpsT, const time_point& dt) {
    // count seconds from Jan. 6, 1980 00:00:00 (gps_time)
    //  - counts nanoseconds and converts to seconds.frac_seconds
    double gps_time_s =
            std::chrono::duration<double, std::nano>(dt - gps_reference_epoch).count() * 1e-9;

    // GPS week and time of week (second)
    gpsT.week = static_cast<int>(gps_time_s / 604800.0);
    gpsT.second = std::fmod(gps_time_s, 604800.0);
};
GpsTime date2gpsTime(const time_point& dt) {
    GpsTime gpsT;
    date2gpsTime(gpsT, dt);
    return gpsT;
};

//! === CHECKTIME ===
template <typename Float>
Float CheckTime(Float t) {
    if (t > HALF_WEEK<Float>) {
        t -= WEEK<Float>;
    } else if (t < -HALF_WEEK<Float>) {
        t += WEEK<Float>;
    }
    return t;
};

//* ===== Broadcast Timing Corrections ========================================================= *//

template <typename Float = double>
class Clock {
  public:
    virtual ~Clock() = default;

    //! === B ===
    /// @brief      calculates the bias provided by the clock model
    /// @param tow      time of week in seconds
    /// @param t        clock bias correction
    virtual void B(const Float& tow, Float& t) = 0;

    //! === BD ===
    /// @brief      calculates the bias and drift provided by the clock model
    /// @param tow      time of week in seconds
    /// @param t        clock bias correction
    /// @param t_dot    clock drift correction
    virtual void BD(const Float& tow, Float& t, Float& t_dot) = 0;

    //! === BDR ===
    /// @brief      calculates the bias, drift, and drift rate provided by the clock model
    /// @param tow      time of week in seconds
    /// @param t        clock bias correction
    /// @param t_dot    clock drift correction
    /// @param t_ddot   clock drift rate correction
    virtual void BDR(const Float& tow, Float& t, Float& t_dot, Float& t_ddot) = 0;

    //! === PRINT ===
    /// @brief      prints out each of the clock parameters
    virtual void print() const = 0;

  protected:
};

template <typename Float = double>
class PolyClock : public Clock<Float> {
  public:
    PolyClock() = default;
    // ~PolyClock() = default;

    Float a_f2{0.0};  // polynomial clock coefficients
    Float a_f1{0.0};  //
    Float a_f0{0.0};  //
    Float T_GD{0.0};  // Estimated group delay differential
    Float IODC{0.0};  // Issue of data, clock
    Float t_oc{0.0};  // clock data reference time

    //! === B ===
    void B(const Float& tow, Float& t) {
        polyfit<false, false>(tow);
        t = t_;
    }

    //! === BD ===
    void BD(const Float& tow, Float& t, Float& t_dot) {
        polyfit<true, false>(tow);
        t = t_;
        t_dot = t_dot_;
    }

    //! === BDR ===
    void BDR(const Float& tow, Float& t, Float& t_dot, Float& t_ddot) {
        polyfit<true, true>(tow);
        t = t_;
        t_dot = t_dot_;
        t_ddot = t_ddot_;
    }

    //! === PRINT ===
    void print() const {
        std::cout << "--- Polynomial Clock Fit ---" << '\n'
                  << "a_f2:  " << a_f2 << '\n'
                  << "a_f1:  " << a_f1 << '\n'
                  << "a_f0:  " << a_f0 << '\n'
                  << "T_GD:  " << T_GD << '\n'
                  << "t_oc:  " << t_oc << '\n'
                  << "IODC:  " << IODC << '\n'
                  << "----------------------------" << '\n';
    }

  protected:
    Float t_{0.0};
    Float t_dot_{0.0};
    Float t_ddot_{0.0};

    //! === POLYFIT ===
    template <bool CalcDrift, bool CalcRate>
    void polyfit(const Float& tow) {
        Float dt = CheckTime(tow - t_oc);
        t_ = a_f0 + (a_f1 * dt) + (a_f2 * dt * dt) - T_GD;
        if constexpr (CalcDrift) {
            t_dot_ = a_f1 + (2.0 * a_f2 * dt);

            if constexpr (CalcRate) {
                t_ddot_ = 2.0 * a_f2;
            }
        }
    }
};

}  // namespace satutils

#endif
