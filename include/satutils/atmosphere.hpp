/**
|======================================== atmosphere.hpp ==========================================|
|                                                                                                  |
|   @file     include/satutils/atmosphere.hpp                                                      |
|   @brief    GNSS atmospheric corrections.                                                        |
|   @date     July 2024                                                                            |
|                                                                                                  |
|==================================================================================================|
*/

#ifndef SATUTILS_ATMOSPHERE_HPP
#define SATUTILS_ATMOSPHERE_HPP

#include <iostream>
#include <cmath>

#include <navtools/math.hpp>
#include <satutils/common.hpp>
#include <satutils/constants.hpp>

// TODO: add models here

namespace satutils {

//* ===== Broadcast Ionosphere Corrections ===================================================== *//

class Ionosphere {
public:
  virtual ~Ionosphere() = default;
  virtual void print() = 0;
};

template <typename Float>
struct Klobuchar : Ionosphere
{
public:
  Klobuchar() = default;
  // ~Klobuchar() = default;

  Float alpha_0{0.0};  // polynomial coefficients for ionospheric correction
  Float alpha_1{0.0};
  Float alpha_2{0.0};
  Float alpha_3{0.0};
  Float beta_0{0.0};
  Float beta_1{0.0};
  Float beta_2{0.0};
  Float beta_3{0.0};

  //! === PRINT ===
  void print() {
      std::cout << "--- Klobuchar Parameters ---" << '\n'
                << "alpha_0:  " << alpha_0 << '\n'
                << "alpha_1:  " << alpha_1 << '\n'
                << "alpha_2:  " << alpha_2 << '\n'
                << "alpha_3:  " << alpha_3 << '\n'
                << "beta_0:   " << beta_0 << '\n'
                << "beta_1:   " << beta_1 << '\n'
                << "beta_2:   " << beta_2 << '\n'
                << "beta_3:   " << beta_3 << '\n'
                << "----------------------------" << '\n';
  }

  // lat_semi - geodetic latitude (semicircles)
  // long_semi - geodetic longitude (semicircles)
  // az_rad - azimuth (radians)
  // el_semi - elevation angle (semicircles)
  Float GpsL1DelayImpl(const Float& gps_tow,
                       const Float& lat_semi, const Float& long_semi,
                       const Float& az_rad, const Float& el_semi) const
  {
    Float psi = ( static_cast<Float>(0.0137) / (el_semi + static_cast<Float>(0.11)) ) - static_cast<Float>(0.022); // semicircles
    Float phi_i = lat_semi + psi * std::cos(az_rad); // semicircles

    constexpr Float phi_thresh = static_cast<Float>(0.416);
    if (phi_i > phi_thresh)
      phi_i = phi_thresh;
    else if (phi_i < -phi_thresh)
      phi_i = -phi_thresh;

    Float lambda = long_semi
                 + ( (psi * std::sin(az_rad)) / std::cos(phi_i * PI<Float>) ); // semicircles
    Float phi_m = phi_i + static_cast<Float>(0.064)
                * std::cos( PI<Float> * (lambda - static_cast<Float>(1.617)) ); // semicircles

    Float t = (static_cast<Float>(43200) * lambda) + gps_tow;
    CircMod(t, static_cast<Float>(86400));

    Float amplitude = alpha_0 + phi_m * (alpha_1 + phi_m * (alpha_2 + phi_m * alpha_3));
    if (amplitude < static_cast<Float>(0)) amplitude = static_cast<Float>(0);

    Float period = beta_0 + phi_m * (beta_1 + phi_m * (beta_2 + phi_m * beta_3));
    if (period < static_cast<Float>(72000)) period = static_cast<Float>(72000);

    Float X_I = TWO_PI<Float> * (t - static_cast<Float>(50400)) / period;
    Float slant_factor = static_cast<Float>(1) + static_cast<Float>(16) * std::pow(0.53 - el_semi,3.0);

    if (X_I >= static_cast<Float>(1.57)) {
      return static_cast<Float>(5.0e-9) * slant_factor;
    }
    else {
      return slant_factor * (static_cast<Float>(5.0e-9) + (amplitude * (
             static_cast<Float>(1) - (static_cast<Float>(0.5) * std::pow(X_I,2.0)) + (std::pow(X_I,4.0) / static_cast<Float>(24.0))
        )));
    }
  }

  Float GpsL1Delay(const Float& gps_tow,
                   // const Eigen::Ref<const Vec3<Float>>& rx_ecef_pos,
                   // const Eigen::Ref<const Vec3<Float>>& tx_ecef_pos)
                   const Vec3<Float>& rx_ecef_pos,
                   const Vec3<Float>& tx_ecef_pos) const
  {
    Vec3<Float> rx_lla = ecef2lla(rx_ecef_pos);
    Vec3<Float> tx_aer = ecef2aer(rx_ecef_pos, tx_ecef_pos);
    return GpsL1DelayImpl(gps_tow, (rx_lla(0) / PI<Float>), (rx_lla(1) / PI<Float>), tx_aer(0), (tx_aer(1) / PI<Float>));
  }
  
  Float Delay(const Float& gps_tow,
                   // const Eigen::Ref<const Vec3<Float>>& rx_ecef_pos,
                   // const Eigen::Ref<const Vec3<Float>>& tx_ecef_pos,
                   const Vec3<Float>& rx_ecef_pos,
                   const Vec3<Float>& tx_ecef_pos,
                   const Float center_frequency) const
  {
    return std::pow(GPS_L1_FREQUENCY<Float> / center_frequency, 2.0)
           * GpsL1Delay(gps_tow, rx_ecef_pos, tx_ecef_pos);
  }
};

}  // namespace satutils

#endif
