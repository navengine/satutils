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

#include <initializer_list>
#include <iostream>
#include <cmath>
#include <array>

#include <navtools/math.hpp>
#include <satutils/common.hpp>
#include <satutils/constants.hpp>

namespace satutils {

//* ===== Broadcast Ionosphere Corrections ===================================================== *//

class Ionosphere {
public:
  virtual ~Ionosphere() = default;
  virtual void print() const = 0;
};

template <typename Float>
struct Klobuchar : public Ionosphere
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
  void print() const {
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

template<typename Float>
class TroposphericParameter
{
private:
  const Float (&values_)[5][2];
  static constexpr std::array<Float,5> lats_ = {
    deg2rad<Float>(15.0),
    deg2rad<Float>(30.0),
    deg2rad<Float>(45.0),
    deg2rad<Float>(60.0),
    deg2rad<Float>(75.0)
  };

  Float CalcParam(const Float& param_0, const Float& del_param, const Float& D, bool is_north) const
  {
    Float D_min = is_north ? static_cast<Float>(28) : static_cast<Float>(211);
    return param_0 - del_param * std::cos(TWO_PI<Float> * (D - D_min) / static_cast<Float>(365.25));
  }

  Float Interp(const std::size_t& i1, const std::size_t& i2, const Float& lat_mag) const
  {
    Float t = (lat_mag - lats_[i1]) / (lats_[i1+1] - lats_[i1]);
    return values_[i1][i2] + (values_[i1+1][i2] - values_[i1][i2]) * t;
  }

public:
  constexpr TroposphericParameter(const Float (&list)[5][2])
    : values_{list}
  {}

  Float operator()(const Float& latitude, const Float& day_of_year) const
  {
    bool is_north = (latitude >= static_cast<Float>(0));
    Float lat_mag = is_north ? latitude : -latitude;
    Float D_min = is_north ? static_cast<Float>(28) : static_cast<Float>(211);
    if (lat_mag <= lats_[0]) {
      return CalcParam(values_[0][0], values_[0][1], day_of_year, D_min);
    }
    else if (lat_mag < lats_[1]) {
      return CalcParam(Interp(0,0,lat_mag), Interp(0,1,lat_mag), day_of_year, is_north);
    }
    else if (lat_mag < lats_[2]) {
      return CalcParam(Interp(1,0,lat_mag), Interp(1,1,lat_mag), day_of_year, is_north);
    }
    else if (lat_mag < lats_[3]) {
      return CalcParam(Interp(2,0,lat_mag), Interp(2,1,lat_mag), day_of_year, is_north);
    }
    else if (lat_mag < lats_[4]) {
      return CalcParam(Interp(3,0,lat_mag), Interp(3,1,lat_mag), day_of_year, is_north);
    }
    else {
      return CalcParam(values_[4][0], values_[4][1], day_of_year, is_north);
    }
  }

  Float operator()(Float& p0, const Float& latitude, const Float& day_of_year) const
  {
    bool is_north = (latitude >= static_cast<Float>(0));
    Float lat_mag = is_north ? latitude : -latitude;
    if (lat_mag <= lats_[0]) {
      p0 = values_[0][0];
      return CalcParam(p0, values_[0][1], day_of_year, is_north);
    }
    else if (lat_mag < lats_[1]) {
      p0 = Interp(0,0,lat_mag);
      return CalcParam(p0, Interp(0,1,lat_mag), day_of_year, is_north);
    }
    else if (lat_mag < lats_[2]) {
      p0 = Interp(1,0,lat_mag);
      return CalcParam(p0, Interp(1,1,lat_mag), day_of_year, is_north);
    }
    else if (lat_mag < lats_[3]) {
      p0 = Interp(2,0,lat_mag);
      return CalcParam(p0, Interp(2,1,lat_mag), day_of_year, is_north);
    }
    else if (lat_mag < lats_[4]) {
      p0 = Interp(3,0,lat_mag);
      return CalcParam(p0, Interp(3,1,lat_mag), day_of_year, is_north);
    }
    else {
      p0 = values_[4][0];
      return CalcParam(p0, values_[4][1], day_of_year, is_north);
    }
  }

  Float Evaluate(const Float latitude, const Float day_of_year) const
  {
    return this->operator()(latitude, day_of_year);
  }
};

template<typename Float>
static constexpr Float tropospheric_pressures[5][2] = {
    {1013.25, 0.0},
    {1017.25, -3.75},
    {1015.75, -2.25},
    {1011.75, -1.75},
    {1013.00, -0.5}
  };

template<typename Float>
static constexpr Float tropospheric_temperatures[5][2] = {
    {299.65, 0.0},
    {294.15, 7.0},
    {283.15, 11.0},
    {272.15, 15.0},
    {263.65, 14.5}
  };

template<typename Float>
static constexpr Float tropospheric_vapor_pressures[5][2] = {
    {26.31, 0.0},
    {21.79, 8.85},
    {11.66, 7.24},
    {6.78, 5.36},
    {4.11, 3.39}
  };

template<typename Float>
static constexpr Float tropospheric_temp_lapse_rates[5][2] = {
    {6.3e-3, 0.0},
    {6.05e-3, 0.25e-3},
    {5.58e-3, 0.32e-3},
    {5.39e-3, 0.81e-3},
    {4.53e-3, 0.62e-3}
  };

template<typename Float>
static constexpr Float tropospheric_vapor_lapse_rates[5][2] = {
    {2.77, 0.0},
    {3.15, 0.33},
    {2.57, 0.46},
    {1.81, 0.74},
    {1.55, 0.3}
  };

template<typename Float>
static constexpr TroposphericParameter<Float> TropPressure(tropospheric_pressures<Float>);

template<typename Float>
static constexpr TroposphericParameter<Float> TropTemperature(tropospheric_temperatures<Float>);

template<typename Float>
static constexpr TroposphericParameter<Float> TropVaporPressure(tropospheric_vapor_pressures<Float>);

template<typename Float>
static constexpr TroposphericParameter<Float> TropTempLapseRate(tropospheric_temp_lapse_rates<Float>);

template<typename Float>
static constexpr TroposphericParameter<Float> TropVaporLapseRate(tropospheric_vapor_lapse_rates<Float>);

// day_of_year is the number of days into the year (from January 1)
// height is meters above mean sea level
// output is in meters
// valid for frequencies up to about 15 GHz
// Based on the UNB3 model introduced in (Collins, J., 1999. Assessment and Development of a Tropospheric Delay Model for Aircraft Users of the Global Positioning System)
// Source: https://gssc.esa.int/navipedia/index.php/Tropospheric_Delay
template<typename Float>
Float TroposphericDelay(const Float& lat_rad, const Float& height, const Float& el_rad, const Float& day_of_year)
{
  constexpr Float k1 = static_cast<Float>(77.604);
  constexpr Float k2 = static_cast<Float>(382000);
  constexpr Float Rd = static_cast<Float>(287.054);
  constexpr Float gm = static_cast<Float>(9.784);
  constexpr Float g = static_cast<Float>(9.80665);

  Float T0;
  Float T = TropTemperature<Float>(T0, lat_rad, day_of_year);
  Float beta = TropTempLapseRate<Float>(lat_rad, day_of_year);
  if (height > (T0 / beta))
    return static_cast<Float>(0);

  Float sin_e = std::sin(el_rad);
  Float M = 1.001 / std::sqrt(0.002001 + (sin_e * sin_e));
  Float P = TropPressure<Float>(lat_rad, day_of_year);
  Float e = TropVaporPressure<Float>(lat_rad, day_of_year);
  Float lambda_term = TropVaporLapseRate<Float>(lat_rad, day_of_year) + static_cast<Float>(1);

  Float height_factor = static_cast<Float>(1) - ((beta * height) / T);
  Float T_dry = ( (1.0e-6 * k1 * Rd * P) / gm )
              * std::pow( height_factor, g / (Rd * beta) );
  Float T_wet = ( (1.0e-6 * k2 * Rd) / ((lambda_term * gm) - (beta * Rd)) * (e / T) )
              * std::pow( height_factor, ((lambda_term * g) / (Rd * beta)) - 1.0 );
  return M * (T_dry + T_wet);
}

}  // namespace satutils
#endif
