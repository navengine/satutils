/**
 * *atmosphere.hpp*
 *
 * =======  ========================================================================================
 * @file    satutils/atmosphere.hpp
 * @brief   GNSS atmospheric corrections.
 * @date    January 2025
 * @author  Daniel Sturdivant <Auburn University GAVLAB>
 * @author  Blake Baker <Auburn University GAVLAB>
 * @ref     1. "IS-GPS-200N", 2022
 *          2. "Klobuchar Ionospheric Model" - Navipedia
 *          3. "Tropospheric Delay" - Navipedia
 *          4. "Mapping of Niell" - Navipedia
 *          5. "Principles of GNSS, Inertial, and Multisensor Integrated Navigation Systems", 2nd
 *              Edition, 2013 - Groves
 * =======  ========================================================================================
 */

#ifndef SATUTILS_ATMOSPHERE_HPP
#define SATUTILS_ATMOSPHERE_HPP

#include <cmath>
#include <navtools/constants.hpp>

#include "satutils/gnss-constants.hpp"

namespace satutils {

/**
 * @brief Struct containing polynomial coefficients for ionospheric corrections
 */
template <typename Tp = double>
struct KlobucharElements {
  Tp a0{std::nan("1")};
  Tp a1{std::nan("1")};
  Tp a2{std::nan("1")};
  Tp a3{std::nan("1")};
  Tp b0{std::nan("1")};
  Tp b1{std::nan("1")};
  Tp b2{std::nan("1")};
  Tp b3{std::nan("1")};
};

//! ------------------------------------------------------------------------------------------------

template <typename Tp = double>
class IonoModel : public KlobucharElements<Tp> {
 public:
  IonoModel() = default;
  IonoModel(const KlobucharElements<Tp> &klob) {
    SetKlobuchar(klob);
  };

  /**
   * *=== SetKlobuchar ===*
   * @brief Set the Klobuchar elements
   */
  void SetKlobuchar(const KlobucharElements<Tp> &klob) {
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
   * *=== CalcIonoDelay ===*
   * @brief Estimates the ionospheric delay based on the Klobuchar model
   * @param Iono  ionospheric time delay [s]
   * @param tow   GPS time of week in seconds [s]
   * @param lat   geodetic latitude [rad]
   * @param lon   geodetic longitude [rad]
   * @param az    azimuth angle to satellite [rad]
   * @param el    elevation angle to satellite [rad]
   * @param gamma frequency depending scaling factor
   */
  Tp CalcIonoDelay(
      const Tp &tow,
      const Tp &lat,
      const Tp &lon,
      const Tp &az,
      const Tp &el,
      const Tp &gamma = 1.0) {
    // elevation: radians to semi-circles
    Tp E = el / navtools::PI<Tp>;

    // 9. compute the slant factor
    Tp F = 1.0 + 16.0 * std::pow(0.53 - E, 3);

    Tp Iono;
    if (std::abs(F) <= 1.57) {
      // radians to semi-circles
      Tp phiu = lat / navtools::PI<Tp>;
      Tp lamu = lon / navtools::PI<Tp>;
      // Tp A = az / navtools::PI<Tp>;

      // 1. calculate earth-centered angle
      Tp psi = 0.0137 / (E + 0.11) - 0.022;

      // 2. calculate latitude of ionospheric pierce point
      Tp phiI = phiu + psi * std::cos(az);
      if (phiI > 0.416) {
        phiI = 0.416;
      } else if (phiI < -0.416) {
        phiI = -0.416;
      }

      // 3. compute longitude of ionospheric pierce point
      Tp lamI = lamu + psi * std::sin(az) / std::cos(phiI);

      // 4. find geomagnetic latitude of ionospheric pierce point
      Tp phim = phiI + 0.064 * std::cos(lamI - 1.617);
      Tp phim2 = phim * phim;
      Tp phim3 = phim2 * phim;

      // 5. find local time at the ionospheric pierce point
      Tp t = 43200.0 * lamI + std::fmod(tow, S_PER_DAY<Tp>);
      if (t < 0.0) {
        t += S_PER_DAY<Tp>;
      } else if (t > S_PER_DAY<Tp>) {
        t -= S_PER_DAY<Tp>;
      }

      // 6. compute the amplitude of ionospheric delay
      Tp AI = this->a0 + this->a1 * phim + this->a2 * phim2 + this->a3 * phim3;
      if (AI < 0.0) AI = 0.0;

      // 7. compute the period of ionospheric delay
      Tp PI = this->b0 + this->b1 * phim + this->b2 * phim2 + this->b3 * phim3;
      if (PI < 72000.0) PI = 72000.0;

      // 8. compute the phase of ionospheric delay
      Tp XI = navtools::TWO_PI<Tp> * (t - 50400.0) / PI;
      Tp XI2 = XI * XI;

      // 10. compute the ionospheric time delay
      Iono = (5e-9 + AI * (1.0 - XI2 / 2.0 + XI2 * XI2 / 24.0)) * F;

    } else {
      // 10. compute the ionospheric time delay
      Iono = 5e-9 * F;
    }

    Iono *= gamma;
    return Iono;
  };

  /**
   * *=== GetKlobuchar ===*
   * @returns Current set of ephemerides
   */
  KlobucharElements<Tp> GetKlobuchar() {
    return KlobucharElements<Tp>{
        this->a0, this->a1, this->a2, this->a3, this->b0, this->b1, this->b2, this->b3};
  };
};

template <typename Tp = double>
class TropoModel {
 public:
  /**
   * *=== CalcTropoDelay ===*
   * @brief Estimates the tropospheric delay based on dry and wet delays
   * @param Tropo tropospheric time delay [s]
   * @param DoY   current day of the year (Jan 1 = 0, Dec 31  = 365)
   * @param lat   geodetic latitude [rad]
   * @param h     geodetic altitude [m]
   * @param el    elevation angle to satellite [rad]
   */
  Tp CalcTropoDelay(const Tp &DoY, const Tp &lat, const Tp &h, const Tp &el) {
    // 1. Interpolate parameters
    Tp mag_lat_deg = navtools::RAD2DEG<Tp> * std::abs(lat);
    Eigen::Array<Tp, 1, 5> avg;
    Eigen::Array<Tp, 1, 5> delta;
    Eigen::Array<Tp, 1, 12> niell;
    Tp dx;
    if (mag_lat_deg < static_cast<Tp>(15)) {
      avg = ParamAvg.row(0);
      delta = SeasonalVar.row(0);
      niell = MappingOfNiell.row(0);

    } else if (mag_lat_deg < static_cast<Tp>(30)) {
      dx = mag_lat_deg - static_cast<Tp>(15);
      avg = interp<5>(ParamAvg.row(0), ParamAvg.row(1), dx);
      delta = interp<5>(SeasonalVar.row(0), SeasonalVar.row(1), dx);
      niell = interp<12>(MappingOfNiell.row(0), MappingOfNiell.row(1), dx);

    } else if (mag_lat_deg < static_cast<Tp>(45)) {
      dx = mag_lat_deg - static_cast<Tp>(30);
      avg = interp<5>(ParamAvg.row(1), ParamAvg.row(2), dx);
      delta = interp<5>(SeasonalVar.row(1), SeasonalVar.row(2), dx);
      niell = interp<12>(MappingOfNiell.row(1), MappingOfNiell.row(2), dx);

    } else if (mag_lat_deg < static_cast<Tp>(60)) {
      dx = mag_lat_deg - static_cast<Tp>(45);
      avg = interp<5>(ParamAvg.row(2), ParamAvg.row(3), dx);
      delta = interp<5>(SeasonalVar.row(2), SeasonalVar.row(3), dx);
      niell = interp<12>(MappingOfNiell.row(2), MappingOfNiell.row(3), dx);

    } else if (mag_lat_deg < static_cast<Tp>(75)) {
      dx = mag_lat_deg - static_cast<Tp>(60);
      avg = interp<5>(ParamAvg.row(3), ParamAvg.row(4), dx);
      delta = interp<5>(SeasonalVar.row(3), SeasonalVar.row(4), dx);
      niell = interp<12>(MappingOfNiell.row(3), MappingOfNiell.row(4), dx);

    } else {
      avg = ParamAvg.row(4);
      delta = SeasonalVar.row(4);
      niell = MappingOfNiell.row(4);
    }

    // 2. calculate parameter scale factor
    Tp Dmin = (lat >= static_cast<Tp>(0.0)) ? static_cast<Tp>(28.0) : static_cast<Tp>(211.0);
    Tp sf = std::cos(navtools::TWO_PI<Tp> * (DoY - Dmin) / static_cast<Tp>(365.25));

    // 3. calculate each parameter
    Tp P = avg(0) - delta(0) * sf;
    Tp T = avg(1) - delta(1) * sf;
    Tp e = avg(2) - delta(2) * sf;
    Tp B = avg(3) - delta(3) * sf;
    Tp l = avg(4) - delta(4) * sf + 1.0;

    // 4. calculate zero altitude vertical delay terms
    Tp T0dry = 1e-6 * k1 * Rd * P / gm;
    Tp T0wet = 1e-6 * k2 * Rd / (l * gm - B * Rd) * (e / T);

    // 5. calculate vertical delay terms
    Tp base = 1.0 - (B * h / T);
    Tp power = navtools::GRAVITY<Tp> / (Rd * B);
    Tp Tdry = std::pow(base, power) * T0dry;
    Tp Twet = std::pow(base, (l * power) - 1.0) * T0wet;

    // 6. calculate obliquity factor
    Tp sinE = std::sin(el);
    Tp M = 1.001 / std::sqrt(0.002001 + sinE * sinE);
    // Tp ad = niell(0) - niell(3) * sf;
    // Tp bd = niell(1) - niell(4) * sf;
    // Tp cd = niell(2) - niell(5) * sf;
    // Tp Mdry = niellmap(sinE, ad, bd, cd) +
    //           (1 / sinE - niellmap(sinE, niell(6), niell(7), niell(8))) * 1e-3 * h;
    // Tp Mwet = niellmap(sinE, niell(9), niell(10), niell(11));

    // 7. calculate tropospheric error
    return (Tdry + Twet) * M / navtools::LIGHT_SPEED<>;
    // Tropo = Tdry * Mdry + Twet * Mwet;
  };

 protected:
  /**
   * @brief constant parameters
   */
  inline static constexpr Tp k1 = 77.604;    // K/mbar
  inline static constexpr Tp k2 = 382000.0;  // K^2/mbar
  inline static constexpr Tp Rd = 287.054;   // J/Kg/K
  inline static constexpr Tp gm = 9.784;     // m/s^2
  inline static const Eigen::Array<Tp, 5, 5> ParamAvg{
      // clang-format off
      // P0       T0      e0     B0      l0
      {1013.25, 299.65, 26.31, 6.30e-3, 2.77}, 
      {1017.25, 294.15, 21.79, 6.05e-3, 3.15},
      {1015.75, 283.15, 11.66, 5.58e-3, 2.57}, 
      {1011.75, 272.15,  6.78, 5.39e-3, 1.81}, 
      {1013.00, 263.65,  4.11, 4.53e-3, 1.55},
      // clang-format on
  };
  inline static const Eigen::Array<Tp, 5, 5> SeasonalVar{
      // clang-format off
      // dP    dT    de      dB     dl
      { 0.00,  0.0, 0.00, 0.00   , 0.00},
      {-3.75,  7.0, 8.85, 0.25e-3, 0.33},
      {-2.25, 11.0, 7.24, 0.32e-3, 0.46},
      {-1.75, 15.0, 5.36, 0.81e-3, 0.74},
      {-0.50, 14.5, 3.39, 0.62e-3, 0.30},
      // clang-format on
  };
  inline static const Eigen::Array<Tp, 5, 12> MappingOfNiell{
      // clang-format off
      // 
      {1.2769934e-3, 2.9153695e-3, 62.610505e-3, 0.0         , 0.0         , 0.0         , 2.53e-5, 5.49e-3, 1.14e-3, 5.8021897e-4, 1.4275268e-3, 4.3472961e-2},
      {1.2683230e-3, 2.9152299e-3, 62.837393e-3, 1.2709626e-5, 2.1414949e-5, 9.0128400e-5, 2.53e-5, 5.49e-3, 1.14e-3, 5.6794847e-4, 1.5138625e-3, 4.6729510e-2},
      {1.2465397e-3, 2.9288445e-3, 63.721774e-3, 2.6523662e-5, 3.0160779e-5, 4.3497037e-5, 2.53e-5, 5.49e-3, 1.14e-3, 5.8118019e-4, 1.4572752e-3, 4.3908931e-2},
      {1.2196049e-3, 2.9022565e-3, 63.824265e-3, 3.4000452e-5, 7.2562722e-5, 84.795348e-5, 2.53e-5, 5.49e-3, 1.14e-3, 5.9727542e-4, 1.5007428e-3, 4.4626982e-2},
      {1.2045996e-3, 2.9024912e-3, 64.258455e-3, 4.1202191e-5, 11.723375e-5, 170.37206e-5, 2.53e-5, 5.49e-3, 1.14e-3, 6.1641692e-4, 1.7599082e-3, 5.4736038e-2},
      // clang-format on
  };

  /**
   * *=== interp ===*
   * @brief Linear interpolation
   */
  template <int S>
  inline constexpr Eigen::Array<Tp, 1, S> interp(
      const Eigen::Array<Tp, 1, S> &y0, const Eigen::Array<Tp, 1, S> &y1, const Tp &dx) {
    return y0 + dx * (y1 - y0) / static_cast<Tp>(15);
  }

  /**
   * *=== niellmap ===
   * @brief parameter map for mapping of niell
   */
  inline constexpr Tp niellmap(const Tp &sinE, const Tp &a, const Tp &b, const Tp &c) {
    return (1.0 + a / (1.0 + b / (1.0 + c))) / (sinE + a / (sinE + b / (sinE + c)));
  };
};

}  // namespace satutils

#endif
