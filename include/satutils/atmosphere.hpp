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
template <typename N = double>
struct KlobucharElements {
  N a0{std::nan("1")};
  N a1{std::nan("1")};
  N a2{std::nan("1")};
  N a3{std::nan("1")};
  N b0{std::nan("1")};
  N b1{std::nan("1")};
  N b2{std::nan("1")};
  N b3{std::nan("1")};
};

//! ------------------------------------------------------------------------------------------------

template <typename N = double>
class IonoModel : public KlobucharElements<N> {
 public:
  IonoModel<N>() = default;
  IonoModel<N>(const KlobucharElements<N> &klob) {
    SetKlobuchar(klob);
  };

  /**
   * *=== SetKlobuchar ===*
   * @brief Set the Klobuchar elements
   */
  void SetKlobuchar(const KlobucharElements<N> &klob) {
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
   * @param Iono  ionospheric time delay [m]
   * @param tow   GPS time of week in seconds [s]
   * @param lat   geodetic latitude [rad]
   * @param lon   geodetic longitude [rad]
   * @param az    azimuth angle to satellite [rad]
   * @param el    elevation angle to satellite [rad]
   */
  void CalcIonoDelay(N &Iono, const N &tow, const N &lat, const N &lon, const N &az, const N &el) {
    // elevation: radians to semi-circles
    N E = el / navtools::PI<N>;

    // 9. compute the slant factor
    N F = 1.0 + 16.0 * std::pow(0.53 - E, 3);

    if (std::abs(F) <= 1.57) {
      // radians to semi-circles
      N phiu = lat / navtools::PI<N>;
      N lamu = lon / navtools::PI<N>;
      // N A = az / navtools::PI<N>;

      // 1. calculate earth-centered angle
      N psi = 0.0137 / (E + 0.11) - 0.022;

      // 2. calculate latitude of ionospheric pierce point
      N phiI = phiu + psi * std::cos(az);
      if (phiI > 0.416) {
        phiI = 0.416;
      } else if (phiI < -0.416) {
        phiI = -0.416;
      }

      // 3. compute longitude of ionospheric pierce point
      N lamI = lamu + psi * std::sin(az) / std::cos(phiI);

      // 4. find geomagnetic latitude of ionospheric pierce point
      N phim = phiI + 0.064 * std::cos(lamI - 1.617);
      N phim2 = phim * phim;
      N phim3 = phim2 * phim;

      // 5. find local time at the ionospheric pierce point
      N t = 43200.0 * lamI + std::fmod(tow, S_PER_DAY<N>);
      if (t < 0.0) {
        t += S_PER_DAY<N>;
      } else if (t > S_PER_DAY<N>) {
        t -= S_PER_DAY<N>;
      }

      // 6. compute the amplitude of ionospheric delay
      N AI = this->a0 + this->a1 * phim + this->a2 * phim2 + this->a3 * phim3;
      if (AI < 0.0) AI = 0.0;

      // 7. compute the period of ionospheric delay
      N PI = this->b0 + this->b1 * phim + this->b2 * phim2 + this->b3 * phim3;
      if (PI < 72000.0) PI = 72000.0;

      // 8. compute the phase of ionospheric delay
      N XI = navtools::TWO_PI<N> * (t - 50400.0) / PI;
      N XI2 = XI * XI;

      // 10. compute the ionospheric time delay
      Iono = (5e-9 + AI * (1.0 - XI2 / 2.0 + XI2 * XI2 / 24.0)) * F;

    } else {
      // 10. compute the ionospheric time delay
      Iono = 5e-9 * F;
    }

    Iono *= navtools::LIGHT_SPEED<N>;
  };

  /**
   * *=== GetKlobuchar ===*
   * @returns Current set of ephemerides
   */
  KlobucharElements<N> GetKlobuchar() {
    return KlobucharElements<N>{
        this->a0, this->a1, this->a2, this->a3, this->b0, this->b1, this->b2, this->b3};
  };
};

template <typename N = double>
class TropoModel {
 public:
  /**
   * *=== CalcTropoDelay ===*
   * @brief Estimates the tropospheric delay based on the Klobuchar model
   * @param Tropo tropospheric time delay [m]
   * @param DoY   current day of the year (Jan 1 = 0, Dec 31  = 365)
   * @param lat   geodetic latitude [rad]
   * @param h     geodetic altitude [m]
   * @param el    elevation angle to satellite [rad]
   */
  void CalcTropoDelay(N &Tropo, const N &DoY, const N &lat, const N &h, const N &el) {
    // 1. Interpolate parameters
    N mag_lat_deg = navtools::RAD2DEG<N> * std::abs(lat);
    Eigen::Array<N, 1, 5> avg;
    Eigen::Array<N, 1, 5> delta;
    Eigen::Array<N, 1, 12> niell;
    N dx;
    if (mag_lat_deg < static_cast<N>(15)) {
      avg = ParamAvg.row(0);
      delta = SeasonalVar.row(0);
      niell = MappingOfNiell.row(0);

    } else if (mag_lat_deg < static_cast<N>(30)) {
      dx = mag_lat_deg - static_cast<N>(15);
      avg = interp<5>(ParamAvg.row(0), ParamAvg.row(1), dx);
      delta = interp<5>(SeasonalVar.row(0), SeasonalVar.row(1), dx);
      niell = interp<12>(MappingOfNiell.row(0), MappingOfNiell.row(1), dx);

    } else if (mag_lat_deg < static_cast<N>(45)) {
      dx = mag_lat_deg - static_cast<N>(30);
      avg = interp<5>(ParamAvg.row(1), ParamAvg.row(2), dx);
      delta = interp<5>(SeasonalVar.row(1), SeasonalVar.row(2), dx);
      niell = interp<12>(MappingOfNiell.row(1), MappingOfNiell.row(2), dx);

    } else if (mag_lat_deg < static_cast<N>(60)) {
      dx = mag_lat_deg - static_cast<N>(45);
      avg = interp<5>(ParamAvg.row(2), ParamAvg.row(3), dx);
      delta = interp<5>(SeasonalVar.row(2), SeasonalVar.row(3), dx);
      niell = interp<12>(MappingOfNiell.row(2), MappingOfNiell.row(3), dx);

    } else if (mag_lat_deg < static_cast<N>(75)) {
      dx = mag_lat_deg - static_cast<N>(60);
      avg = interp<5>(ParamAvg.row(3), ParamAvg.row(4), dx);
      delta = interp<5>(SeasonalVar.row(3), SeasonalVar.row(4), dx);
      niell = interp<12>(MappingOfNiell.row(3), MappingOfNiell.row(4), dx);

    } else {
      avg = ParamAvg.row(4);
      delta = SeasonalVar.row(4);
      niell = MappingOfNiell.row(4);
    }

    // 2. calculate parameter scale factor
    N Dmin = (lat >= static_cast<N>(0.0)) ? static_cast<N>(28.0) : static_cast<N>(211.0);
    N sf = std::cos(navtools::TWO_PI<N> * (DoY - Dmin) / static_cast<N>(365.25));

    // 3. calculate each parameter
    N P = avg(0) - delta(0) * sf;
    N T = avg(1) - delta(1) * sf;
    N e = avg(2) - delta(2) * sf;
    N B = avg(3) - delta(3) * sf;
    N l = avg(4) - delta(4) * sf + 1.0;

    // 4. calculate zero altitude vertical delay terms
    N T0dry = 1e-6 * k1 * Rd * P / gm;
    N T0wet = 1e-6 * k2 * Rd / (l * gm - B * Rd) * (e / T);

    // 5. calculate vertical delay terms
    N base = 1.0 - (B * h / T);
    N power = navtools::GRAVITY<N> / (Rd * B);
    N Tdry = std::pow(base, power) * T0dry;
    N Twet = std::pow(base, (l * power) - 1.0) * T0wet;

    // 6. calculate obliquity factor
    N sinE = std::sin(el);
    // N M = 1.001 / std::sqrt(0.002001 + sinE * sinE);
    N ad = niell(0) - niell(3) * sf;
    N bd = niell(1) - niell(4) * sf;
    N cd = niell(2) - niell(5) * sf;
    N Mdry = niellmap(sinE, ad, bd, cd) +
             (1 / sinE - niellmap(sinE, niell(6), niell(7), niell(8))) * 1e-3 * h;
    N Mwet = niellmap(sinE, niell(9), niell(10), niell(11));

    // 7. calculate tropospheric error
    // Tropo = (Tdry + Twet) * M;
    Tropo = Tdry * Mdry + Twet * Mwet;
  };

 protected:
  /**
   * @brief constant parameters
   */
  inline static constexpr N k1 = 77.604;    // K/mbar
  inline static constexpr N k2 = 382000.0;  // K^2/mbar
  inline static constexpr N Rd = 287.054;   // J/Kg/K
  inline static constexpr N gm = 9.784;     // m/s^2
  inline static const Eigen::Array<N, 5, 5> ParamAvg{
      // clang-format off
      // P0       T0      e0     B0      l0
      {1013.25, 299.65, 26.31, 6.30e-3, 2.77}, 
      {1017.25, 294.15, 21.79, 6.05e-3, 3.15},
      {1015.75, 283.15, 11.66, 5.58e-3, 2.57}, 
      {1011.75, 272.15,  6.78, 5.39e-3, 1.81}, 
      {1013.00, 263.65,  4.11, 4.53e-3, 1.55},
      // clang-format on
  };
  inline static const Eigen::Array<N, 5, 5> SeasonalVar{
      // clang-format off
      // dP    dT    de      dB     dl
      { 0.00,  0.0, 0.00, 0.00   , 0.00},
      {-3.75,  7.0, 8.85, 0.25e-3, 0.33},
      {-2.25, 11.0, 7.24, 0.32e-3, 0.46},
      {-1.75, 15.0, 5.36, 0.81e-3, 0.74},
      {-0.50, 14.5, 3.39, 0.62e-3, 0.30},
      // clang-format on
  };
  inline static const Eigen::Array<N, 5, 12> MappingOfNiell{
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
  inline constexpr Eigen::Array<N, 1, S> interp(
      const Eigen::Array<N, 1, S> &y0, const Eigen::Array<N, 1, S> &y1, const N &dx) {
    return y0 + dx * (y1 - y0) / static_cast<N>(15);
  }

  /**
   * *=== niellmap ===
   * @brief parameter map for mapping of niell
   */
  inline constexpr N niellmap(const N &sinE, const N &a, const N &b, const N &c) {
    return (1.0 + a / (1.0 + b / (1.0 + c))) / (sinE + a / (sinE + b / (sinE + c)));
  };
};

}  // namespace satutils

#endif
