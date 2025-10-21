/**
 * *ephemeris.hpp*
 *
 * =======  ========================================================================================
 * @file    satutils/ephemeris.hpp
 * @brief   Satellite ephemeris structures.
 * @date    January 2025
 * @author  Daniel Sturdivant <Auburn University GAVLAB>
 * @author  Blake Baker <Auburn University GAVLAB>
 * @ref     1. "IS-GPS-200N", 2022
 *          2. "A Software-Defined GPS and Galileo Receiver: A Single-Frequency Approach", 2007
 *              - Borre, Akos, Bertelsen, Rinder, Jensen
 * =======  ========================================================================================
 */

// TODO: add base ECEF based elements for Glonass constellation
// TODO: test tle navigation states output

#ifndef SATUTILS_EPHEMERIS_HPP
#define SATUTILS_EPHEMERIS_HPP

#include <Eigen/Dense>
#include <cmath>
#include <random>

#include <navtools/core/constants.hpp>
#include "satutils/time.hpp"

namespace satutils {

/**
 * @brief Ephemerides based on Keplerian orbital elements
 */
template <typename T = double>
struct KeplerElements
{
  T iode{std::nan("1")};      // Issue of data Ephemeris
  T iodc{std::nan("1")};      // Issue of data Clock
  T toe{std::nan("1")};       // Time of Ephemeris
  T toc{std::nan("1")};       // Time of Clock
  T tgd{std::nan("1")};       // Group delay
  T af2{std::nan("1")};       // 2nd order clock correction coef.
  T af1{std::nan("1")};       // 1st order clock correction coef.
  T af0{std::nan("1")};       // 0th order clock correction coef.
  T e{std::nan("1")};         // Eccentricity
  T sqrtA{std::nan("1")};     // Square root of semi-major axis
  T deltan{std::nan("1")};    // Mean motion difference
  T m0{std::nan("1")};        // Mean anomaly
  T omega0{std::nan("1")};    // Longitude of ascending node
  T omega{std::nan("1")};     // Argument of perigee
  T omegaDot{std::nan("1")};  // Rate of right ascension
  T i0{std::nan("1")};        // Inclination angle
  T iDot{std::nan("1")};      // Rate of inclination angle
  T cuc{std::nan("1")};       // Cos-harmonic correction coef. to the argument of latitude
  T cus{std::nan("1")};       // Sin-harmonic correction coef. to the argument of latitude
  T cic{std::nan("1")};       // Cos-harmonic correction coef. to the angle of inclination
  T cis{std::nan("1")};       // Sin-harmonic correction coef. to the angle of inclination
  T crc{std::nan("1")};       // Cos-harmonic correction coef. to the orbit radius
  T crs{std::nan("1")};       // Sin-harmonic correction coef. to the orbit radius
  T ura{std::nan("1")};       // Estimated accuracy
  T health{std::nan("1")};    // Satellite health
};


/**
 * @brief Ephemerides based on SGP4 orbital elements
 */
template <typename T = double>
struct Sgp4Elements {
  T catalog_id{std::nan("1")};  // Satellite catalog number
  T week{std::nan("1")};        // GPS week number of the TLE epoch
  T toe{std::nan("1")};         // Second of GPS week of the TLE epoch
  T Bstar{std::nan("1")};       // (BSTAR) drag/radian pressure coefficient
  T e0{std::nan("1")};          // (E0) eccentricity
  T omega0{std::nan("1")};      // (OMEGA0) right ascension of the ascending node
  T omega{std::nan("1")};       // (XNODE0) argument of perigee
  T i0{std::nan("1")};          // (XINCL) inclination angle
  T n0{std::nan("1")};          // (XN0) mean motion
  T nDot{std::nan("1")};        // (XNDT20) 1st derivative of mean motion
  T nDDot{std::nan("1")};       // (XNDD60) 2nd derivative of mean motion
  T m0{std::nan("1")};          // (XM0) mean anomaly
};


//! ------------------------------------------------------------------------------------------------

template <typename T = double>
class KeplerEphem : public KeplerElements<T> {
 public:
  KeplerEphem() {};
  KeplerEphem(const KeplerElements<T> &eph) {
    SetEphemerides(eph);
  };

  /**
   * *=== SetEphemerides ===*
   * @brief set the ephemeris elements
   */
  void SetEphemerides(const KeplerElements<T> &eph) {
    this->iode = eph.iode;
    this->iodc = eph.iodc;
    this->toe = eph.toe;
    this->toc = eph.toc;
    this->tgd = eph.tgd;
    this->af2 = eph.af2;
    this->af1 = eph.af1;
    this->af0 = eph.af0;
    this->e = eph.e;
    this->sqrtA = eph.sqrtA;
    this->deltan = eph.deltan;
    this->m0 = eph.m0;
    this->omega0 = eph.omega0;
    this->omega = eph.omega;
    this->omegaDot = eph.omegaDot;
    this->i0 = eph.i0;
    this->iDot = eph.iDot;
    this->cuc = eph.cuc;
    this->cus = eph.cus;
    this->cic = eph.cic;
    this->cis = eph.cis;
    this->crc = eph.crc;
    this->crs = eph.crs;
    this->ura = eph.ura;
    this->health = eph.health;
  };

  /**
   * *=== CalcNavStates ===
   * @brief Calculates satellite position, velocity, and acceleration using ephemeris
   * @param clk           Satellite clock corrections vector
   * @param pos           Satellite position vector
   * @param vel           Satellite velocity vector
   * @param acc           Satellite acceleration vector
   * @param transmit_time GPS system/transmitter time (TOW) of the satellite (accounting for transit
   *                      time from satellite to receiver) [gps seconds]
   * @return True|False based on success
   */
  template <bool calc_acc = false>
  void CalcNavStates(
      Eigen::Ref<Eigen::Vector<T, 3>> clk,
      Eigen::Ref<Eigen::Vector<T, 3>> pos,
      Eigen::Ref<Eigen::Vector<T, 3>> vel,
      Eigen::Ref<Eigen::Vector<T, 3>> acc,
      const T &transmit_time) {
    if (!initialized_) {
      init();
      initialized_ = true;
    }

    // satellite clock correction (sv time)
    T dt = CheckGpsSecond(transmit_time - this->toc);          // time from clock epoch
    T dt_sv = this->af0 + dt * (this->af1 + dt * this->af2);   // group delay depends on frequency
    T tk = CheckGpsSecond(transmit_time - dt_sv - this->toe);  // corrected time difference

    // mean anomaly
    T Mk = std::fmod(this->m0 + n_ * tk + nt::TWO_PI<T>, nt::TWO_PI<T>);

    // calculate eccentric anomaly
    T COSE, SINE, dE;
    T Ek = Mk;
    for (int i = 0; i < 10; i++) {
      COSE = std::cos(Ek);  // cosine of eccentric anomaly
      SINE = std::sin(Ek);  // sine of eccentric anomaly
      dE = (Mk - Ek + this->e * SINE) / (1.0 - this->e * COSE);
      if (std::abs(dE) < 1e-15) {
        break;
      }
      Ek += dE;
    }
    Ek = std::fmod(Ek + nt::TWO_PI<T>, nt::TWO_PI<T>);
    T DEN = 1.0 - this->e * COSE;  // common denominator

    // true anomaly
    // double vk = 2.0 * np.atan2(np.sqrt((1.0 + e) / (1.0 - e)) * np.tan(0.5 * Ek), 1.0);
    T vk = std::atan2(SQ1ME2_ * SINE, COSE - this->e);

    // argument of latitude
    T Phik = std::fmod(vk + this->omega, nt::TWO_PI<T>);
    T COS2PHI = std::cos(2.0 * Phik);
    T SIN2PHI = std::sin(2.0 * Phik);

    // corrections
    T uk = Phik + (this->cus * SIN2PHI + this->cuc * COS2PHI);      // argument of latitude
    T rk = A_ * DEN + (this->crs * SIN2PHI + this->crc * COS2PHI);  // radius
    T ik = this->i0 + this->iDot * tk + (this->cis * SIN2PHI + this->cic * COS2PHI);  // inclination
    T wk = std::fmod(
        this->omega0 +
            tk *
                (this->omegaDot - nt::WGS84_OMEGA<T>)-(nt::WGS84_OMEGA<T> * this->toe) +
            nt::TWO_PI<T>,
        nt::TWO_PI<T>);  // longitude of ascending node - (omega == w)
    T COSU = std::cos(uk);
    T SINU = std::sin(uk);
    T COSI = std::cos(ik);
    T SINI = std::sin(ik);
    T COSW = std::cos(wk);
    T SINW = std::sin(wk);

    // derivatives
    T EDotk = n_ / DEN;               // eccentric anomaly rate
    T vDotk = EDotk * SQ1ME2_ / DEN;  // true anomaly rate
    T iDotk = this->iDot +
              2.0 * vDotk * (this->cis * COS2PHI - this->cic * SIN2PHI);  // inclination angle rate
    T uDotk =
        vDotk *
        (1.0 + 2.0 * (this->cus * COS2PHI - this->cuc * SIN2PHI));  // argument of latitude rate
    T rDotk = (this->e * A_ * EDotk * SINE) +
              2.0 * vDotk * (this->crs * COS2PHI - this->crc * SIN2PHI);  // radius rate
    T wDotk = this->omegaDot - nt::WGS84_OMEGA<T>;  // longitude of ascending node rate

    // position calculations
    T xk_orb = rk * COSU;  // x-position in orbital frame
    T yk_orb = rk * SINU;  // y-position in orbital frame
    pos(0) = xk_orb * COSW - yk_orb * COSI * SINW;
    pos(1) = xk_orb * SINW + yk_orb * COSI * COSW;
    pos(2) = yk_orb * SINI;

    // velocity calculations
    T xDotk_orb = rDotk * COSU - rk * uDotk * SINU;  // x-velocity in orbital frame
    T yDotk_orb = rDotk * SINU + rk * uDotk * COSU;  // y-velocity in orbital frame
    vel(0) = -(xk_orb * wDotk * SINW) + (xDotk_orb * COSW) - (yDotk_orb * SINW * COSI) -
             (yk_orb * (wDotk * COSW * COSI - iDotk * SINW * SINI));
    vel(1) = (xk_orb * wDotk * COSW) + (xDotk_orb * SINW) + (yDotk_orb * COSW * COSI) -
             (yk_orb * (wDotk * SINW * COSI + iDotk * COSW * SINI));
    vel(2) = (yDotk_orb * SINI) + (yk_orb * iDotk * COSI);

    // relativistic clock calculations (user time)
    T FESQA = nt::WGS84_REL_F<T> * this->e * this->sqrtA;  // relativistic time factor
    clk(0) = dt_sv + (FESQA * SINE);
    // clk(0) = dt_sv - 2.0 * pos.dot(vel) / (nt::LIGHT_SPEED<T> * nt::LIGHT_SPEED<T>);
    clk(1) = this->af1 + (2.0 * this->af2 * dt) + (n_ * FESQA * COSE / DEN);

    if constexpr (calc_acc) {
      T F = -1.5 * nt::WGS84_J2<T> * (nt::WGS84_GM<T> / (rk * rk)) *
            std::pow(nt::WGS84_A<T> / rk, 2);
      T TMP1 = -nt::WGS84_GM<T> / (rk * rk * rk);
      T TMP2 = 5.0 * std::pow(pos(2) / rk, 2);
      T TMP3 = nt::WGS84_OMEGA<T> * nt::WGS84_OMEGA<T>;

      // state
      acc(0) = TMP1 * pos(0) + F * (1.0 - TMP2) * (pos(0) / rk) +
               2.0 * vel(1) * nt::WGS84_OMEGA<T> + pos(0) * TMP3;
      acc(1) = TMP1 * pos(1) + F * (1.0 - TMP2) * (pos(1) / rk) -
               2.0 * vel(0) * nt::WGS84_OMEGA<T> + pos(1) * TMP3;
      acc(2) = TMP1 * pos(2) + F * (3.0 - TMP2) * (pos(2) / rk);

      // clock
      clk(2) = 2.0 * this->af2 - (n_ * n_ * FESQA * SINE / (DEN * DEN));
    }
  };

  /**
   * *=== init ===*
   * @brief Initialize additional ephemeris constants
   */
  void init()
  {
    A_ = this->sqrtA * this->sqrtA;
    n0_ = std::sqrt(nt::WGS84_GM<T> / (A_ * A_ * A_));  // computed mean motion
    n_ = n0_ + this->deltan;                                  // corrected mean motion
    SQ1ME2_ = std::sqrt(1.0 - (this->e * this->e));                            // common eccentricity factor
  };

  /**
   * *=== GetEphemerides ===*
   * @returns Current set of ephemerides
   */
  KeplerElements<T> GetEphemerides() {
    return KeplerElements<T>{this->iode,   this->iodc, this->toe,    this->toc,   this->tgd,
                             this->af2,    this->af1,  this->af0,    this->e,     this->sqrtA,
                             this->deltan, this->m0,   this->omega0, this->omega, this->omegaDot,
                             this->i0,     this->iDot, this->cuc,    this->cus,   this->cic,
                             this->cis,    this->crc,  this->crs,    this->ura,   this->health};
  };

 protected:
  /**
   * @brief Extra constants used by navigation processor
   */
  bool initialized_{false};
  T A_;       // semi-major axis
  T n0_;      // computed mean motion
  T n_;       // corrected mean motion
  T SQ1ME2_;  // common eccentricity factor
};


/**
 * @brief Dedicated solely to calculating Kepler orbits with nothing extra
 */
template<typename T = double>
class KeplerOrbit
{
protected:
  T toe_{std::nan("1")};       // Reference Time
  T e_{std::nan("1")};         // Eccentricity
  T A_{std::nan("1")};         // Square root of semi-major axis
  T n_{std::nan("1")};         // Mean motion
  T m0_{std::nan("1")};        // Mean anomaly
  T omega0_{std::nan("1")};    // Longitude of ascending node
  T omega_{std::nan("1")};     // Argument of perigee
  T omega_dot_{std::nan("1")}; // Rate of right ascension
  T i0_{std::nan("1")};        // Inclination angle
  T i_dot_{std::nan("1")};     // Rate of inclination angle
  T cuc_{std::nan("1")};       // Cos-harmonic correction coef. to the argument of latitude
  T cus_{std::nan("1")};       // Sin-harmonic correction coef. to the argument of latitude
  T cic_{std::nan("1")};       // Cos-harmonic correction coef. to the angle of inclination
  T cis_{std::nan("1")};       // Sin-harmonic correction coef. to the angle of inclination
  T crc_{std::nan("1")};       // Cos-harmonic correction coef. to the orbit radius
  T crs_{std::nan("1")};       // Sin-harmonic correction coef. to the orbit radius
  T sq1me2_{std::nan("1")};

  static constexpr int E_ITERATIONS = 7;
public:
  typedef Eigen::Vector<T, 3> Vec;

  KeplerOrbit()
  {}

  void Randomize()
  {
    std::random_device rd; // random seed
    std::mt19937 gen(rd());
    std::uniform_real_distribution<T> dist(0.0,1.0);
    toe_ = dist(gen) * T(604784);
    e_ = dist(gen) * T(0.03);
    A_ = dist(gen) * T(60707964) + T(6400900);
    n_ = std::sqrt(nt::WGS84_GM<T> / (A_ * A_ * A_)) + (dist(gen) * std::pow(2.0,-27.0) - std::pow(2.0,-26.0));
    m0_ = (dist(gen) * T(2)) - T(1);
    omega0_ = (dist(gen) * T(2)) - T(1);
    omega_ = (dist(gen) * T(2)) - T(1);
    omega_dot_ = ((dist(gen) * T(2)) - T(1)) * std::pow(2.0,-20.0);
    i0_ = (dist(gen) * T(2)) - T(1);
    i_dot_ = ((dist(gen) * T(2)) - T(1)) * std::pow(2.0,-30.0);
    cuc_ = ((dist(gen) * T(2)) - T(1)) * std::pow(2.0,-14.0);
    cus_ = ((dist(gen) * T(2)) - T(1)) * std::pow(2.0,-14.0);
    cic_ = ((dist(gen) * T(2)) - T(1)) * std::pow(2.0,-14.0);
    cis_ = ((dist(gen) * T(2)) - T(1)) * std::pow(2.0,-14.0);
    crc_ = ((dist(gen) * T(2)) - T(1)) * std::pow(2.0,10.0);
    crs_ = ((dist(gen) * T(2)) - T(1)) * std::pow(2.0,10.0);
    sq1me2_ = std::sqrt(1.0 - (e_ * e_));
  }

  void Set(const KeplerElements<T>& elms)
  {
    toe_ = elms.toe;
    e_ = elms.e;
    A_ = elms.sqrtA * elms.sqrtA;
    n_ = std::sqrt(nt::WGS84_GM<T> / (A_ * A_ * A_)) + elms.deltan;
    m0_ = elms.m0;
    omega0_ = elms.omega0;
    omega_ = elms.omega;
    omega_dot_ = elms.omegaDot;
    i0_ = elms.i0;
    i_dot_ = elms.iDot;
    cuc_ = elms.cuc;
    cus_ = elms.cus;
    cic_ = elms.cic;
    cis_ = elms.cis;
    crc_ = elms.crc;
    crs_ = elms.crs;
    sq1me2_ = std::sqrt(1.0 - (e_ * e_));
  }

  KeplerOrbit(const KeplerElements<T>& elms)
  {
    this->Set(elms);
  }

  template <bool CalcPos = true, bool CalcVel = false, bool CalcAccel = false>
  void CalculatePVA(const T& transmit_time, Vec* pos, Vec* vel, Vec* acc) const
  {
    // satellite clock correction (sv time)
    T tk = CheckGpsSecond(transmit_time - toe_);

    // mean anomaly
    T Mk = std::fmod(m0_ + (n_ * tk) + nt::TWO_PI<T>, nt::TWO_PI<T>);

    // calculate eccentric anomaly
    T cos_E, sin_E, dE;
    T Ek = Mk;
    for (int i = 0; i < E_ITERATIONS; i++) {
      cos_E = std::cos(Ek);  // cosine of eccentric anomaly
      sin_E = std::sin(Ek);  // sine of eccentric anomaly
      dE = (Mk - Ek + (e_ * sin_E)) / (1.0 - (e_ * cos_E));
      if (std::abs(dE) < 1e-15) {
        break;
      }
      Ek += dE;
    }
    Ek = std::fmod(Ek + nt::TWO_PI<T>, nt::TWO_PI<T>);
    T r_scale = 1.0 - (e_ * cos_E);  // common denominator

    // true anomaly
    // double vk = 2.0 * np.atan2(np.sqrt((1.0 + e_) / (1.0 - e_)) * np.tan(0.5 * Ek), 1.0);
    T vk = std::atan2(sq1me2_ * sin_E, cos_E - e_);

    // argument of latitude
    T Phik = std::fmod(vk + omega_, nt::TWO_PI<T>);
    T COS2PHI = std::cos(2.0 * Phik);
    T SIN2PHI = std::sin(2.0 * Phik);

    // corrections
    T uk = Phik + (cus_ * SIN2PHI + cuc_ * COS2PHI);      // argument of latitude
    T rk = A_ * r_scale + (crs_ * SIN2PHI + crc_ * COS2PHI);  // radius
    T ik = i0_ + i_dot_ * tk + (cis_ * SIN2PHI + cic_ * COS2PHI);  // inclination
    T wk = std::fmod(
              omega0_ + tk * (omega_dot_ - nt::WGS84_OMEGA<T>)-(nt::WGS84_OMEGA<T> * toe_)
                + nt::TWO_PI<T>,
              nt::TWO_PI<T>
           );  // corrected longitude of ascending node
    T COSU = std::cos(uk);
    T SINU = std::sin(uk);
    T COSI = std::cos(ik);
    T SINI = std::sin(ik);
    T COSW = std::cos(wk);
    T SINW = std::sin(wk);

    // position calculations
    T xk_orb = rk * COSU;  // x-position in orbital frame
    T yk_orb = rk * SINU;  // y-position in orbital frame

    if constexpr (CalcPos || CalcAccel) {
      pos->operator()(0) = xk_orb * COSW - yk_orb * COSI * SINW;
      pos->operator()(1) = xk_orb * SINW + yk_orb * COSI * COSW;
      pos->operator()(2) = yk_orb * SINI;
    }

    // derivatives
    if constexpr (CalcVel || CalcAccel) {
      T EDotk = n_ / r_scale;               // eccentric anomaly rate
      T vDotk = EDotk * sq1me2_ / r_scale;  // true anomaly rate
      T iDotk = this->i_dot_ +
                2.0 * vDotk * (this->cis_ * COS2PHI - this->cic_ * SIN2PHI);  // inclination angle rate
      T uDotk =
          vDotk *
          (1.0 + 2.0 * (this->cus_ * COS2PHI - this->cuc_ * SIN2PHI));  // argument of latitude rate
      T rDotk = (this->e_ * A_ * EDotk * sin_E) +
                2.0 * vDotk * (this->crs_ * COS2PHI - this->crc_ * SIN2PHI);  // radius rate
      T wDotk = this->omega_dot_ - nt::WGS84_OMEGA<T>;  // longitude of ascending node rate

      // velocity calculations
      T xDotk_orb = rDotk * COSU - rk * uDotk * SINU;  // x-velocity in orbital frame
      T yDotk_orb = rDotk * SINU + rk * uDotk * COSU;  // y-velocity in orbital frame
      vel->operator()(0) = -(xk_orb * wDotk * SINW) + (xDotk_orb * COSW) - (yDotk_orb * SINW * COSI) -
               (yk_orb * (wDotk * COSW * COSI - iDotk * SINW * SINI));
      vel->operator()(1) = (xk_orb * wDotk * COSW) + (xDotk_orb * SINW) + (yDotk_orb * COSW * COSI) -
               (yk_orb * (wDotk * SINW * COSI + iDotk * COSW * SINI));
      vel->operator()(2) = (yDotk_orb * SINI) + (yk_orb * iDotk * COSI);

      if constexpr (CalcAccel) {
        T F = -1.5 * nt::WGS84_J2<T> * (nt::WGS84_GM<T> / (rk * rk)) *
              std::pow(nt::WGS84_A<T> / rk, 2);
        T TMP1 = -nt::WGS84_GM<T> / (rk * rk * rk);
        T TMP2 = 5.0 * std::pow(pos->operator()(2) / rk, 2);
        T TMP3 = nt::WGS84_OMEGA<T> * nt::WGS84_OMEGA<T>;

        // state
        acc->operator()(0) = TMP1 * pos->operator()(0) + F * (1.0 - TMP2) * (pos->operator()(0) / rk) +
                 2.0 * vel->operator()(1) * nt::WGS84_OMEGA<T> + pos->operator()(0) * TMP3;
        acc->operator()(1) = TMP1 * pos->operator()(1) + F * (1.0 - TMP2) * (pos->operator()(1) / rk) -
                 2.0 * vel->operator()(0) * nt::WGS84_OMEGA<T> + pos->operator()(1) * TMP3;
        acc->operator()(2) = TMP1 * pos->operator()(2) + F * (3.0 - TMP2) * (pos->operator()(2) / rk);
      }
    }
  }

  Vec P(const T& time)
  {
    Vec result;
    this->CalculatePVA<true,false,false>(time,&result,nullptr,nullptr);
    return result;
  }

  void PV(const T& time, Vec* pos, Vec* vel)
  {
    this->CalculatePVA<true,true,false>(time,pos,vel,nullptr);
  }
};


template <typename T = double>
class Sgp4Ephem : public Sgp4Elements<T> {
 public:
  /**
   * *=== SetEphemerides ===*
   * @brief set the ephemeris elements
   */
  void SetEphemerides(const Sgp4Elements<T> &eph) {
    this->catalog_id = eph.catalog_id;
    this->week = eph.week;
    this->toe = eph.toe;
    this->Bstar = eph.Bstar;
    this->e0 = eph.e0;
    this->omega0 = eph.omega0;
    this->omega = eph.omega;
    this->i0 = eph.i0;
    this->n0 = eph.n0;
    this->nDot = eph.nDot;
    this->nDDot = eph.nDDot;
    this->m0 = eph.m0;
  };

  /**
   * *=== GetEphemerides ===*
   * @returns Current set of ephemerides
   */
  Sgp4Elements<T> GetEphemerides() {
    return Sgp4Elements<T>{
        this->catalog_id,
        this->week,
        this->toe,
        this->Bstar,
        this->e0,
        this->omega0,
        this->omega,
        this->i0,
        this->n0,
        this->nDot,
        this->nDDot,
        this->m0};
  };

  /**
   * *=== CalcNavStates ===
   * @brief Calculates satellite position, velocity, and acceleration using ephemeris
   * @param pos           Satellite position vector
   * @param vel           Satellite velocity vector
   * @param transmit_time GPS system/transmitter time (TOW) of the satellite (accounting for
   * transit time from satellite to receiver) [gps seconds]
   * @return True|False based on success
   */
  void CalcNavStates(Eigen::Vector<T, 3> &pos, Eigen::Vector<T, 3> &vel, const T &transmit_time) {
    // https://apps.dtic.mil/sti/pdfs/ADA093554.pdf
    // https://celestrak.org/NORAD/documentation/spacetrk.pdf
    if (!initialized_) {
      init();
      initialized_ = true;
    }

    T TEMP, TEMP1, TEMP2, TEMP3, TEMP4, TEMP5, TEMP6;
    T dt = CheckGpsSecond(transmit_time - this->toe);

    // update for secular gravity and atmospheric drag
    T XMDF = this->m0 + XMDOT_ * dt;
    T OMGADF = this->omega0 + OMGDOT_ * dt;
    T XNODDF = this->omega + XNODOT_ * dt;
    T OMEGA = OMGADF;
    T XMP = XMDF;
    T TSQ = dt * dt;
    T XNODE = XNODDF + XNODCF_ * TSQ;
    T TEMPA = 1.0 - C1_ * dt;
    T TEMPE = this->Bstar * C4_ * dt;
    T TEMPL = T2COF_ * TSQ;
    if (!sgp4_simple_) {
      T C1SQ = C1_ * C1_;
      T D2 = 4.0 * AODP_ * TSI_ * C1SQ;
      TEMP = D2 * TSI_ * C1_ / 3.0;
      T D3 = (17.0 * AODP_ + S4_) * TEMP;
      T D4 = 0.5 * TEMP * AODP_ * TSI_ * (221.0 * AODP_ + 31.0 * S4_) * C1_;
      T T3COF = D2 + 2.0 * C1SQ;
      T T4COF = 0.25 * (3.0 * D3 + C1_ * (12.0 * D2 + 10.0 * C1SQ));
      T T5COF =
          0.2 * (3.0 * D4 + 12.0 * C1_ * D3 + 6.0 * D2 * D2 + 15.0 * C1SQ * (2.0 * D2 + C1SQ));
      T DELOMG = OMGCOF_ * dt;
      T DELM = XMCOF_ * (std::pow(1.0 + ETA_ * std::cos(XMDF), 3) - DELMO_);
      TEMP = DELOMG + DELM;
      XMP = XMDF + TEMP;
      OMEGA = OMGADF - TEMP;
      T TCUBE = TSQ * dt;
      T TFOUR = dt * TCUBE;
      TEMPA -= (D2 * TSQ - D3 * TCUBE - D4 * TFOUR);
      TEMPE -= (TEMPE + this->Bstar * C5_ * (std::sin(XMP) - SINMO_));
      TEMPL -= (TEMPL + T3COF * TCUBE + TFOUR * (T4COF + dt * T5COF));
    }
    T A = AODP_ * TEMPA * TEMPA;
    T E = this->e0 - TEMPE;
    T XL = XMP + OMEGA + XNODE + XNODP_ * TEMPL;
    T BETA = std::sqrt(1.0 - E * E);
    T XN = SGP_XKE<T> / std::pow(A, 1.5);

    // long period periodics
    T AXN = E * std::cos(OMEGA);
    TEMP = 1.0 / (A * BETA * BETA);
    T XLL = TEMP * XLCOF_ * AXN;
    T AYNL = TEMP * AYCOF_;
    T XLT = XL + XLL;
    T AYN = E * std::sin(OMEGA) + AYNL;

    // solve keplers equation
    T CAPU = std::fmod(XLT - XNODE, nt::TWO_PI<T>);
    TEMP2 = CAPU;
    T COSPW, SINPW, EPW;
    for (int i = 0; i < 10; i++) {
      SINPW = std::sin(TEMP2);
      COSPW = std::cos(TEMP2);
      TEMP3 = AXN * SINPW;
      TEMP4 = AYN * COSPW;
      TEMP5 = AXN * COSPW;
      TEMP6 = AYN * SINPW;
      EPW = (CAPU - TEMP4 + TEMP3 - TEMP2) / (1.0 - TEMP5 - TEMP6) + TEMP2;
      if (std::fabs(EPW - TEMP2) <= 1e-6) {
        break;
      }
      TEMP2 = EPW;
    }

    // short period preliminary quantities
    T ECOSE = TEMP5 + TEMP6;
    T ESINE = TEMP3 - TEMP4;
    T ELSQ = AXN * AXN + AYN * AYN;
    TEMP = 1.0 - ELSQ;
    T PL = A * TEMP;
    T R = A * (1.0 - ECOSE);
    TEMP1 = 1.0 / R;
    TEMP2 = A * TEMP1;
    T BETAL = std::sqrt(TEMP);
    TEMP3 = 1.0 / (1.0 + BETAL);
    T COSU = TEMP2 * (COSPW - AXN + AYN * ESINE * TEMP3);
    T SINU = TEMP2 * (SINPW - AYN - AXN * ESINE * TEMP3);
    T U = std::atan2(SINU, COSU);
    T SIN2U = 2.0 * SINU * COSU;
    T COS2U = 2.0 * COSU * COSU - 1.0;
    TEMP = 1.0 / PL;
    TEMP1 = SGP_CK2<T> * TEMP;
    TEMP2 = TEMP1 * TEMP;

    // update for short periodics
    T RK = R * (1.0 - 1.5 * TEMP2 * BETAL * X3THM1_) + 0.5 * TEMP1 * X1MTH2_ * COS2U;
    T UK = U - 0.25 * TEMP2 * X7THM1_ * SIN2U;
    T XNODEK = XNODE + 1.5 * TEMP2 * COSIO_ * SIN2U;
    T XINCK = this->i0 + 1.5 * TEMP2 * COSIO_ * SINIO_ * COS2U;

    // orientation vectors
    T SINUK = std::sin(UK);
    T COSUK = std::cos(UK);
    T SINIK = std::sin(XINCK);
    T COSIK = std::cos(XINCK);
    T SINNOK = std::sin(XNODEK);
    T COSNOK = std::cos(XNODEK);

    T XMX = -SINNOK * COSIK;
    T XMY = COSNOK * COSIK;
    T UX = XMX * SINUK + COSNOK * COSUK;
    T UY = XMY * SINUK + SINNOK * COSUK;
    T UZ = SINIK * SINUK;

    // position
    pos(0) = RK * UX;
    pos(1) = RK * UY;
    pos(2) = RK * UZ;
    pos *= (1000.0 * SGP_XKMPER<T> / SGP_AE<T>);  // [earth-radii] -> [m]

    // velocity
    // short period preliminary quantities
    T RDOT = SGP_XKE<T> * std::sqrt(A) * ESINE * TEMP1;
    T RFDOT = SGP_XKE<T> * std::sqrt(PL) * TEMP1;
    // update for short periodics
    T RDOTK = RDOT - XN * TEMP1 * X1MTH2_ * SIN2U;
    T RFDOTK = RFDOT + XN * TEMP1 * (X1MTH2_ * COS2U + 1.5 * X3THM1_);
    // orientation vectors
    T VX = XMX * COSUK - COSNOK * SINUK;
    T VY = XMY * COSUK - SINNOK * SINUK;
    T VZ = SINIK * COSUK;

    vel(0) = RDOTK * UX + RFDOTK * VX;
    vel(1) = RDOTK * UY + RFDOTK * VY;
    vel(2) = RDOTK * UZ + RFDOTK * VZ;
    vel *= (1000.0 * SGP_XKMPER<T> / SGP_AE<T> / 60.0);
  };

  /**
   * *=== init ===*
   * @brief initializes constants one time only
   */
  void init() {
    // TLE specific constants
    T E0SQ = this->e0 * this->e0;
    T BETA02 = 1.0 - E0SQ;
    T BETA0 = std::sqrt(BETA02);
    SINMO_ = std::sin(this->m0);
    COSIO_ = std::cos(this->i0);
    SINIO_ = std::sin(this->i0);
    T THETA2 = COSIO_ * COSIO_;
    T THETA4 = THETA2 * THETA2;
    X1MTH2_ = -THETA2 + 1.0;
    X3THM1_ = -1.0 + 3.0 * THETA2;
    sgp4_simple_ = false;

    // recover original mean motion and semi-major axis from elements (Constant per TLE)
    T a1 = std::pow((SGP_XKE<T> / this->n0), 2.0 / 3.0);
    T del1 = 1.5 * SGP_CK2<T> * X3THM1_ / (a1 * a1 * BETA0 * BETA02);
    T a0 = a1 * (1.0 - del1 * (1.0 / 3.0 + del1 * (1.0 + 134.0 / 81.0 * del1)));
    T del0 = 1.5 * SGP_CK2<T> * X3THM1_ / (a0 * a0 * BETA0 * BETA02);
    XNODP_ = this->n0 / (1.0 + del0);
    AODP_ = a0 / (1.0 - del0);

    // for perigee less than 220 km, equations are truncated to linear
    if ((AODP_ * (1.0 - this->e0) / SGP_AE<T>) < (220.0 / SGP_XKMPER<T> + SGP_AE<T>)) {
      sgp4_simple_ = true;
    }

    // initialize
    S4_ = SGP_S<T>;
    QOMS24_ = SGP_QOMS2T<T>;
    T perigee = (AODP_ * (1.0 - this->e0) - SGP_AE<T>)*SGP_XKMPER<T>;

    if (perigee < 156.0) {
      // For perigee between 98-156 km, the value of the constant s used in SGP4 is:
      S4_ = perigee - 78.0;

      if (perigee <= 98.0) {
        // For perigee below 98 km, the value of s is:
        S4_ = 20.0;
      }

      // If s is changed, (q0 - s*)^4 becomes:
      QOMS24_ = std::pow((120.0 - S4_) * SGP_AE<T> / SGP_XKMPER<T>, 4.0);
      S4_ /= (SGP_XKMPER<T> + SGP_AE<T>);
    }

    // constants given appropriate values of s_star and (q0 - s_star)^4
    T PINVSQ = 1.0 / (AODP_ * AODP_ * BETA02 * BETA02);
    TSI_ = 1.0 / (AODP_ - S4_);
    ETA_ = AODP_ * this->e0 * TSI_;
    T ETASQ = ETA_ * ETA_;
    T EETA = this->e0 * ETA_;
    T PSISQ = std::fabs(1.0 - ETASQ);
    T COEF = QOMS24_ * std::pow(TSI_, 4.0);
    T COEF1 = COEF / std::pow(PSISQ, 3.5);

    C2_ = COEF1 * XNODP_ *
          (AODP_ * (1.0 + 1.5 * ETASQ + EETA * (4.0 + ETASQ)) +
           0.75 * SGP_CK2<T> * TSI_ / PSISQ * X3THM1_ * (8.0 + 3.0 * ETASQ * (8.0 + ETASQ)));
    C1_ = this->Bstar * C2_;
    C3_ = COEF * TSI_ * SGP_A3OVK2<T> * XNODP_ * SGP_AE<T> * SINIO_ / this->e0;
    C4_ = 2.0 * XNODP_ * COEF1 * AODP_ * BETA02 *
          (ETA_ * (2.0 + 0.5 * ETASQ) + this->e0 * (0.5 + 2.0 * ETASQ) -
           2.0 * SGP_CK2<T> * TSI_ / (AODP_ * PSISQ) *
               (-3.0 * X3THM1_ * (1.0 - 2.0 * EETA + ETASQ * (1.5 + 0.5 * EETA)) +
                0.75 * X1MTH2_ * (2.0 * ETASQ - EETA * (1.0 + ETASQ)) *
                    std::cos(2.0 * this->omega0)));
    C5_ = 2.0 * COEF1 * AODP_ * BETA02 * (1.0 + 2.75 * (ETASQ + EETA) + EETA * ETASQ);

    T TEMP1 = 3.0 * SGP_CK2<T> * PINVSQ * XNODP_;
    T TEMP2 = TEMP1 * SGP_CK2<T> * PINVSQ;
    T TEMP3 = 1.25 * SGP_CK4<T> * PINVSQ * PINVSQ * XNODP_;
    XMDOT_ = XNODP_ + 0.5 * TEMP1 * BETA0 * X3THM1_ +
             0.0625 * TEMP2 * BETA0 * (13.0 - 78.0 * THETA2 + 137.0 * THETA4);
    T X1M5TH = -5.0 * THETA2 + 1.0;
    OMGDOT_ = 0.5 * TEMP1 * X1M5TH + 0.0625 * TEMP2 * (7.0 - 114.0 * THETA2 + 395.0 * THETA4);
    T XHDOT1 = TEMP1 * COSIO_;
    XNODOT_ = XHDOT1 +
              (0.5 * TEMP2 * (4.0 - 19.0 * THETA2) + 2.0 * TEMP3 * (3.0 - 7.0 * THETA2)) * COSIO_;
    OMGCOF_ = this->Bstar * C3_ * std::cos(this->omega0);
    XMCOF_ = -2.0 / 3.0 * COEF * this->Bstar * SGP_AE<T> / EETA;
    XNODCF_ = 3.5 * BETA02 * XHDOT1 * C1_;
    T2COF_ = 1.5 * C1_;
    XLCOF_ = 0.125 * SGP_A3OVK2<T> * SINIO_ * (3.0 + 5.0 * COSIO_) / (1.0 + COSIO_);
    AYCOF_ = 0.25 * SGP_A3OVK2<T> * SINIO_;
    DELMO_ = std::pow(1.0 + ETA_ * std::cos(this->m0), 3.0);
    X7THM1_ = 7.0 * THETA2 - 1.0;
  };

 protected:
  bool initialized_{false};
  bool sgp4_simple_{false};
  T tR_{0.0};
  T tR_dot_{0.0};
  T tR_ddot_{0.0};
  T COSIO_{0.0};
  T SINIO_{0.0};
  T SINMO_{0.0};
  T DELMO_{0.0};
  T X3THM1_{0.0};
  T X1MTH2_{0.0};
  T X7THM1_{0.0};
  T AODP_{0.0};
  T XNODP_{0.0};
  T S4_{0.0};
  T QOMS24_{0.0};
  T TSI_{0.0};
  T ETA_{0.0};
  T XMDOT_{0.0};
  T OMGDOT_{0.0};
  T XNODOT_{0.0};
  T XNODCF_{0.0};
  T T2COF_{0.0};
  T XMCOF_{0.0};
  T XLCOF_{0.0};
  T OMGCOF_{0.0};
  T AYCOF_{0.0};
  T C1_{0.0};
  T C2_{0.0};
  T C3_{0.0};
  T C4_{0.0};
  T C5_{0.0};
};

}  // namespace satutils

#endif
