/**
|========================================= ephemeris.hpp ==========================================|
|                                                                                                  |
|   @file     include/satutils/ephemeris.hpp                                                       |
|   @brief    GNSS Ephemeris and TLE parameters.                                                   |
|   @date     July 2024                                                                            |
|                                                                                                  |
|==================================================================================================|
*/

//! TODO:
//!   - Add Glonass, Beidou, QZSS, and NavIC Ephemerides/Satellite objects
//!   - Add TLE ephemeris object
//!   - Add Iridium, Orbcomm, GlobalStar, and StarLink satellite objects

#ifndef SATUTILS_EPHEMERIS_HPP
#define SATUTILS_EPHEMERIS_HPP

#include <Eigen/Dense>
#include <chrono>
#include <cmath>
#include <string>

#include "navtools/constants.hpp"
#include "satutils/clock.hpp"

namespace satutils {

using namespace navtools;

//* ===== Ephemeris ============================================================================ *//

template <typename Float = double>
class Ephemeris {
  public:
    virtual ~Ephemeris() = default;

    //! === PVA ===
    /// @brief      calculates the position provided by the ephemeris model
    /// @param tow  time of week in seconds
    /// @param p    ECEF position [m]
    virtual void P(const Float& tow, Vec3<Float>& p) = 0;

    //! === PVA ===
    /// @brief      calculates the position and velocity provided by the ephemeris model
    /// @param tow  time of week in seconds
    /// @param p    ECEF position [m]
    /// @param v    ECEF velocity [m/s]
    virtual void PV(const Float& tow, Vec3<Float>& p, Vec3<Float>& v) = 0;

    //! === PVA ===
    /// @brief      calculates the position, velocity, and acceleration provided by the ephemeris
    ///             model
    /// @param tow  time of week in seconds
    /// @param p    ECEF position [m]
    /// @param v    ECEF velocity [m/s]
    /// @param a    ECEF acceleration [m/s^2]
    virtual void PVA(const Float& tow, Vec3<Float>& p, Vec3<Float>& v, Vec3<Float>& a) = 0;

    //! === B ===
    /// @brief      calculates the relativistic bias provided by the clock model
    /// @param tow      time of week in seconds
    /// @param t        relativistic clock bias correction
    virtual void B(const Float& tow, Float& t) = 0;

    //! === BD ===
    /// @brief      calculates the relativistic bias and drift provided by the clock model
    /// @param tow      time of week in seconds
    /// @param t        relativistic clock bias correction
    /// @param t_dot    relativistic clock drift correction
    virtual void BD(const Float& tow, Float& t, Float& t_dot) = 0;

    //! === BDR ===
    /// @brief      calculates the relativistic bias, drift, and drift rate provided by the clock
    ///             model
    /// @param tow      time of week in seconds
    /// @param t        relativistic clock bias correction
    /// @param t_dot    relativistic clock drift correction
    /// @param t_ddot   relativistic clock drift rate correction
    virtual void BDR(const Float& tow, Float& t, Float& t_dot, Float& t_ddot) = 0;

    //! === PRINT ===
    /// @brief      prints out each of the ephemeris parameters
    virtual void print() = 0;

  protected:
};

//* ===== Keplerian Ephemeris ================================================================== *//

// // GPS Ephemeris
// int leap_seconds{0};
// int week{0};
// bool l2c_flag{false};
// bool l2p_flag{false};
// Float accuracy{0.0};
// uint8_t health{0u};

// // Galileo Ephemeris
// int leap_seconds{0};
// int week{0};
// Float accuracy{0.0};           // SISA
// uint8_t health{0u};             //
// uint16_t data_source_flag{0u};  //
// Float BGD_e5a_e1{0.0};         //
// Float BGD_e5b_e1{0.0};         //

template <typename Float = double>
class KeplerEphemeris : public Ephemeris<Float> {
  public:
    KeplerEphemeris() = default;
    // ~KeplerEphemeris() = default;

    int week{0};           //
    Float M_0{0.0};        // mean anomaly at reference time
    Float delta_n{0.0};    // mean motion difference from computed value
    Float e{0.0};          // eccentricity
    Float A{0.0};          // semi major axis
    Float OMEGA_0{0.0};    // longitude of ascending node of orbital plane at weekly epoch
    Float i_0{0.0};        // inclination angle at reference time
    Float omega{0.0};      // argument of perigee
    Float OMEGA_DOT{0.0};  // rate of right ascension
    Float IDOT{0.0};       // rate of inclination angle
    Float C_uc{0.0};  // amplitude of the cos-harmonic correction term to the argument of latitude
    Float C_us{0.0};  // amplitude of the sin-harmonic correction term to the argument of latitude
    Float C_rc{0.0};  // amplitude of the cos-harmonic term to the orbit radius
    Float C_rs{0.0};  // amplitude of the sin-harmonic term to the orbit radius
    Float C_ic{0.0};  // amplitude of the cos-harmonic term to the angle of inclination
    Float C_is{0.0};  // amplitude of the sin-harmonic term to the angle of inclination
    Float t_oe{0.0};  // reference time, ephemeris
    Float IODE{0.0};  // issue of data, ephemeris

    //! === P ===
    void P(const Float& tow, Vec3<Float>& p) {
        kepler<false, false>(tow);
        p = r_eb_e_;
    }

    //! === PV ===
    void PV(const Float& tow, Vec3<Float>& p, Vec3<Float>& v) {
        kepler<true, false>(tow);
        p = r_eb_e_;
        v = v_eb_e_;
    }

    //! === PVA ===
    void PVA(const Float& tow, Vec3<Float>& p, Vec3<Float>& v, Vec3<Float>& a) {
        kepler<true, true>(tow);
        p = r_eb_e_;
        v = v_eb_e_;
        a = a_eb_e_;
    }

    //! === B ===
    void B(const Float& tow, Float& tR) {
        RelativeClock<false, false>(tow);
        tR = tR_;
    }

    //! === BD ===
    void BD(const Float& tow, Float& tR, Float& tR_dot) {
        RelativeClock<true, false>(tow);
        tR = tR_;
        tR_dot = tR_dot_;
    }

    //! === BDR ===
    void BDR(const Float& tow, Float& tR, Float& tR_dot, Float& tR_ddot) {
        RelativeClock<true, true>(tow);
        tR = tR_;
        tR_dot = tR_dot_;
        tR_ddot = tR_ddot_;
    }

    //! === PRINT ==
    void print() {
        std::cout << "--- Kepler Ephemerides ---" << std::endl
                  << "M_0:       " << M_0 << std::endl
                  << "delta_n:   " << delta_n << std::endl
                  << "e:         " << e << std::endl
                  << "A:         " << A << std::endl
                  << "OMEGA_0:   " << OMEGA_0 << std::endl
                  << "i_0:       " << i_0 << std::endl
                  << "omega:     " << omega << std::endl
                  << "OMEGA_DOT: " << OMEGA_DOT << std::endl
                  << "IDOT:      " << IDOT << std::endl
                  << "C_uc:      " << C_uc << std::endl
                  << "C_us:      " << C_us << std::endl
                  << "C_rc:      " << C_rc << std::endl
                  << "C_rs:      " << C_rs << std::endl
                  << "C_ic:      " << C_ic << std::endl
                  << "C_is:      " << C_is << std::endl
                  << "t_oe:      " << t_oe << std::endl
                  << "IODE:      " << IODE << std::endl
                  << "--------------------------" << std::endl;
    }

  protected:
    Vec3<Float> r_eb_e_{Vec3<Float>::Zero()};  // ecef position
    Vec3<Float> v_eb_e_{Vec3<Float>::Zero()};  // ecef velocity
    Vec3<Float> a_eb_e_{Vec3<Float>::Zero()};  // ecef acceleration
    Float tR_{0.0};
    Float tR_dot_{0.0};
    Float tR_ddot_{0.0};
    bool initialized_{false};
    Float n_{0.0};
    Float SQ1ME2_{0.0};
    Float FESQA_{0.0};
    Float cosE_{0.0};
    Float sinE_{0.0};
    Float den_{0.0};

    //! === INIT ===
    void init() {
        n_ = std::sqrt(GM<Float> / (A * A * A)) + delta_n;
        SQ1ME2_ = std::sqrt(1.0 - WGS84_E2<Float>);
        FESQA_ = F<Float> * e * std::sqrt(A);
    }

    //! === ECCENTRICANOMALY ===
    Float EccentricAnomalyM(const Float& M_k) {
        Float dE = 1.0;
        Float E_k = M_k;
        for (size_t i = 0; i < 10; i++) {
            dE = (M_k - E_k + (e * std::sin(E_k))) / (1.0 - (e * std::cos(E_k)));
            E_k += dE;
            if (dE < 1e-15) break;
        }
        return std::fmod(E_k + GPS_TWOPI<Float>, GPS_TWOPI<Float>);
    }
    Float EccentricAnomalyT(const Float& t_k) {
        Float M_k = std::fmod(M_0 + n_ * t_k + GPS_TWOPI<Float>, GPS_TWOPI<Float>);
        return EccentricAnomalyM(M_k);
    }

    //! === RELATIVECLOCK ===
    template <bool CalcDrift, bool CalcRate>
    void RelativeClock(const Float& tow) {
        Float t_k = CheckTime(tow - t_oe);
        Float E_k = EccentricAnomalyT(t_k);
        cosE_ = std::cos(E_k);
        sinE_ = std::sin(E_k);
        den_ = 1.0 - e * cosE_;

        // relativistic clock corrections
        tR_ = FESQA_ * sinE_;
        if constexpr (CalcDrift) {
            tR_dot_ = n_ * FESQA_ * cosE_ / den_;
            if constexpr (CalcRate) {
                tR_ddot_ = -n_ * n_ * FESQA_ * sinE_ / (den_ * den_);
            }
        }
    }

    //! === KEPLER ===
    template <bool CalcVel, bool CalcAccel>
    void kepler(const Float& tow) {
        if (!initialized_) {
            init();
            initialized_ = true;
        }

        // time since ephemeris
        Float t_k = CheckTime(tow - t_oe);

        // mean anomaly
        Float M_k = std::fmod(M_0 + n_ * t_k + GPS_TWOPI<Float>, GPS_TWOPI<Float>);

        // eccentric anomaly
        Float E_k = EccentricAnomalyM(M_k);
        cosE_ = std::cos(E_k);
        sinE_ = std::sin(E_k);
        den_ = 1.0 - e * cosE_;

        // true anomaly
        Float v_k = 2.0 * std::atan2(std::sqrt((1.0 + e) / (1.0 - e)) * std::tan(0.5 * E_k), 1.0);

        // argument of latitude
        Float Phi_k = std::fmod(v_k + omega, GPS_TWOPI<Float>);
        Float cos2Phi = std::cos(2.0 * Phi_k);
        Float sin2Phi = std::sin(2.0 * Phi_k);

        // corrected argument latitude, radius, inclination, and latitude of ascending node
        Float u_k = Phi_k + (C_us * sin2Phi + C_uc * cos2Phi);
        Float r_k = A * den_ + (C_rs * sin2Phi + C_rc * cos2Phi);
        Float i_k = i_0 + IDOT * t_k + (C_is * sin2Phi + C_ic * sin2Phi);
        Float OMEGA_k = std::fmod(
                OMEGA_0 + (OMEGA_DOT - WGS84_OMEGA<Float>)*t_k - (WGS84_OMEGA<Float> * t_oe) +
                        GPS_TWOPI<Float>,
                GPS_TWOPI<Float>);
        Float cosu = std::cos(u_k);
        Float sinu = std::sin(u_k);
        Float cosi = std::cos(i_k);
        Float sini = std::sin(i_k);
        Float cosOMEGA = std::cos(OMEGA_k);
        Float sinOMEGA = std::sin(OMEGA_k);

        // position
        Float x_k_prime = r_k * cosu;  // x-position in orbital frame
        Float y_k_prime = r_k * sinu;  // y-position in orbital frame
        r_eb_e_(1) = x_k_prime * cosOMEGA - y_k_prime * cosi * sinOMEGA;
        r_eb_e_(2) = x_k_prime * sinOMEGA + y_k_prime * cosi * cosOMEGA;
        r_eb_e_(3) = y_k_prime * sini;

        if constexpr (CalcVel) {
            // derivatives
            Float Edot_k = n_ / den_;
            Float vdot_k = Edot_k * SQ1ME2_ / den_;
            Float idot_k = IDOT + 2.0 * vdot_k * (C_is * cos2Phi - C_ic * sin2Phi);
            Float udot_k = vdot_k * (1.0 + 2.0 * (C_us * cos2Phi - C_uc * sin2Phi));
            Float rdot_k =
                    (e * A * Edot_k * sinE_) + 2.0 * vdot_k * (C_rs * cos2Phi - C_rc * sin2Phi);
            Float OMEGAdot_k = OMEGA_DOT - WGS84_OMEGA<Float>;

            // velocity
            Float xdot_k_prime =
                    rdot_k * cosu - r_k * udot_k * sinu;  // x-velocity in orbital frame
            Float ydot_k_prime =
                    rdot_k * sinu + r_k * udot_k * cosu;  // y-velocity in orbital frame
            v_eb_e_(1) = -(x_k_prime * OMEGAdot_k * sinOMEGA) + (xdot_k_prime * cosOMEGA) -
                         (ydot_k_prime * sinOMEGA * cosi) -
                         (y_k_prime * (OMEGAdot_k * cosOMEGA * cosi - idot_k * sinOMEGA * sini));
            v_eb_e_(2) = (x_k_prime * OMEGAdot_k * cosOMEGA) + (xdot_k_prime * sinOMEGA) +
                         (ydot_k_prime * cosOMEGA * cosi) -
                         (y_k_prime * (OMEGAdot_k * sinOMEGA * cosi + idot_k * cosOMEGA * sini));
            v_eb_e_(3) = (ydot_k_prime * sini) + (y_k_prime * idot_k * cosi);

            if constexpr (CalcAccel) {
                Float tmp = WGS84_R0<Float> / r_k;
                Float r2 = r_k * r_k;
                Float r3 = r2 * r_k;
                Float F = -1.5 * J2<Float> * (GM<Float> / r2) * (tmp * tmp);
                Float OMEGA_DOT_E2 = WGS84_OMEGA<Float> * WGS84_OMEGA<Float>;
                tmp = 5.0 * std::pow(r_eb_e_(2) / r_k, 2);

                a_eb_e_(1) = -GM<Float> * (r_eb_e_(0) / r3) +
                             F * ((1.0 - tmp) * (r_eb_e_(0) / r_k)) +
                             (2.0 * v_eb_e_(1) * WGS84_OMEGA<Float>)+(r_eb_e_(0) * OMEGA_DOT_E2);
                a_eb_e_(2) = -GM<Float> * (r_eb_e_(1) / r3) +
                             F * ((1.0 - tmp) * (r_eb_e_(1) / r_k)) -
                             (2.0 * v_eb_e_(0) * WGS84_OMEGA<Float>)+(r_eb_e_(1) * OMEGA_DOT_E2);
                a_eb_e_(3) =
                        -GM<Float> * (r_eb_e_(2) / r3) + F * ((3.0 - tmp) * (r_eb_e_(2) / r_k));
            }
        }
    }
};

//* ===== SGP4 Ephemeris ======================================================================= *//

template <typename Float = double>
class SGP4Ephemeris : public Ephemeris<Float> {
  public:
    SGP4Ephemeris() = default;
    // ~SGP4Ephemeris() = default;

    int catalog_id{0};   // Satellite catalog number
    int week{0};         // GPS week number of the TLE epoch
    Float t_oe{0.0};     // Second of GPS week of the TLE epoch
    Float n_dot{0.0};    // (XNDT20) 1st derivative of mean motion
    Float n_ddot{0.0};   // (XNDD60) 2nd derivative of mean motion
    Float B_star{0.0};   // (BSTAR) drag/radian pressure coefficient
    Float i_0{0.0};      // (XINCL) inclination angle
    Float OMEGA_0{0.0};  // (OMEGA0) right ascension of the ascending node
    Float e_0{0.0};      // (E0) eccentricity
    Float omega{0.0};    // (XNODE0) argument of perigee
    Float M_0{0.0};      // (XM0) mean anomaly
    Float n_0{0.0};      // (XN0) mean motion

    //! === P ===
    void P(const Float& tow, Vec3<Float>& p) {
        sgp4<false, false>(tow);
        p = r_eb_e_;
    }

    //! === PV ===
    void PV(const Float& tow, Vec3<Float>& p, Vec3<Float>& v) {
        sgp4<true, false>(tow);
        p = r_eb_e_;
        v = v_eb_e_;
    }

    //! === PVA ===
    void PVA(const Float& tow, Vec3<Float>& p, Vec3<Float>& v, Vec3<Float>& a) {
        sgp4<true, true>(tow);
        p = r_eb_e_;
        v = v_eb_e_;
        a = a_eb_e_;
    }

    //! === B ===
    void B(const Float& tow, Float& tR) {
        // RelativeClock<false, false>(tow);
        tR = tR_;
    }

    //! === BD ===
    void BD(const Float& tow, Float& tR, Float& tR_dot) {
        // RelativeClock<true, false>(tow);
        tR = tR_;
        tR_dot = tR_dot_;
    }

    //! === BDR ===
    void BDR(const Float& tow, Float& tR, Float& tR_dot, Float& tR_ddot) {
        // RelativeClock<true, true>(tow);
        tR = tR_;
        tR_dot = tR_dot_;
        tR_ddot = tR_ddot_;
    }

    //! === PRINT ==
    void print() {
        std::cout << "--- SGP4 Ephemerides ---" << std::endl
                  << "t_oe:    " << M_0 << std::endl
                  << "n_dot:   " << n_dot << std::endl
                  << "n_ddot:  " << n_ddot << std::endl
                  << "B_star:  " << B_star << std::endl
                  << "i_0:     " << i_0 << std::endl
                  << "OMEGA_0: " << OMEGA_0 << std::endl
                  << "e_0:     " << e_0 << std::endl
                  << "omega:   " << omega << std::endl
                  << "M_0:     " << M_0 << std::endl
                  << "n_0:     " << n_0 << std::endl
                  << "------------------------" << std::endl;
    }

  protected:
    Vec3<Float> r_eb_e_ = Vec3<Float>::Zero();  // ecef position
    Vec3<Float> v_eb_e_ = Vec3<Float>::Zero();  // ecef velocity
    Vec3<Float> a_eb_e_ = Vec3<Float>::Zero();  // ecef acceleration
    Float tR_{0.0};
    Float tR_dot_{0.0};
    Float tR_ddot_{0.0};
    bool initialized_{false};
    bool sgp4_simple_{false};
    Float COSIO_{0.0};
    Float SINIO_{0.0};
    Float SINMO_{0.0};
    Float DELMO_{0.0};
    Float X3THM1_{0.0};
    Float X1MTH2_{0.0};
    Float X7THM1_{0.0};
    Float AODP_{0.0};
    Float XNODP_{0.0};
    Float S4_{0.0};
    Float QOMS24_{0.0};
    Float TSI_{0.0};
    Float ETA_{0.0};
    Float XMDOT_{0.0};
    Float OMGDOT_{0.0};
    Float XNODOT_{0.0};
    Float XNODCF_{0.0};
    Float T2COF_{0.0};
    Float XMCOF_{0.0};
    Float XLCOF_{0.0};
    Float OMGCOF_{0.0};
    Float AYCOF_{0.0};
    Float C1_{0.0};
    Float C2_{0.0};
    Float C3_{0.0};
    Float C4_{0.0};
    Float C5_{0.0};

    //! === init ===
    void init() {
        // TLE specific constants
        Float E0SQ = e_0 * e_0;
        Float BETA02 = 1.0 - E0SQ;
        Float BETA0 = std::sqrt(BETA02);
        SINMO_ = std::sin(M_0);
        COSIO_ = std::cos(i_0);
        SINIO_ = std::sin(i_0);
        Float THETA2 = COSIO_ * COSIO_;
        Float THETA4 = THETA2 * THETA2;
        X1MTH2_ = -THETA2 + 1.0;
        X3THM1_ = -1.0 + 3.0 * THETA2;
        sgp4_simple_ = false;

        // recover original mean motion and semi-major axis from elements (Constant per TLE)
        Float a1 = std::pow((SGP_XKE<Float> / n_0), 2.0 / 3.0);
        Float del1 = 1.5 * SGP_CK2<Float> * X3THM1_ / (a1 * a1 * BETA0 * BETA02);
        Float a0 = a1 * (1.0 - del1 * (1.0 / 3.0 + del1 * (1.0 + 134.0 / 81.0 * del1)));
        Float del0 = 1.5 * SGP_CK2<Float> * X3THM1_ / (a0 * a0 * BETA0 * BETA02);
        XNODP_ = n_0 / (1.0 + del0);
        AODP_ = a0 / (1.0 - del0);

        // for perigee less than 220 km, equations are truncated to linear
        if ((AODP_ * (1.0 - e_0) / SGP_AE<Float>) < (220.0 / SGP_XKMPER<Float> + SGP_AE<Float>)) {
            sgp4_simple_ = true;
        }

        // initialize
        S4_ = SGP_S<Float>;
        QOMS24_ = SGP_QOMS2T<Float>;
        Float perigee = (AODP_ * (1.0 - e_0) - SGP_AE<Float>)*SGP_XKMPER<Float>;

        if (perigee < 156.0) {
            // For perigee between 98-156 km, the value of the constant s used in SGP4 is:
            S4_ = perigee - 78.0;

            if (perigee <= 98.0) {
                // For perigee below 98 km, the value of s is:
                S4_ = 20.0;
            }

            // If s is changed, (q0 - s*)^4 becomes:
            QOMS24_ = std::pow((120.0 - S4_) * SGP_AE<Float> / SGP_XKMPER<Float>, 4.0);
            S4_ /= (SGP_XKMPER<Float> + SGP_AE<Float>);
        }

        // constants given appropriate values of s_star and (q0 - s_star)^4
        Float PINVSQ = 1.0 / (AODP_ * AODP_ * BETA02 * BETA02);
        TSI_ = 1.0 / (AODP_ - S4_);
        ETA_ = AODP_ * e_0 * TSI_;
        Float ETASQ = ETA_ * ETA_;
        Float EETA = e_0 * ETA_;
        Float PSISQ = std::fabs(1.0 - ETASQ);
        Float COEF = QOMS24_ * std::pow(TSI_, 4.0);
        Float COEF1 = COEF / std::pow(PSISQ, 3.5);

        C2_ = COEF1 * XNODP_ *
              (AODP_ * (1.0 + 1.5 * ETASQ + EETA * (4.0 + ETASQ)) +
               0.75 * SGP_CK2<Float> * TSI_ / PSISQ * X3THM1_ *
                       (8.0 + 3.0 * ETASQ * (8.0 + ETASQ)));
        C1_ = B_star * C2_;
        C3_ = COEF * TSI_ * SGP_A3OVK2<Float> * XNODP_ * SGP_AE<Float> * SINIO_ / e_0;
        C4_ = 2.0 * XNODP_ * COEF1 * AODP_ * BETA02 *
              (ETA_ * (2.0 + 0.5 * ETASQ) + e_0 * (0.5 + 2.0 * ETASQ) -
               2.0 * SGP_CK2<Float> * TSI_ / (AODP_ * PSISQ) *
                       (-3.0 * X3THM1_ * (1.0 - 2.0 * EETA + ETASQ * (1.5 + 0.5 * EETA)) +
                        0.75 * X1MTH2_ * (2.0 * ETASQ - EETA * (1.0 + ETASQ)) *
                                std::cos(2.0 * OMEGA_0)));
        C5_ = 2.0 * COEF1 * AODP_ * BETA02 * (1.0 + 2.75 * (ETASQ + EETA) + EETA * ETASQ);

        Float TEMP1 = 3.0 * SGP_CK2<Float> * PINVSQ * XNODP_;
        Float TEMP2 = TEMP1 * SGP_CK2<Float> * PINVSQ;
        Float TEMP3 = 1.25 * SGP_CK4<Float> * PINVSQ * PINVSQ * XNODP_;
        XMDOT_ = XNODP_ + 0.5 * TEMP1 * BETA0 * X3THM1_ +
                 0.0625 * TEMP2 * BETA0 * (13.0 - 78.0 * THETA2 + 137.0 * THETA4);
        Float X1M5TH = -5.0 * THETA2 + 1.0;
        OMGDOT_ = 0.5 * TEMP1 * X1M5TH + 0.0625 * TEMP2 * (7.0 - 114.0 * THETA2 + 395.0 * THETA4);
        Float XHDOT1 = TEMP1 * COSIO_;
        XNODOT_ =
                XHDOT1 +
                (0.5 * TEMP2 * (4.0 - 19.0 * THETA2) + 2.0 * TEMP3 * (3.0 - 7.0 * THETA2)) * COSIO_;
        OMGCOF_ = B_star * C3_ * std::cos(OMEGA_0);
        XMCOF_ = -2.0 / 3.0 * COEF * B_star * SGP_AE<Float> / EETA;
        XNODCF_ = 3.5 * BETA02 * XHDOT1 * C1_;
        T2COF_ = 1.5 * C1_;
        XLCOF_ = 0.125 * SGP_A3OVK2<Float> * SINIO_ * (3.0 + 5.0 * COSIO_) / (1.0 + COSIO_);
        AYCOF_ = 0.25 * SGP_A3OVK2<Float> * SINIO_;
        DELMO_ = std::pow(1.0 + ETA_ * std::cos(M_0), 3.0);
        X7THM1_ = 7.0 * THETA2 - 1.0;
    }

    //! === SGP4 ===
    template <bool CalcVel, bool CalcAccel>
    void sgp4(const Float& tow) {
        // https://apps.dtic.mil/sti/pdfs/ADA093554.pdf
        // https://celestrak.org/NORAD/documentation/spacetrk.pdf

        if (!initialized_) {
            init();
            initialized_ = true;
        }

        Float TEMP, TEMP1, TEMP2, TEMP3, TEMP4, TEMP5, TEMP6;
        Float dt = CheckTime(tow - t_oe);

        // update for secular gravity and atmospheric drag
        Float XMDF = M_0 + XMDOT_ * dt;
        Float OMGADF = OMEGA_0 + OMGDOT_ * dt;
        Float XNODDF = omega + XNODOT_ * dt;
        Float OMEGA = OMGADF;
        Float XMP = XMDF;
        Float TSQ = dt * dt;
        Float XNODE = XNODDF + XNODCF_ * TSQ;
        Float TEMPA = 1.0 - C1_ * dt;
        Float TEMPE = B_star * C4_ * dt;
        Float TEMPL = T2COF_ * TSQ;
        if (!sgp4_simple_) {
            Float C1SQ = C1_ * C1_;
            Float D2 = 4.0 * AODP_ * TSI_ * C1SQ;
            TEMP = D2 * TSI_ * C1_ / 3.0;
            Float D3 = (17.0 * AODP_ + S4_) * TEMP;
            Float D4 = 0.5 * TEMP * AODP_ * TSI_ * (221.0 * AODP_ + 31.0 * S4_) * C1_;
            Float T3COF = D2 + 2.0 * C1SQ;
            Float T4COF = 0.25 * (3.0 * D3 + C1_ * (12.0 * D2 + 10.0 * C1SQ));
            Float T5COF = 0.2 * (3.0 * D4 + 12.0 * C1_ * D3 + 6.0 * D2 * D2 +
                                 15.0 * C1SQ * (2.0 * D2 + C1SQ));
            Float DELOMG = OMGCOF_ * dt;
            Float DELM = XMCOF_ * (std::pow(1.0 + ETA_ * std::cos(XMDF), 3) - DELMO_);
            TEMP = DELOMG + DELM;
            XMP = XMDF + TEMP;
            OMEGA = OMGADF - TEMP;
            Float TCUBE = TSQ * dt;
            Float TFOUR = dt * TCUBE;
            TEMPA -= (D2 * TSQ - D3 * TCUBE - D4 * TFOUR);
            TEMPE -= (TEMPE + B_star * C5_ * (std::sin(XMP) - SINMO_));
            TEMPL -= (TEMPL + T3COF * TCUBE + TFOUR * (T4COF + dt * T5COF));
        }
        Float A = AODP_ * TEMPA * TEMPA;
        Float E = e_0 - TEMPE;
        Float XL = XMP + OMEGA + XNODE + XNODP_ * TEMPL;
        Float BETA = std::sqrt(1.0 - E * E);
        Float XN = SGP_XKE<Float> / std::pow(A, 1.5);

        // long period periodics
        Float AXN = E * std::cos(OMEGA);
        TEMP = 1.0 / (A * BETA * BETA);
        Float XLL = TEMP * XLCOF_ * AXN;
        Float AYNL = TEMP * AYCOF_;
        Float XLT = XL + XLL;
        Float AYN = E * std::sin(OMEGA) + AYNL;

        // solve keplers equation
        Float CAPU = std::fmod(XLT - XNODE, TWO_PI<Float>);
        TEMP2 = CAPU;
        Float COSPW, SINPW, EPW;
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
        Float ECOSE = TEMP5 + TEMP6;
        Float ESINE = TEMP3 - TEMP4;
        Float ELSQ = AXN * AXN + AYN * AYN;
        TEMP = 1.0 - ELSQ;
        Float PL = A * TEMP;
        Float R = A * (1.0 - ECOSE);
        TEMP1 = 1.0 / R;
        TEMP2 = A * TEMP1;
        Float BETAL = std::sqrt(TEMP);
        TEMP3 = 1.0 / (1.0 + BETAL);
        Float COSU = TEMP2 * (COSPW - AXN + AYN * ESINE * TEMP3);
        Float SINU = TEMP2 * (SINPW - AYN - AXN * ESINE * TEMP3);
        Float U = std::atan2(SINU, COSU);
        Float SIN2U = 2.0 * SINU * COSU;
        Float COS2U = 2.0 * COSU * COSU - 1.0;
        TEMP = 1.0 / PL;
        TEMP1 = SGP_CK2<Float> * TEMP;
        TEMP2 = TEMP1 * TEMP;

        // update for short periodics
        Float RK = R * (1.0 - 1.5 * TEMP2 * BETAL * X3THM1_) + 0.5 * TEMP1 * X1MTH2_ * COS2U;
        Float UK = U - 0.25 * TEMP2 * X7THM1_ * SIN2U;
        Float XNODEK = XNODE + 1.5 * TEMP2 * COSIO_ * SIN2U;
        Float XINCK = i_0 + 1.5 * TEMP2 * COSIO_ * SINIO_ * COS2U;

        // orientation vectors
        Float SINUK = std::sin(UK);
        Float COSUK = std::cos(UK);
        Float SINIK = std::sin(XINCK);
        Float COSIK = std::cos(XINCK);
        Float SINNOK = std::sin(XNODEK);
        Float COSNOK = std::cos(XNODEK);

        Float XMX = -SINNOK * COSIK;
        Float XMY = COSNOK * COSIK;
        Float UX = XMX * SINUK + COSNOK * COSUK;
        Float UY = XMY * SINUK + SINNOK * COSUK;
        Float UZ = SINIK * SINUK;

        // position
        r_eb_e_(0) = RK * UX;
        r_eb_e_(1) = RK * UY;
        r_eb_e_(2) = RK * UZ;
        r_eb_e_ *= (1000.0 * SGP_XKMPER<Float> / SGP_AE<Float>);  // [earth-radii] -> [m]

        // velocity
        if constexpr (CalcVel) {
            // short period preliminary quantities
            Float RDOT = SGP_XKE<Float> * std::sqrt(A) * ESINE * TEMP1;
            Float RFDOT = SGP_XKE<Float> * std::sqrt(PL) * TEMP1;
            // update for short periodics
            Float RDOTK = RDOT - XN * TEMP1 * X1MTH2_ * SIN2U;
            Float RFDOTK = RFDOT + XN * TEMP1 * (X1MTH2_ * COS2U + 1.5 * X3THM1_);
            // orientation vectors
            Float VX = XMX * COSUK - COSNOK * SINUK;
            Float VY = XMY * COSUK - SINNOK * SINUK;
            Float VZ = SINIK * COSUK;

            v_eb_e_(0) = RDOTK * UX + RFDOTK * VX;
            v_eb_e_(1) = RDOTK * UY + RFDOTK * VY;
            v_eb_e_(2) = RDOTK * UZ + RFDOTK * VZ;
            v_eb_e_ *=
                    (1000.0 * SGP_XKMPER<Float> / SGP_AE<Float> /
                     60.0);  // [earth-radii/min] -> [m/s]

            // acceleration
            if constexpr (CalcAccel) {
            }
        }
    }
};

}  // namespace satutils

#endif