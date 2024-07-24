#ifndef NAVENGINE_INCLUDE_SATUTILS_GPS_EPHEMERIS_HPP
#define NAVENGINE_INCLUDE_SATUTILS_GPS_EPHEMERIS_HPP

#include <cstdint>
#include <cmath>
#include <numbers>

#include <Eigen/Dense>

#include <navtools/math.hpp>
#include <satutils/constants.hpp>

namespace satutils {

// struct GpsClockDataBits
// {
//   int8_t T_GD {0};
//   uint16_t t_oc {0};
//   int32_t a_f0 {0};
//   int16_t a_f1 {0};
//   int8_t a_f2 {0};
//   uint16_t IODC {0};
// };

template<typename Float>
struct GpsClockData
{
  Float T_GD {0.0};
  Float t_oc {0.0};
  Float a_f0 {0.0};
  Float a_f1 {0.0};
  Float a_f2 {0.0};
  uint16_t IODC {0};

  Float Offset(const Float gps_time) const
  {
    Float dt = gps_time - t_oc;
    return a_f0 + (a_f1 * dt) + (a_f2 * dt * dt);
  }
  
  Float OffsetRate(const Float gps_time) const
  {
    return a_f1 + (2.0 * (gps_time - t_oc) * a_f2);
  }
  
  Float OffsetRateRate() const
  {
    return 2.0 * a_f2;
  }

  struct Storage
  {
    Float T_GD {0.0};
    Float t_oc {0.0};
    Float a_f0 {0.0};
    Float a_f1 {0.0};
    Float a_f2 {0.0};
    uint16_t IODC {0};
  };

  static constexpr Storage LowerLimits =
  {
    -std::pow(2.0,7-31),
    0,
    -std::pow(2.0,21-31),
    -std::pow(2.0,15-43),
    -std::pow(2.0,7-55),
    0 
  };

  static constexpr Storage UpperLimits = 
  {
    127 * std::pow(2.0,-31),
    604784,
    (std::pow(2.0,21) - 1.0) * std::pow(2.0,-31),
    (std::pow(2.0,15) - 1.0) * std::pow(2.0,-43),
    127 * std::pow(2.0,-55),
    1023
  };

  static constexpr Storage ScaleFactors =
  {
    std::pow(2.0,-31),
    16,
    std::pow(2.0,-31),
    std::pow(2.0,-43),
    std::pow(2.0,-55),
    1
  };

  // TODO
  // void Randomize();
};

struct KeplerEphemeris
{
  double M_0 {0.0};
  double del_n {0.0};
  double e {0.0};
  double sqrtA {0.0};
  double Omega_0 {0.0};
  double i_0 {0.0};
  double omega {0.0};
  double Omega_dot {0.0};
  double IDOT {0.0};
  
  double C_uc {0.0};
  double C_us {0.0};
  double C_rc {0.0};
  double C_rs {0.0};
  double C_ic {0.0};
  double C_is {0.0};

  double t_oe {0.0};
  uint8_t IODE {0};

  constexpr static double J2 = 0.0010826262;
  constexpr static double RELETIVISTIC_F = -4.442807633e-10;
  constexpr static double WGS84_MU = 3.986005e14;
  constexpr static double WGS84_EARTH_RATE = 7.2921151467e-5;
  constexpr static double WGS84_EQUAT_RADIUS = 6378137.0;

  inline double EfromAnomaly(const double M_k, const unsigned int iterations) const
  {
    double E_k = M_k;
    double del_E;
    do {
      del_E = (M_k - E_k + (e * std::sin(E_k))) / (1.0 - (e * std::cos(E_k)));
      E_k += del_E;
    }
    while (del_E > 1.0e-15);
    return E_k;
  }

  inline double EfromTime(const double gps_time, const unsigned int iterations) const
  {
    double A = std::pow(sqrtA,2.0);

    double n_0 = std::sqrt(WGS84_MU / std::pow(A,3.0));
    double t_k = gps_time - t_oe;
    if (t_k > 302400.0) {
      t_k -= 604800.0;
    } else if (t_k < -302400) {
      t_k += 604800.0;
    }

    double n = n_0 + del_n;
    return EfromAnomaly(M_0 + (n * t_k),5);
  }

  template<bool CalcVel, bool CalcAccel>
  void CalcPVA(const double gps_time, Eigen::Vector3d& pos, Eigen::Vector3d& vel,
    Eigen::Vector3d& accel) const
  {
    double A = std::pow(sqrtA,2.0);
    
    double n_0 = std::sqrt(WGS84_MU / std::pow(A,3.0));
    double t_k = gps_time - t_oe;
    if (t_k > 302400.0) {
      t_k -= 604800.0;
    } else if (t_k < -302400) {
      t_k += 604800.0;
    }
    
    double n = n_0 + del_n;
    double M_k = M_0 + (n * t_k);
    M_k = std::fmod(M_k + GPS_2PI<double>, GPS_2PI<double>);

    double E_k = EfromAnomaly(M_k,5);
    E_k = std::fmod(E_k + GPS_2PI<double>, GPS_2PI<double>);

    // double v_k = std::sqrt((1.0+e)/(1.0-e)) * std::tan(E_k / 2.0);
    // v_k = 2.0 * std::atan(v_k);
    double cv_k = (std::cos(E_k) - e) / (1.0 - (e * std::cos(E_k)));
    double sv_k = (std::sqrt(1.0 - (e*e)) * std::sin(E_k)) / (1.0 - (e * std::cos(E_k)));
    double v_k = std::atan2(sv_k,cv_k);

    double Phi_k = v_k + omega;
    navtools::CircMod(Phi_k,GPS_2PI<double>);

    double Phi2 = 2.0 * Phi_k;
    double du_k = (C_us * std::sin(Phi2)) + (C_uc * std::cos(Phi2));
    double dr_k = (C_rs * std::sin(Phi2)) + (C_rc * std::cos(Phi2));
    double di_k = (C_is * std::sin(Phi2)) + (C_ic * std::cos(Phi2));
    
    double u_k = Phi_k + du_k;
    double r_k = A * (1.0 - (e * std::cos(E_k))) + dr_k;
    double i_k = i_0 + di_k + (IDOT * t_k);

    double x_orb = r_k * std::cos(u_k);
    double y_orb = r_k * std::sin(u_k);

    double Omega_k = Omega_0 + ((Omega_dot - WGS84_EARTH_RATE) * t_k) - (WGS84_EARTH_RATE * t_oe);

    pos(0) = (x_orb * std::cos(Omega_k)) - (y_orb * std::cos(i_k) * std::sin(Omega_k));
    pos(1) = (x_orb * std::sin(Omega_k)) + (y_orb * std::cos(i_k) * std::cos(Omega_k));
    pos(2) = y_orb * std::sin(i_k);

    // Velocity Terms
    if constexpr (CalcVel) {
      double Ed_k = n / (1.0 - (e * std::cos(E_k)));
      double vd_k = Ed_k * std::sqrt(1.0 - (e*e)) / (1.0 - (e * std::cos(E_k)));
      double id_k = IDOT + ( 2.0 * vd_k * ((C_is * std::cos(Phi2)) - (C_ic * std::sin(Phi2))) );
      double ud_k = vd_k + ( 2.0 * vd_k * ((C_us * std::cos(Phi2)) - (C_uc * std::sin(Phi2))) );
      double rd_k = (e * A * Ed_k * std::sin(E_k))
                    + ( 2.0 * vd_k * ((C_rs * std::cos(Phi2)) - (C_rc * std::sin(Phi2))) );
      double Omega_dot_k = Omega_dot - WGS84_EARTH_RATE;
      
      double xd_orb = (rd_k * std::cos(u_k)) - (r_k * ud_k * std::sin(u_k));
      double yd_orb = (rd_k * std::sin(u_k)) + (r_k * ud_k * std::cos(u_k));
      
      vel(0) = (-x_orb * Omega_dot_k * std::sin(Omega_k))
              + (xd_orb * std::cos(Omega_k))
              - (yd_orb * std::sin(Omega_k) * std::cos(i_k))
              - ( y_orb * ((Omega_dot_k * std::cos(Omega_k) * std::cos(i_k)) 
                        - (id_k * std::sin(Omega_k * std::sin(i_k)))) ); 
    
      vel(1) = (x_orb * Omega_dot_k * std::cos(Omega_k))
              + (xd_orb * std::sin(Omega_k))
              + (yd_orb * std::cos(Omega_k) * std::cos(i_k))
              - ( y_orb * ((Omega_dot_k * std::sin(Omega_k) * std::cos(i_k)) 
                        + (id_k * std::cos(Omega_k) * std::sin(i_k))) ); 
    
      vel(2) = (yd_orb * std::sin(i_k)) + (y_orb * id_k * std::cos(i_k));
      
    }

    // Acceleration Terms
    if constexpr (CalcAccel) {  
      double r2 = r_k * r_k;
      double r3 = r2 * r_k;
      double F = -1.5 * J2 * (WGS84_MU / r2) * std::pow(WGS84_EQUAT_RADIUS / r_k, 2.0);
      double F_term = F * ( 1.0 - (5.0 * std::pow(pos(2) / r_k, 2.0)) );  
      double omega_e2 = std::pow(WGS84_EARTH_RATE, 2.0);

      accel(0) = (-WGS84_MU * pos(0) / r3) + (F_term * pos(0) / r_k)
                + (2.0 * vel(1) * WGS84_EARTH_RATE) + (pos(0) * omega_e2);
    
      accel(1) = (-WGS84_MU * pos(1) / r3) + (F_term * pos(1) / r_k)
                - (2.0 * vel(0) * WGS84_EARTH_RATE) + (pos(1) * omega_e2);
    
      accel(2) = (-WGS84_MU * pos(2) / r3)
                + ( F * (3.0 - (5.0 * std::pow(pos(2) / r_k, 2.0))) * pos(2) / r_k );  
    }
  }

  void P(const double gps_time, Eigen::Vector3d& pos) const
  {
    Eigen::Vector3d filler;
    CalcPVA<false,false>(gps_time, pos, filler, filler);
  }

  void PV(const double gps_time, Eigen::Vector3d& pos, Eigen::Vector3d& vel) const
  {
    Eigen::Vector3d filler;
    CalcPVA<true,false>(gps_time, pos, vel, filler);
  }

  void PVA(const double gps_time, Eigen::Vector3d& pos, Eigen::Vector3d& vel,
    Eigen::Vector3d& accel) const
  {
    CalcPVA<true,true>(gps_time, pos, vel, accel);
  }

  double RelTime(const double gps_time) const
  {
    return RELETIVISTIC_F * std::numbers::e * sqrtA * std::sin(EfromTime(gps_time,5));
  }

  double RelTimeRate(const double gps_time) const
  {
    double A = std::pow(sqrtA,2.0);
    double n_0 = std::sqrt(WGS84_MU / std::pow(A,3.0));
    double t_k = gps_time - t_oe;
    if (t_k > 302400.0) {
      t_k -= 604800.0;
    } else if (t_k < -302400) {
      t_k += 604800.0;
    }  
    double n = n_0 + del_n;
    double e_cos_E = std::numbers::e * std::cos( EfromAnomaly(M_0 + (n * t_k),5) );

    return (n * RELETIVISTIC_F * sqrtA * e_cos_E) / (1.0 - e_cos_E);
  }

  double RelTimeRateRate(const double gps_time) const
  {
    double A = std::pow(sqrtA,2.0);
    double n_0 = std::sqrt(WGS84_MU / std::pow(A,3.0));
    double t_k = gps_time - t_oe;
    if (t_k > 302400.0) {
      t_k -= 604800.0;
    } else if (t_k < -302400) {
      t_k += 604800.0;
    }  
    double n = n_0 + del_n;
    double E_k = EfromAnomaly(M_0 + (n * t_k),5);

    return ( n * n * RELETIVISTIC_F * std::numbers::e * sqrtA * std::sin(E_k) )
            / std::pow(1.0 - std::numbers::e * std::cos(E_k), 2);
  }

  void Randomize(); // TODO

  void Print() const
  {
    std::cout << "M_0:       " << M_0 << '\n';
    std::cout << "del_n:     " << del_n << '\n';
    std::cout << "e:         " << e << '\n';
    std::cout << "sqrtA:     " << sqrtA << '\n';
    std::cout << "Omega_0:   " << Omega_0 << '\n';
    std::cout << "i_0:       " << i_0 << '\n';
    std::cout << "omega:     " << omega << '\n';
    std::cout << "Omega_dot: " << Omega_dot << '\n';
    std::cout << "IDOT:      " << IDOT << '\n';
    std::cout << "C_uc:      " << C_uc << '\n';
    std::cout << "C_us:      " << C_us << '\n';
    std::cout << "C_rc:      " << C_rc << '\n';
    std::cout << "C_rs:      " << C_rs << '\n';
    std::cout << "C_ic:      " << C_ic << '\n';
    std::cout << "C_is:      " << C_is << '\n';
    std::cout << "t_oe:      " << t_oe << '\n';
    std::cout << "IODE:      " << static_cast<int>(IODE) << '\n';
  }
};

constexpr KeplerEphemeris GpsEphemerisScaleFactors =
{
  std::pow(2.0,-31),
  std::pow(2.0,-43),
  std::pow(2.0,-33),
  std::pow(2.0,-19),
  std::pow(2.0,-31),
  std::pow(2.0,-31),
  std::pow(2.0,-31),
  std::pow(2.0,-43),
  std::pow(2.0,-43),
  std::pow(2.0,-29),
  std::pow(2.0,-29),
  std::pow(2.0,-5),
  std::pow(2.0,-5),
  std::pow(2.0,-29),
  std::pow(2.0,-29),
  16,
  1
};

constexpr KeplerEphemeris GpsEphemerisLowerLimits = 
{
  -std::pow(2.0,31-31),
  -std::pow(2.0,15-43),
  0.0,
  2530.0,
  -std::pow(2.0,31-31),
  -std::pow(2.0,31-31),
  -std::pow(2.0,31-31),
  -6.33e-7,
  -std::pow(2.0,13-43), // IDOT
  -std::pow(2.0,15-29), // C_uc
  -std::pow(2.0,15-29), // C_us
  -std::pow(2.0,15-5), // C_rc
  -std::pow(2.0,15-5), // C_rs
  -std::pow(2.0,15-29), // C_ic
  -std::pow(2.0,15-29), // C_is
  0,
  0
};

constexpr KeplerEphemeris GpsEphemerisUpperLimits = 
{
  (std::pow(2.0,31) - 1.0) * std::pow(2.0,-31),
  (std::pow(2.0,15) - 1.0) * std::pow(2.0,-43),
  0.03,
  8192.0,
  (std::pow(2.0,31) - 1.0) * std::pow(2.0,-31),
  (std::pow(2.0,31) - 1.0) * std::pow(2.0,-31),
  (std::pow(2.0,31) - 1.0) * std::pow(2.0,-31),
  0.0,
  (std::pow(2.0,13) - 1.0) * std::pow(2.0,-43), // IDOT
  (std::pow(2.0,15) - 1.0) * std::pow(2.0,-29), // C_uc
  (std::pow(2.0,15) - 1.0) * std::pow(2.0,-29), // C_us
  (std::pow(2.0,15) - 1.0) * std::pow(2.0,-5), // C_rc
  (std::pow(2.0,15) - 1.0) * std::pow(2.0,-5), // C_rs
  (std::pow(2.0,15) - 1.0) * std::pow(2.0,-29), // C_ic
  (std::pow(2.0,15) - 1.0) * std::pow(2.0,-29), // C_is
  604784,
  255
};

} // namespace satutils
#endif
