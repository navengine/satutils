#ifndef SATUTILS_INCLUDE_GPS_EPHEMERIS_HPP
#define SATUTILS_INCLUDE_GPS_EPHEMERIS_HPP

#include <cstdint>

#include <navtools/binary-ops.hpp>
#include <satutils/ephemeris.hpp>
#include "ephemeris.hpp"

namespace satutils {

// CAUTION: max val allowable for Pow depends on the bit size of intmax_t. It should be at least 64 on most systems.
template<int Pow, typename Float = double>
constexpr Float PowerOfTwo()
{
  if constexpr (Pow == 0) return Float(1);
  else if constexpr (Pow < 0) {
    return Float(1) / Float(1 << -Pow);
  }
  else {
    return Float(intmax_t(1) << Pow);
  }
}


template<typename Float>
constexpr uint32_t ParamToBinary(const Float param, const Float scale_factor)
{
  Float quotient = param / scale_factor;
  return static_cast<uint32_t>(quotient < 0 ? quotient - 0.5 : quotient + 0.5);
}


double ParamFromBinary(uint32_t data, const double scale_factor, const uint8_t num_bits,
  const bool twos_comp)
{
  assert((num_bits >= 1) && (num_bits <= 32));
  assert(scale_factor > 0.0);

  static const uint32_t one = 1; // just avoiding the differing default size of integer literals

  if (num_bits != 32)
    data &= ((one << num_bits) - 1); // sets all MSBs past (num_bits-1) position to zero
  if (twos_comp && (bitVal<true>(data, num_bits-1))) {
    uint32_t signed_bit = (one << (num_bits-1));
    data &= (~signed_bit);
    data = (signed_bit - data); // un-signs the integer
    return -scale_factor * static_cast<double>(data);
  } else {
    return scale_factor * static_cast<double>(data);
  }
}


struct GpsL1Ephemeris
{
  int32_t M_0 {0};
  int16_t delta_n {0};
  uint32_t e {0};
  uint32_t sqrtA {0};
  int32_t OMEGA_0 {0};
  int32_t i_0 {0};
  int32_t omega {0};
  int32_t OMEGA_DOT {0};
  int16_t IDOT {0};
  int16_t C_uc {0};
  int16_t C_us {0};
  int16_t C_rc {0};
  int16_t C_rs {0};
  int16_t C_ic {0};
  int16_t C_is {0};
  uint16_t t_oe {0};
  uint8_t IODE {0};

  template<typename T>
  struct Info
  {
    T M_0 {0};
    T delta_n {0};
    T e {0};
    T sqrtA {0};
    T OMEGA_0 {0};
    T i_0 {0};
    T omega {0};
    T OMEGA_DOT {0};
    T IDOT {0};
    T C_uc {0};
    T C_us {0};
    T C_rc {0};
    T C_rs {0};
    T C_ic {0};
    T C_is {0};
    T t_oe {0};
    T IODE {0};
  };

  template<typename Float>
  static constexpr Info<Float> scale_factors = 
  {
    PowerOfTwo<-31,Float>(), /* M_0 */ 
    PowerOfTwo<-43,Float>(), /* delta_n */ 
    PowerOfTwo<-33,Float>(), /* e */ 
    PowerOfTwo<-19,Float>(), /* sqrtA */ 
    PowerOfTwo<-31,Float>(), /* OMEGA_0 */ 
    PowerOfTwo<-31,Float>(), /* i_0 */ 
    PowerOfTwo<-31,Float>(), /* omega */ 
    PowerOfTwo<-43,Float>(), /* OMEGA_DOT */ 
    PowerOfTwo<-43,Float>(), /* IDOT */ 
    PowerOfTwo<-29,Float>(), /* C_uc */ 
    PowerOfTwo<-29,Float>(), /* C_us */ 
    PowerOfTwo<-5,Float>(), /* C_rc */ 
    PowerOfTwo<-5,Float>(), /* C_rs */ 
    PowerOfTwo<-29,Float>(), /* C_ic */ 
    PowerOfTwo<-29,Float>(), /* C_is */ 
    PowerOfTwo<4,Float>(), /* t_oe */ 
    PowerOfTwo<0,Float>()    /* IODE */ 
  };

  static constexpr Info<uint8_t> num_bits =
  {
    .M_0 = 32,
    .delta_n = 16,
    .e = 32,
    .sqrtA = 32,
    .OMEGA_0 = 32,
    .i_0 = 32,
    .omega = 32,
    .OMEGA_DOT = 24,
    .IDOT = 14,
    .C_uc = 16,
    .C_us = 16,
    .C_rc = 16,
    .C_rs = 16,
    .C_ic = 16,
    .C_is = 16,
    .t_oe = 16,
    .IODE = 8
  };

  static constexpr Info<bool> signage =
  {
    .M_0 = true,
    .delta_n = true,
    .e = false,
    .sqrtA = false,
    .OMEGA_0 = true,
    .i_0 = true,
    .omega = true,
    .OMEGA_DOT = true,
    .IDOT = true,
    .C_uc = true,
    .C_us = true,
    .C_rc = true,
    .C_rs = true,
    .C_ic = true,
    .C_is = true,
    .t_oe = false,
    .IODE = false
  };

  GpsL1Ephemeris()
  {}

  template<typename Float>
  void Store(const KeplerEphemeris<Float>& ephem)
  {
    M_0 = ParamToBinary(ephem.M_0, scale_factors<Float>.M_0);
    delta_n = ParamToBinary(ephem.delta_n, scale_factors<Float>.delta_n);
    e = ParamToBinary(ephem.e, scale_factors<Float>.e);
    sqrtA = ParamToBinary(ephem.sqrtA, scale_factors<Float>.sqrtA);
    OMEGA_0 = ParamToBinary(ephem.OMEGA_0, scale_factors<Float>.OMEGA_0);
    i_0 = ParamToBinary(ephem.i_0, scale_factors<Float>.i_0);
    omega = ParamToBinary(ephem.omega, scale_factors<Float>.omega);
    OMEGA_DOT = ParamToBinary(ephem.OMEGA_DOT, scale_factors<Float>.OMEGA_DOT);
    IDOT = ParamToBinary(ephem.IDOT, scale_factors<Float>.IDOT);
    C_uc = ParamToBinary(ephem.C_uc, scale_factors<Float>.C_uc);
    C_us = ParamToBinary(ephem.C_us, scale_factors<Float>.C_us);
    C_rc = ParamToBinary(ephem.C_rc, scale_factors<Float>.C_rc);
    C_rs = ParamToBinary(ephem.C_rs, scale_factors<Float>.C_rs);
    C_ic = ParamToBinary(ephem.C_ic, scale_factors<Float>.C_ic);
    C_is = ParamToBinary(ephem.C_is, scale_factors<Float>.C_is);
    t_oe = ParamToBinary(ephem.t_oe, scale_factors<Float>.t_oe);
    IODE = uint8_t(ephem.IODE);
  }

  template<typename Float>
  GpsL1Ephemeris(const KeplerEphemeris<Float>& ephem)
  {
    Store(ephem);
  }

  template<typename Float>
  Float ParamValue(uint8_t index) const
  {
    switch(index) {
      case 0:
        return ParamFromBinary(M_0, scale_factors<Float>.M_0, num_bits.M_0, signage.M_0);
      case 1:
        return ParamFromBinary(delta_n, scale_factors<Float>.delta_n, num_bits.delta_n, signage.delta_n);
      case 2:
        return ParamFromBinary(e, scale_factors<Float>.e, num_bits.e, signage.e);
      case 3:
        return ParamFromBinary(sqrtA, scale_factors<Float>.sqrtA, num_bits.sqrtA, signage.sqrtA);
      case 4:
        return ParamFromBinary(OMEGA_0, scale_factors<Float>.OMEGA_0, num_bits.OMEGA_0, signage.OMEGA_0);
      case 5:
        return ParamFromBinary(i_0, scale_factors<Float>.i_0, num_bits.i_0, signage.i_0);
      case 6:
        return ParamFromBinary(omega, scale_factors<Float>.omega, num_bits.omega, signage.omega);
      case 7:
        return ParamFromBinary(OMEGA_DOT, scale_factors<Float>.OMEGA_DOT, num_bits.OMEGA_DOT, signage.OMEGA_DOT);
      case 8:
        return ParamFromBinary(IDOT, scale_factors<Float>.IDOT, num_bits.IDOT, signage.IDOT);
      case 9:
        return ParamFromBinary(C_uc, scale_factors<Float>.C_uc, num_bits.C_uc, signage.C_uc);
      case 10:
        return ParamFromBinary(C_us, scale_factors<Float>.C_us, num_bits.C_us, signage.C_us);
      case 11:
        return ParamFromBinary(C_rc, scale_factors<Float>.C_rc, num_bits.C_rc, signage.C_rc);
      case 12:
        return ParamFromBinary(C_rs, scale_factors<Float>.C_rs, num_bits.C_rs, signage.C_rs);
      case 13:
        return ParamFromBinary(C_ic, scale_factors<Float>.C_ic, num_bits.C_ic, signage.C_ic);
      case 14:
        return ParamFromBinary(C_is, scale_factors<Float>.C_is, num_bits.C_is, signage.C_is);
      case 15:
        return ParamFromBinary(t_oe, scale_factors<Float>.t_oe, num_bits.t_oe, signage.t_oe);
      case 16:
        return static_cast<Float>(IODE);
      default:
        return 0;
    }
  }

  template<typename Float>
  CreateKepler(KeplerEphemeris<Float>& ephem) const
  {
    ephem.M_0 = ParamValue<Float>(0); 
    ephem.delta_n = ParamValue<Float>(1);
    ephem.e = ParamValue<Float>(2);
    ephem.sqrtA = ParamValue<Float>(3);
    ephem.OMEGA_0 = ParamValue<Float>(4);
    ephem.i_0 = ParamValue<Float>(5);
    ephem.omega = ParamValue<Float>(6);
    ephem.OMEGA_DOT = ParamValue<Float>(7);
    ephem.IDOT = ParamValue<Float>(8);
    ephem.C_uc = ParamValue<Float>(9);
    ephem.C_us = ParamValue<Float>(10);
    ephem.C_rc = ParamValue<Float>(11);
    ephem.C_rs = ParamValue<Float>(12);
    ephem.C_ic = ParamValue<Float>(13);
    ephem.C_is = ParamValue<Float>(14);
    ephem.t_oe = ParamValue<Float>(15);
    ephem.IODE = ParamValue<Float>(16);
  }

  template<typename Float>
  KeplerEphemeris<Float> CreateKepler() const
  {
    KeplerEphemeris<Float> result;
    CreateKepler(result);
    return result;
  }

};


} // namespace satutils
#endif
