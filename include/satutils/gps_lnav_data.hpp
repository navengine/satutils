#ifndef SATUTILS_INCLUDE_GPS_LNAV_DATA_HPP
#define SATUTILS_INCLUDE_GPS_LNAV_DATA_HPP

#include <cstdint>
#include <string>
#include <sstream>

#include <navtools/constants.hpp>
#include <navtools/binary-ops.hpp>
#include <satutils/ephemeris.hpp>


namespace satutils {


class LnavData
{
public:
  class Word
  {
  private:
    // MSB is aligned with first chronological bit
    uint32_t data_ {0};

  public:
    Word() {}
    Word(const uint32_t& val) : data_{val} {}
    explicit Word(const Word& word) : data_{word.data_} {}

    uint32_t& data()
    { return data_; }

    // start_idx and end_idx: 0 is MSB (first bit transmitted chronologically)
    uint32_t GetSegment(uint8_t start_idx, uint8_t end_idx) const
    {
      return navtools::GetBits<false>(data_, start_idx, end_idx);
    }
    
    /*
    Gets a segment of the word's data and stores the segment in a particular range of a given value
    -mainly used for getting parameter values out of words
    val_start_idx: 0 is LSB
    start_idx and end idx: 0 is first chronological bit (MSB)
    */
    void StoreSegment(uint32_t& val, uint8_t val_start_idx, uint8_t start_idx, uint8_t end_idx) const
    {
      static constexpr uint32_t one_32 = static_cast<uint32_t>(1);
      val &= 
           ~(((one_32 << (end_idx - start_idx + 1)) - one_32) << val_start_idx);
      val |= (GetSegment(start_idx, end_idx) << val_start_idx);
    }

    /*
    val_start and val_end: 0 is LSB
    start_word_idx: 0 is MSB
    */
    void SetSegment(uint32_t val, uint8_t start_word_idx, uint8_t val_start, uint8_t val_end)
    {
      static constexpr uint32_t one_32 = static_cast<uint32_t>(1);
      uint8_t end_idx = start_word_idx + (val_end - val_start);
      data_ &= 
           ~(((one_32 << (val_end - val_start + 1)) - one_32) << (31 - end_idx));
      data_ |= (navtools::GetBits<true>(val, val_start, val_end) << (31 - end_idx));
    }
    
    bool operator()(const uint8_t& i) const
    {
      assert((i >= 0) && (i < 30));
      return navtools::GetBit<false>(data_,i);
    }

    void SetBit(const uint8_t& i, const bool val)
    {
      assert((i >= 0) && (i < 30));
      navtools::SetBitTo(data_,i,val);
    }

    std::string str() const
    {
      std::stringstream stream;
      for (int i = 0; i < 30; i++) {
        stream << navtools::GetBit<false>(data_,i);
      }
      return stream.str();
    }
  };

  template<typename T>
  static uint32_t ParamToBinary(const T param, const T scale_factor)
  {
    T quotient = param / scale_factor;
    return static_cast<uint32_t>(quotient < 0 ? quotient - 0.5 : quotient + 0.5);
  9

  template<typename T>
  static T ParamFromBinary(uint8_t data, const T scale_factor,
    const uint8_t num_bits, const bool twos_comp)
  {
    assert((num_bits >= 1) && (num_bits <= 8));
    assert(scale_factor > T(0));

    static const uint8_t one = 1; // just avoiding the differing default size of integer literals

    if (num_bits != 8)
      data &= ((one << num_bits) - 1); // sets all MSBs past (num_bits-1) position to zero
    if (twos_comp && (navtools::GetBit<true>(data, num_bits-1))) {
      uint32_t signed_bit = (one << (num_bits-1));
      data &= (~signed_bit);
      data = (signed_bit - data); // un-signs the integer
      return -scale_factor * static_cast<T>(data);
    } else {
      return scale_factor * static_cast<T>(data);
    }
  }
  
  template<typename T>
  static T ParamFromBinary(uint16_t data, const T scale_factor,
    const uint8_t num_bits, const bool twos_comp)
  {
    assert((num_bits >= 1) && (num_bits <= 16));
    assert(scale_factor > T(0));

    static const uint16_t one = 1; // just avoiding the differing default size of integer literals

    if (num_bits != 16)
      data &= ((one << num_bits) - 1); // sets all MSBs past (num_bits-1) position to zero
    if (twos_comp && (navtools::GetBit<true>(data, num_bits-1))) {
      uint32_t signed_bit = (one << (num_bits-1));
      data &= (~signed_bit);
      data = (signed_bit - data); // un-signs the integer
      return -scale_factor * static_cast<T>(data);
    } else {
      return scale_factor * static_cast<T>(data);
    }
  }

  template<typename T>
  static T ParamFromBinary(uint32_t data, const T scale_factor,
    const uint8_t num_bits, const bool twos_comp)
  {
    assert((num_bits >= 1) && (num_bits <= 32));
    assert(scale_factor > T(0));

    static const uint32_t one = 1; // just avoiding the differing default size of integer literals

    if (num_bits != 32)
      data &= ((one << num_bits) - 1); // sets all MSBs past (num_bits-1) position to zero
    if (twos_comp && (navtools::GetBit<true>(data, num_bits-1))) {
      uint32_t signed_bit = (one << (num_bits-1));
      data &= (~signed_bit);
      data = (signed_bit - data); // un-signs the integer
      return -scale_factor * static_cast<T>(data);
    } else {
      return scale_factor * static_cast<T>(data);
    }
  }

  template<typename T>
  struct EphemInfo
  {
    T T_GD {0};
    T t_oc {0};
    T a_f2 {0};
    T a_f1 {0};
    T a_f0 {0};
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
  };

  template<typename Float>
  static constexpr EphemInfo<Float> scale_factors = 
  {
    .T_GD = navtools::PowerOfTwo<-31,Float>(),
    .t_oc = navtools::PowerOfTwo<4,Float>(),
    .a_f2 = navtools::PowerOfTwo<-55,Float>(),
    .a_f1 = navtools::PowerOfTwo<-43,Float>(),
    .a_f0 = navtools::PowerOfTwo<-31,Float>(),
    .M_0 = navtools::PowerOfTwo<-31,Float>(),
    .delta_n = navtools::PowerOfTwo<-43,Float>(),
    .e = navtools::PowerOfTwo<-33,Float>(),
    .sqrtA = navtools::PowerOfTwo<-19,Float>(),
    .OMEGA_0 = navtools::PowerOfTwo<-31,Float>(),
    .i_0 = navtools::PowerOfTwo<-31,Float>(),
    .omega = navtools::PowerOfTwo<-31,Float>(),
    .OMEGA_DOT = navtools::PowerOfTwo<-43,Float>(),
    .IDOT = navtools::PowerOfTwo<-43,Float>(),
    .C_uc = navtools::PowerOfTwo<-29,Float>(),
    .C_us = navtools::PowerOfTwo<-29,Float>(),
    .C_rc = navtools::PowerOfTwo<-5,Float>(),
    .C_rs = navtools::PowerOfTwo<-5,Float>(),
    .C_ic = navtools::PowerOfTwo<-29,Float>(),
    .C_is = navtools::PowerOfTwo<-29,Float>(),
    .t_oe = navtools::PowerOfTwo<4,Float>(),
  };

  static constexpr EphemInfo<uint8_t> num_bits =
  {
    .T_GD = 8,
    .t_oc = 16,
    .a_f2 = 8,
    .a_f1 = 16,
    .a_f0 = 22,
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
  };

  static constexpr EphemInfo<bool> signage =
  {
    .T_GD = true,
    .t_oc = false,
    .a_f2 = true,
    .a_f1 = true,
    .a_f0 = true,
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
    .t_oe = false
  };

  template<typename T>
  static int8_t T_GD(const T tgd)
  { return ParamToBinary<T>(tgd, scale_factors<T>.T_GD); }

  template<typename T>
  static uint16_t t_oc(const T toc)
  { return ParamToBinary<T>(toc, scale_factors<T>.t_oc); }

  template<typename T>
  static int8_t a_f2(const T af2)
  { return ParamToBinary<T>(af2, scale_factors<T>.a_f2); }
  
  template<typename T>
  static int16_t a_f1(const T af1)
  { return ParamToBinary<T>(af1, scale_factors<T>.a_f1); }
  
  template<typename T>
  static int32_t a_f0(const T af0)
  { return ParamToBinary<T>(af0, scale_factors<T>.a_f0); }

  template<typename T>
  static int32_t M_0(const T m0)
  { return ParamToBinary<T>(m0, scale_factors<T>.M_0); }

  template<typename T>
  static int16_t delta_n(const T deltan)
  { return ParamToBinary<T>(deltan, scale_factors<T>.delta_n); }

  template<typename T>
  static int32_t e(const T ec)
  { return ParamToBinary<T>(ec, scale_factors<T>.e); }

  template<typename T>
  static uint32_t sqrtA(const T sqrt)
  { return ParamToBinary<T>(sqrt, scale_factors<T>.sqrtA); }

  template<typename T>
  static int32_t OMEGA_0(const T omega0)
  { return ParamToBinary<T>(omega0, scale_factors<T>.OMEGA_0); }

  template<typename T>
  static int32_t i_0(const T i0)
  { return ParamToBinary<T>(i0, scale_factors<T>.i_0); }

  template<typename T>
  static int32_t omega(const T w)
  { return ParamToBinary<T>(w, scale_factors<T>.omega); }

  template<typename T>
  static int32_t OMEGA_DOT(const T omega_dot)
  { return ParamToBinary<T>(omega_dot, scale_factors<T>.OMEGA_DOT); }

  template<typename T>
  static int16_t IDOT(const T idot)
  { return ParamToBinary<T>(idot, scale_factors<T>.IDOT); }

  template<typename T>
  static int16_t C_uc(const T cuc)
  { return ParamToBinary<T>(cuc, scale_factors<T>.C_uc); }

  template<typename T>
  static int16_t C_us(const T cus)
  { return ParamToBinary<T>(cus, scale_factors<T>.C_us); }

  template<typename T>
  static int16_t C_rc(const T crc)
  { return ParamToBinary<T>(crc, scale_factors<T>.C_rc); }

  template<typename T>
  static int16_t C_rs(const T crs)
  { return ParamToBinary<T>(crs, scale_factors<T>.C_rs); }

  template<typename T>
  static int16_t C_ic(const T cic)
  { return ParamToBinary<T>(cic, scale_factors<T>.C_ic); }

  template<typename T>
  static int16_t C_is(const T cis)
  { return ParamToBinary<T>(cis, scale_factors<T>.C_is); }

  template<typename T>
  static uint16_t t_oe(const T toe)
  { return ParamToBinary<T>(toe, scale_factors<T>.t_oe); }

private:
  Word words_ [10]; // first bit chronologically is the MSB

public:
  // returns false if parity check failed
  //bool LoadParitySubframe(subframe);

  
  // -------------- WORDS 1 and 2 - APPLICABLE TO ALL SUBFRAMES --------------
  bool CheckPreamble() const
  { return (words_[0].GetSegment(0, 7) == 0x8B); }

  // 14 LSBs of result are the TLM message, with LSB transmitted last
  uint16_t GetTlmMessage() const
  { return words_[0].GetSegment(8,21); }
  
  bool GetIntegrityStatusFlag() const
  { return words_[0](22); }

  uint32_t GetTruncatedTow() const
  { return words_[1].GetSegment(0,16); }

  uint32_t GetTow() const
  { return GetTruncatedTow() << 2; }

  bool GetAlertFlag() const
  { return words_[1](17); }
  
  bool GetAntiSpoofFlag() const
  { return words_[1](18); }
  
  uint8_t GetSubframeId() const
  { return words_[1].GetSegment(19,21); }

  bool ValidSubframeId() const
  {
    uint8_t id = GetSubframeId();
    return (id > 0) && (id < 6);
  }


  // ---------------------- GETTING SUBFRAME 1 CONTENTS ----------------------
  template<typename T>
  T T_GD() const
  {
    return ParamFromBinary(words_[6].GetSegment(16,23), scale_factors<T>.T_GD,
                           num_bits.T_GD, signage.T_GD);
  }

  template<typename T>
  T t_oc() const
  {
    return ParamFromBinary(words_[7].GetSegment(8,23), scale_factors<T>.t_oc,
                           num_bits.t_oc, signage.t_oc);
  }

  template<typename T>
  T a_f2() const
  {
    return ParamFromBinary(words_[8].GetSegment(0,7), scale_factors<T>.a_f2,
                           num_bits.a_f2, signage.a_f2);
  }

  template<typename T>
  T a_f1() const
  {
    return ParamFromBinary(words_[8].GetSegment(8,23), scale_factors<T>.a_f1,
                           num_bits.a_f1, signage.a_f1);
  }

  template<typename T>
  T a_f0() const
  {
    return ParamFromBinary(words_[9].GetSegment(0,21), scale_factors<T>.a_f0,
                           num_bits.a_f0, signage.a_f0);
  }

  uint16_t IODC() const
  {
    uint16_t iode = words_[2].GetSegment(22,23) << 8;
    iode |= words_[7].GetSegment(0,7);
    return iode;
  }


  // ---------------------- GETTING SUBFRAME 2 CONTENTS ----------------------
  uint8_t IODE_sf2() const
  { return words_[2].GetSegment(0,7); }

  template<typename T>
  T C_rs() const
  {
    return ParamFromBinary(words_[2].GetSegment(8,23), scale_factors<T>.C_rs,
                           num_bits.C_rs, signage.C_rs);
  }

  template<typename T>
  T delta_n() const
  {
    return ParamFromBinary(words_[3].GetSegment(0,15), scale_factors<T>.delta_n,
                           num_bits.delta_n, signage.delta_n);
  }

  template<typename T>
  T M_0() const
  {
    uint32_t temp = words_[3].GetSegment(16,23) << 24;
    temp |= words_[4].GetSegment(0,23);
    return ParamFromBinary(temp, scale_factors<T>.M_0,
                           num_bits.M_0, signage.M_0);
  }

  template<typename T>
  T C_uc() const
  {
    return ParamFromBinary(words_[5].GetSegment(0,15), scale_factors<T>.C_uc,
                           num_bits.C_uc, signage.C_uc);
  }

  template<typename T>
  T e() const
  {
    uint32_t temp = words_[5].GetSegment(16,23) << 24;
    temp |= words_[6].GetSegment(0,23);
    return ParamFromBinary(temp, scale_factors<T>.e,
                           num_bits.e, signage.e);
  }

  template<typename T>
  T C_us() const
  {
    return ParamFromBinary(words_[7].GetSegment(0,15), scale_factors<T>.C_us,
                           num_bits.C_us, signage.C_us);
  }

  template<typename T>
  T sqrtA() const
  {
    uint32_t temp = words_[7].GetSegment(16,23) << 24;
    temp |= words_[8].GetSegment(0,23);
    return ParamFromBinary(temp, scale_factors<T>.sqrtA,
                           num_bits.sqrtA, signage.sqrtA);
  }

  template<typename T>
  T t_oe() const
  {
    return ParamFromBinary(words_[9].GetSegment(0,15), scale_factors<T>.t_oe,
                           num_bits.t_oe, signage.t_oe);
  }

  bool GetFitIntervalFlag() const
  { return words_[9](16); }

  uint8_t GetAodo() const
  { return words_[9].GetSegment(17,21); }


  // ---------------------- GETTING SUBFRAME 3 CONTENTS ----------------------
  template<typename T>
  T C_ic() const
  {
    return ParamFromBinary(words_[2].GetSegment(0,15), scale_factors<T>.C_ic,
                           num_bits.C_ic, signage.C_ic);
  }

  template<typename T>
  T OMEGA_0() const
  {
    uint32_t temp = words_[2].GetSegment(16,23) << 24;
    temp |= words_[3].GetSegment(0,23);
    return ParamFromBinary(temp, scale_factors<T>.OMEGA_0,
                           num_bits.OMEGA_0, signage.OMEGA_0);
  }

  template<typename T>
  T C_is() const
  {
    return ParamFromBinary(words_[4].GetSegment(0,15), scale_factors<T>.C_is,
                           num_bits.C_is, signage.C_is);
  }

  template<typename T>
  T i_0() const
  {
    uint32_t temp = words_[4].GetSegment(16,23) << 24;
    temp |= words_[5].GetSegment(0,23);
    return ParamFromBinary(temp, scale_factors<T>.i_0,
                           num_bits.i_0, signage.i_0);
  }

  template<typename T>
  T C_rc() const
  {
    return ParamFromBinary(words_[6].GetSegment(0,15), scale_factors<T>.C_rc,
                           num_bits.C_rc, signage.C_rc);
  }

  template<typename T>
  T omega() const
  {
    uint32_t temp = words_[6].GetSegment(16,23) << 24;
    temp |= words_[7].GetSegment(0,23);
    return ParamFromBinary(temp, scale_factors<T>.omega,
                           num_bits.omega, signage.omega);
  }

  template<typename T>
  T OMEGA_DOT() const
  {
    return ParamFromBinary(words_[8].GetSegment(0,23), scale_factors<T>.OMEGA_DOT,
                           num_bits.OMEGA_DOT, signage.OMEGA_DOT);
  }

  uint8_t IODE_sf3() const
  { return words_[9].GetSegment(0,7); }

  template<typename T>
  T IDOT() const
  {
    return ParamFromBinary(words_[9].GetSegment(8,21), scale_factors<T>.IDOT,
                           num_bits.IDOT, signage.IDOT);
  }


  // --------------------- SETTING WORD 1 and 2 CONTENTS ---------------------
  void SetBit(const int word_idx, const int bit_idx, const bool val)
  {
    assert((word_idx >= 0) && (word_idx < 10));
    assert((bit_idx >= 0) && (bit_idx < 30));
    words_[word_idx].SetBit(bit_idx, val);
  }

  void SetBit(int index, bool val)
  {
    assert((index >= 0) && (index < 300));
    this->SetBit(index / 30, index % 30, val);
  }

  void SetPreamble()
  { words_[0].SetSegment(0x8B,0,0,7); }

  void SetTlmMessage(const uint32_t tlm)
  { words_[0].SetSegment(tlm,8,0,13); }

  void SetIntegrityStatusFlag(const bool flag)
  { words_[0].SetBit(22,flag); }

  void SetTruncatedTow(const uint32_t tow_trunc)
  { words_[1].SetSegment(tow_trunc,0,2,18); }

  void SetTow(const uint32_t tow)
  { SetTruncatedTow(tow >> 2); }

  void SetAlertFlag(bool flag)
  { words_[1].SetBit(17,flag); }
  
  void SetAntiSpoofFlag(const bool flag)
  { words_[1].SetBit(18,flag); }

  void SetSubframeId(const uint8_t sf_id)
  { words_[1].SetSegment(sf_id,19,0,2); }


  // ------------------- SETTING SUBFRAME-SPECIFIC CONTENT -------------------
  // Sets all subframe 1 - specific info (words 3 - 10)
  template<typename T>
  void SetSubframe1Params(const KeplerElements<T>& ephems,
                          const uint32_t tow,
                          const uint16_t week_number,
                          const uint8_t l2_flag,
                          const uint8_t ura,
                          const uint8_t health)
  {
    SetPreamble();
    SetTow(tow);
    SetSubframeId(1);

    words_[2].SetSegment(0, week_number, 0, 9);
    words_[2].SetSegment(10, l2_flag, 0, 1);
    words_[2].SetSegment(12, ura, 0, 3);
    words_[2].SetSegment(16, health, 0, 5);
    words_[2].SetSegment(22, ephems.iodc, 8, 9);

    words_[6].SetSegment(16, T_GD(ephems.tgd), 0, 7);

    words_[7].SetSegment(0, ephems.iodc, 0, 7);
    words_[7].SetSegment(8, t_oc(ephems.toc), 0, 15);
    
    words_[8].SetSegment(0, a_f2(ephems.af2), 0, 7);
    words_[8].SetSegment(8, a_f1(ephems.af1), 0, 15);

    words_[9].SetSegment(0, a_f0(ephems.af0), 0, 21);
  }

  template<typename T>
  void SetSubframe2Params(const KeplerElements<T>& ephems,
                          const uint32_t tow,
                          const bool fit_interval_flag,
                          const uint8_t aodo)
  {
    SetPreamble();
    SetTow(tow);
    SetSubframeId(1);

    words_[2].SetSegment(0, (uint8_t) ephems.iode, 0, 7);
    words_[2].SetSegment(8, C_rs(ephems.crs), 0, 15);

    words_[3].SetSegment(0, delta_n(ephems.deltan), 0, 15);
    int32_t M = M_0(ephems.m0);
    words_[3].SetSegment(16, M, 24, 31);
    
    words_[4].SetSegment(0, M, 0, 23);

    words_[5].SetSegment(0, C_uc(ephems.cuc), 0, 15);
    int32_t e_bin = e(ephems.e);
    words_[5].SetSegment(16, e_bin, 24, 31);

    words_[6].SetSegment(0, e_bin, 0, 23);
  
    words_[7].SetSegment(0, C_us(ephems.cus), 0, 15);
    int32_t sqrtA_bin = sqrtA(ephems.sqrtA);
    words_[7].SetSegment(16, sqrtA_bin, 24, 31);

    words_[8].SetSegment(0, sqrtA_bin, 0, 23);

    words_[9].SetSegment(0, t_oe(ephems.toe), 0, 15);
    words_[9].SetBit(16, fit_interval_flag);
    words_[9].SetSegment(17, aodo, 0, 4);
  }

  template<typename T>
  void SetSubframe3Params(const KeplerElements<T>& ephems, const uint32_t tow)
  {
    SetPreamble();
    SetTow(tow);
    SetSubframeId(1);

    words_[2].SetSegment(0, C_ic(ephems.cic), 0, 15);
    int32_t omega0_bin = OMEGA_0(ephems.omega0);
    words_[2].SetSegment(16, omega0_bin, 24, 31);

    words_[3].SetSegment(0, omega0_bin, 0, 23);

    words_[4].SetSegment(0, C_is(ephems.cis), 0, 15);
    int32_t i0_bin = i_0(ephems.i0);
    words_[4].SetSegment(16, i0_bin, 24, 31);

    words_[5].SetSegment(0, i0_bin, 0, 23);

    words_[6].SetSegment(0, C_rc(ephems.crc), 0, 15);
    int32_t omega_bin = omega(ephems.omega);
    words_[6].SetSegment(16, omega_bin, 24, 31);
    
    words_[7].SetSegment(0, omega_bin, 0, 23);
    
    words_[8].SetSegment(0, OMEGA_DOT(ephems.omegaDot), 0, 23);

    words_[9].SetSegment(0, ephems.iode, 0, 7);
    words_[9].SetSegment(8, IDOT(ephems.iDot), 0, 13);
  }

};


#include <array>

// Prototype
class GpsSatellite
{
private:
  unsigned int id_;
  //TODO Ephemerides
  LnavData lnav_data [2];

  // ---------------------------- LNAV DATA STUFF ----------------------------
  // LNAV info - current values will be written when a new subframe is generated
  uint16_t tlm_message_ {0xAAAA};
  bool integrity_status_flag_ {false};
  bool alert_flag_ {false};
  bool anti_spoof_flag_ {false};
  uint8_t l2_flag_ {2};
  uint8_t ura_ {0}; // user range accuracy
  uint8_t health_ {0};
  // bool l2p_flag_ {false};
  bool fit_interval_flag_ {true};
  uint8_t aodo_ {0};

  uint8_t lnav_almanac_page_ {1};
  
  // latest last two values of subframe - for parity of next word
  bool d29_ {false};
  bool d30_ {false};

public:
  std::array<bool,1023> ca_code_;

  

};


} // namespace satutils
#endif
