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

// TODO: add models here

namespace satutils {

//* ===== Broadcast Ionosphere Corrections ===================================================== *//

template <typename Float>
class Ionosphere {
  public:
    virtual ~Ionosphere() = default;
    virtual void print() = 0;
};

template <typename Float>
struct Klobuchar : Ionosphere<Float> {
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
        std::cout << "--- Klobuchar Parameters ---" << std::endl
                  << "alpha_0:  " << alpha_0 << std::endl
                  << "alpha_1:  " << alpha_1 << std::endl
                  << "alpha_2:  " << alpha_2 << std::endl
                  << "alpha_3:  " << alpha_3 << std::endl
                  << "beta_0:   " << beta_0 << std::endl
                  << "beta_1:   " << beta_1 << std::endl
                  << "beta_2:   " << beta_2 << std::endl
                  << "beta_3:   " << beta_3 << std::endl
                  << "----------------------------" << std::endl;
    }
};

}  // namespace satutils

#endif