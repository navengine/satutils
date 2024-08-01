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

// TODO: add models here

namespace satutils {

//* ===== Broadcast Ionosphere Corrections ===================================================== *//

class Ionosphere {
  public:
    virtual ~Ionosphere() = default;
    virtual void print() = 0;
};

template <typename Float>
struct Klobuchar : Ionosphere {
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
};

}  // namespace satutils

#endif
