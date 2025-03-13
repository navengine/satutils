
#include <iostream>

#include "satutils/atmosphere.hpp"
#include "satutils/ephemeris.hpp"

int main() {
  satutils::KeplerEphem<double> eph;
  eph.toe = 0.0;

  double DoY = 123;
  double el = 2.8;
  double az = 3.0;
  double h = 215.0;
  double lat = 0.75;
  double lon = -1.5;
  double ToW = 412600;

  satutils::IonoModel<double> iono;
  iono.a0 = 2.6768e-08;
  iono.a1 = 4.4914e-09;
  iono.a2 = -3.2658e-07;
  iono.a3 = -5.2153e-07;
  iono.b0 = 1.3058e05;
  iono.b1 = -1.1203e05;
  iono.b2 = -7.0416e05;
  iono.b3 = -6.4865e06;
  double err1 = iono.CalcIonoDelay(ToW, lat, lon, az, el);
  std::cout << "err1 = " << err1 << std::endl;

  satutils::TropoModel<double> trop;
  double err2 = trop.CalcTropoDelay(DoY, lat, h, el);
  std::cout << "err2 = " << err2 << std::endl;

  return 0;
}