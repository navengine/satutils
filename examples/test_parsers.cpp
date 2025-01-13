
#include "satutils/rinex-parser.hpp"
#include "satutils/tle-parser.hpp"

int main() {
  // Initialize
  std::cout << "#####* TESTING TLE PARSER *#####" << std::endl;
  std::string filename = "celestrak_aug_23_2023_iridium-next.tle";
  std::cout << "File: " << filename << std::endl;

  // write TLE string to file
  std::string tle_string =
      "IRIDIUM 123             \n"
      "1 42804U 17039B   23234.63803740  .00000091  00000+0  25259-4 0  9999\n"
      "2 42804  86.3925  63.8997 0001981  98.2724 261.8696 14.34218517322388";
  std::ofstream fid(filename);
  if (!fid.is_open() || fid.bad()) {
    std::cout << "Failed to create TLE file! Ending 'test_parser' script!" << std::endl;
    return -1;
  }
  fid << tle_string;
  fid.close();

  // parse the tle file
  std::map<std::string, satutils::Sgp4Elements<double>> my_tle_map =
      satutils::TleParser<double>(filename);
  std::map<std::string, satutils::Sgp4Elements<double>>::iterator it = my_tle_map.begin();

  std::cout << "sv_id: " << it->first << "\ncatalog_id: " << it->second.catalog_id
            << "\nweek: " << it->second.week << "\ntoe: " << it->second.toe
            << "\nBstar: " << it->second.Bstar << "\ne0: " << it->second.e0
            << "\nomega0: " << it->second.omega0 << "\noemga: " << it->second.omega
            << "\ni0: " << it->second.i0 << "\nn0: " << it->second.n0
            << "\nnDot: " << it->second.nDot << "\nnDDot: " << it->second.nDDot
            << "\nm0: " << it->second.m0 << "\n\n";

  //! ==============================================================================================

  std::cout << "#####* TESTING RINEX PARSER *#####" << std::endl;
  filename = "fair_march_16_2023_gps_and_galileo.rnx";
  std::cout << "File: " << filename << std::endl;

  // write Rinex string to file
  std::string rnx_string =
      "     3.04           N: GNSS NAV DATA    E: MIXED            RINEX VERSION / TYPE\n"
      "sbf2rin-15.6.1                          20230317 000803 UTC PGM / RUN BY / DATE\n"
      "GPSA   2.6077E-08  7.4506E-09 -1.1921E-07  0.0000E+00       IONOSPHERIC CORR\n"
      "GPSB   1.2902E+05  0.0000E+00 -2.6214E+05  1.3107E+05       IONOSPHERIC CORR\n"
      "GAL    1.6700E+02 -3.1250E-02  7.5989E-03  0.0000E+00       IONOSPHERIC CORR\n"
      "    18                                                      LEAP SECONDS\n"
      "                                                            COMMENT\n"
      "                                                            COMMENT\n"
      "FAIR                                    MARKER NAME         COMMENT\n"
      "40408M001                               MARKER NUMBER       COMMENT\n"
      " -2281621.7717 -1453595.9493  5756961.9444                  COMMENT\n"
      "This data is provided as a public service by NASA/JPL.      COMMENT\n"
      "No warranty is expressed or implied regarding suitability   COMMENT\n"
      "for use.  For further information, contact:                 COMMENT\n"
      "ggnops at jpl dot nasa dot gov                              COMMENT\n"
      "                                                            END OF HEADER\n"
      "G07 2023 03 16 04 00 00 2.222126349807E-04-7.730704965070E-12 0.000000000000E+00\n"
      "     2.800000000000E+01 1.228125000000E+01 4.771270171125E-09-2.295500017430E+00\n"
      "     4.433095455170E-07 1.688047742937E-02 7.208436727524E-06 5.153755500793E+03\n"
      "     3.600000000000E+05 1.303851604462E-07 4.892990281257E-01-1.378357410431E-07\n"
      "     9.507183618301E-01 2.328437500000E+02-2.214618724925E+00-7.730679156674E-09\n"
      "     1.153619481453E-10 1.000000000000E+00 2.253000000000E+03 0.000000000000E+00\n"
      "     2.000000000000E+00 0.000000000000E+00-1.117587089539E-08 2.800000000000E+01\n"
      "     3.567180000000E+05 4.000000000000E+00\n"
      "E13 2023 03 16 04 00 00-1.774553675205E-05-3.552713678801E-13 0.000000000000E+00\n"
      "     8.800000000000E+01 1.440625000000E+01 2.720827619106E-09 8.080442274283E-01\n"
      "     5.997717380524E-07 2.912484342232E-04 1.465901732445E-06 5.440612743378E+03\n"
      "     3.600000000000E+05 3.725290298462E-08-2.713627018920E+00-6.891787052155E-08\n"
      "     9.994776247577E-01 3.293750000000E+02-4.282817581000E-01-5.644520831235E-09\n"
      "     1.096474243982E-10 5.170000000000E+02 2.253000000000E+03\n"
      "     3.120000000000E+00 0.000000000000E+00 6.053596735001E-09 6.286427378654E-09\n"
      "     3.606650000000E+05\n"
      "G30 2023 03 16 06 00 00-5.267057567835E-04 1.477928890381E-12 0.000000000000E+00\n"
      "     7.300000000000E+01 1.431250000000E+01 5.215217234731E-09-1.402269337324E+00\n"
      "     6.090849637985E-07 6.240330054425E-03 7.018446922302E-06 5.153616279602E+03\n"
      "     3.672000000000E+05 1.173466444016E-07 4.964383594732E-01-1.303851604462E-07\n"
      "     9.358414665820E-01 2.291250000000E+02-2.635181888783E+00-8.166411592393E-09\n"
      "     -4.964492505326E-11 1.000000000000E+00 2.253000000000E+03 0.000000000000E+00\n"
      "     2.000000000000E+00 0.000000000000E+00 3.725290298462E-09 7.300000000000E+01\n"
      "     3.610620000000E+05 4.000000000000E+00";
  fid.open(filename);
  if (!fid.is_open() || fid.bad()) {
    std::cout << "Failed to create Rinex file! Ending 'test_parser' script!" << std::endl;
    return -1;
  }
  fid << rnx_string;
  fid.close();

  // parse the rinex file
  auto my_rnx_map = satutils::RinexParser<double>(filename);
  for (auto &it : my_rnx_map) {
    std::cout << "sv_id: " << it.first << "\niode: " << it.second.second.iode
              << "\niodc: " << it.second.second.iodc << "\ntoe: " << it.second.second.toe
              << "\ntoc: " << it.second.second.toc << "\ntgd: " << it.second.second.tgd
              << "\noaf2: " << it.second.second.af2 << "\naf1: " << it.second.second.af1
              << "\naf0: " << it.second.second.af0 << "\ne: " << it.second.second.e
              << "\nsqrtA: " << it.second.second.sqrtA << "\ndeltan: " << it.second.second.deltan
              << "\nm0: " << it.second.second.m0 << "\nomega0: " << it.second.second.omega0
              << "\nomega: " << it.second.second.omega
              << "\nomegaDot: " << it.second.second.omegaDot << "\ni0: " << it.second.second.i0
              << "\niDot: " << it.second.second.iDot << "\ncuc: " << it.second.second.cuc
              << "\ncus: " << it.second.second.cus << "\ncic: " << it.second.second.cic
              << "\ncis: " << it.second.second.cis << "\ncrc: " << it.second.second.crc
              << "\ncrs: " << it.second.second.crs << "\nura: " << it.second.second.ura
              << "\nhealth: " << it.second.second.health << "\nalpha0: " << it.second.first.a0
              << "\nalpha1: " << it.second.first.a1 << "\nalpha2: " << it.second.first.a2
              << "\nalpha3: " << it.second.first.a3 << "\nbeta0: " << it.second.first.b0
              << "\nbeta1: " << it.second.first.b1 << "\nbeta2: " << it.second.first.b2
              << "\nbeta3: " << it.second.first.b3 << "\n\n";
  }
  return 0;
}