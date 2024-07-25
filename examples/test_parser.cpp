#include <cmath>
#include <filesystem>
#include <fstream>
#include <iostream>
#include <string>
#include <vector>

#include "navtools/constants.hpp"
#include "satutils/parsers.hpp"

// TODO: parsers no longer return satellite identifiers

const std::string RST = "\033[0m";  // reset
// const std::string RED = "\033[0;31m";   // red
// const std::string GRN = "\033[0;32m";   // green
// const std::string YEL = "\033[0;33m";   // yellow
// const std::string BLU = "\033[0;34m";   // blue
// const std::string MAG = "\033[0;35m";   // magenta
// const std::string CYN = "\033[0;36m";   // cyan
// const std::string WHT = "\033[0;37m";   // white
const std::string BRED = "\033[1;31m";  // bold red
const std::string BGRN = "\033[1;32m";  // bold green
const std::string BYEL = "\033[1;33m";  // bold yellow
// const std::string BBLU = "\033[1;34m";  // bold blue
// const std::string BMAG = "\033[1;35m";  // bold magenta
// const std::string BCYN = "\033[1;36m";  // bold cyan
// const std::string BWHT = "\033[1;37m";  // bold white

template <typename T>
bool test(T a, T b) {
    if (a == b) {
        return true;
    }
    return false;
}

int main() {
    // Initialize
    std::cout << std::setprecision(20) << std::endl;
    std::cout << BYEL << "#####* TESTING RINEX PARSER *#####" << std::endl << RST;
    std::filesystem::path cwd = std::filesystem::current_path();
    std::string filename = "fair_march_16_2023_gps_and_galileo.rnx";
    std::string full_file = (cwd / filename).string();
    std::cout << "Path: " << cwd.string() << std::endl;
    std::cout << "File: " << filename << std::endl;

    // Write Rinex string to the file
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
    std::ofstream fid(full_file);
    if (!fid.is_open() || fid.bad()) {
        std::cout << "Failed to create rinex file! Ending 'test_parser' script!" << std::endl;
        return -1;
    }
    fid << rnx_string;
    fid.close();

    // Read Rinex file
    std::istringstream b0("2023 03 16 04 00 30");
    std::chrono::system_clock::time_point t0;
    std::chrono::from_stream(b0, "%Y %m %d %H %M %S", t0);
    // std::tm tm0;
    // b0 >> std::get_time(&tm0, "%Y %m %d %H %M %S");
    // t0 = std::chrono::system_clock::from_time_t(std::mktime(&tm0));

    satutils::RinexParser<double> rp(full_file);
    std::vector<satutils::KeplerEphemeris<double>> gps_eph_vec, galileo_eph_vec;
    std::vector<satutils::PolyClock<double>> gps_clk_vec, galileo_clk_vec;
    std::vector<satutils::Klobuchar<double>> gps_atm_vec, galileo_atm_vec;
    rp.getGps(t0, gps_eph_vec, gps_clk_vec, gps_atm_vec);
    rp.getGalileo(t0, galileo_eph_vec, galileo_clk_vec, galileo_atm_vec);

    std::vector<bool> results;
    // results.push_back(test(gps_eph_vec[0].id, (std::string) "G07"));
    // results.push_back(test(gps_eph_vec[0].leap_seconds, 18));
    // results.push_back(test(gps_eph_vec[0].l2c_flag, true));
    // results.push_back(test(gps_eph_vec[0].l2p_flag, false));
    // results.push_back(test(gps_eph_vec[0].accuracy, 2.000000000000E+00));
    // results.push_back(test(gps_eph_vec[0].health, (uint8_t)0));
    results.push_back(test(gps_eph_vec[0].week, 2253));
    results.push_back(test(gps_eph_vec[0].IODE, 2.800000000000E+01));
    results.push_back(test(gps_eph_vec[0].C_rs, 1.228125000000E+01));
    results.push_back(test(gps_eph_vec[0].delta_n, 4.771270171125E-09));
    results.push_back(test(gps_eph_vec[0].M_0, -2.295500017430E+00));
    results.push_back(test(gps_eph_vec[0].C_uc, 4.433095455170E-07));
    results.push_back(test(gps_eph_vec[0].e, 1.688047742937E-02));
    results.push_back(test(gps_eph_vec[0].C_us, 7.208436727524E-06));
    results.push_back(test(gps_eph_vec[0].A, std::pow(5.153755500793E+03, 2)));
    results.push_back(test(gps_eph_vec[0].t_oe, 3.600000000000E+05));
    results.push_back(test(gps_eph_vec[0].C_ic, 1.303851604462E-07));
    results.push_back(test(gps_eph_vec[0].OMEGA_0, 4.892990281257E-01));
    results.push_back(test(gps_eph_vec[0].C_is, -1.378357410431E-07));
    results.push_back(test(gps_eph_vec[0].i_0, 9.507183618301E-01));
    results.push_back(test(gps_eph_vec[0].C_rc, 2.328437500000E+02));
    results.push_back(test(gps_eph_vec[0].omega, -2.214618724925E+00));
    results.push_back(test(gps_eph_vec[0].OMEGA_DOT, -7.730679156674E-09));
    results.push_back(test(gps_eph_vec[0].IDOT, 1.153619481453E-10));
    results.push_back(test(gps_clk_vec[0].a_f0, 2.222126349807E-04));
    results.push_back(test(gps_clk_vec[0].a_f1, -7.730704965070E-12));
    results.push_back(test(gps_clk_vec[0].a_f2, 0.000000000000E+00));
    results.push_back(test(gps_clk_vec[0].T_GD, -1.117587089539E-08));
    results.push_back(test(gps_clk_vec[0].IODC, 2.800000000000E+01));
    // TODO: time of clock
    results.push_back(test(gps_atm_vec[0].alpha_0, 2.6077E-08));
    results.push_back(test(gps_atm_vec[0].alpha_1, 7.4506E-09));
    results.push_back(test(gps_atm_vec[0].alpha_2, -1.1921E-07));
    results.push_back(test(gps_atm_vec[0].alpha_3, 0.0000E+00));
    results.push_back(test(gps_atm_vec[0].beta_0, 1.2902E+05));
    results.push_back(test(gps_atm_vec[0].beta_1, 0.0000E+00));
    results.push_back(test(gps_atm_vec[0].beta_2, -2.6214E+05));
    results.push_back(test(gps_atm_vec[0].beta_3, 1.3107E+05));

    // results.push_back(test(galileo_eph_vec[0].id, (std::string) "E13"));
    // results.push_back(test(galileo_eph_vec[0].leap_seconds, 18));
    // results.push_back(test(galileo_eph_vec[0].data_source_flag, (uint16_t)517));
    // results.push_back(test(galileo_eph_vec[0].week, 2253));
    // results.push_back(test(galileo_eph_vec[0].accuracy, 3.120000000000E+00));
    // results.push_back(test(galileo_eph_vec[0].health, (uint8_t)0));
    results.push_back(test(galileo_eph_vec[0].IODE, 8.800000000000E+01));
    results.push_back(test(galileo_eph_vec[0].C_rs, 1.440625000000E+01));
    results.push_back(test(galileo_eph_vec[0].delta_n, 2.720827619106E-09));
    results.push_back(test(galileo_eph_vec[0].M_0, 8.080442274283E-01));
    results.push_back(test(galileo_eph_vec[0].C_uc, 5.997717380524E-07));
    results.push_back(test(galileo_eph_vec[0].e, 2.912484342232E-04));
    results.push_back(test(galileo_eph_vec[0].C_us, 1.465901732445E-06));
    results.push_back(test(galileo_eph_vec[0].A, std::pow(5.440612743378E+03, 2)));
    results.push_back(test(galileo_eph_vec[0].t_oe, 3.600000000000E+05));
    results.push_back(test(galileo_eph_vec[0].C_ic, 3.725290298462E-08));
    results.push_back(test(galileo_eph_vec[0].OMEGA_0, -2.713627018920E+00));
    results.push_back(test(galileo_eph_vec[0].C_is, -6.891787052155E-08));
    results.push_back(test(galileo_eph_vec[0].i_0, 9.994776247577E-01));
    results.push_back(test(galileo_eph_vec[0].C_rc, 3.293750000000E+02));
    results.push_back(test(galileo_eph_vec[0].omega, -4.282817581000E-01));
    results.push_back(test(galileo_eph_vec[0].OMEGA_DOT, -5.644520831235E-09));
    results.push_back(test(galileo_eph_vec[0].IDOT, 1.096474243982E-10));
    results.push_back(test(galileo_clk_vec[0].a_f0, -1.774553675205E-05));
    results.push_back(test(galileo_clk_vec[0].a_f1, -3.552713678801E-13));
    results.push_back(test(galileo_clk_vec[0].a_f2, 0.000000000000E+00));
    results.push_back(test(galileo_clk_vec[0].T_GD, 6.053596735001E-09));
    // results.push_back(test(galileo_clk_vec[0].BGD_e5a_e1, 6.053596735001E-09));
    // results.push_back(test(galileo_clk_vec[0].BGD_e5b_e1, 6.286427378654E-09));
    // TODO: time of clock
    results.push_back(test(galileo_atm_vec[0].alpha_0, 1.6700E+02));
    results.push_back(test(galileo_atm_vec[0].alpha_1, -3.1250E-02));
    results.push_back(test(galileo_atm_vec[0].alpha_2, 7.5989E-03));
    results.push_back(test(galileo_atm_vec[0].alpha_3, 0.0000E+00));

    // check rinex results
    size_t total_passed = 0;
    for (std::vector<bool>::iterator it = results.begin(); it != results.end(); ++it) {
        total_passed += *it;
    }
    if (total_passed == results.size()) {
        std::cout << BGRN << total_passed << "/" << results.size() << " Rinex parser tests passed!"
                  << RST << std::endl;
    } else {
        std::cout << BRED << "Rinex parser tests failed! " << total_passed << "/" << results.size()
                  << " tests succeeded." << RST << std::endl;
    }

    // Delete rinex file
    std::cout << std::endl;
    std::remove(full_file.c_str());

    // Initialize
    std::cout << BYEL << "#####* TESTING TLE PARSER *#####" << std::endl << RST;
    cwd = std::filesystem::current_path();
    filename = "celestrak_aug_23_2023_iridium-next.tle";
    full_file = (cwd / filename).string();
    std::cout << "Path: " << cwd.string() << std::endl;
    std::cout << "File: " << filename << std::endl;

    // write TLE string to file
    std::string tle_string =
            "IRIDIUM 123             \n"
            "1 42804U 17039B   23234.63803740  .00000091  00000+0  25259-4 0  9999\n"
            "2 42804  86.3925  63.8997 0001981  98.2724 261.8696 14.34218517322388";
    fid.open(full_file);
    if (!fid.is_open() || fid.bad()) {
        std::cout << "Failed to create TLE file! Ending 'test_parser' script!" << std::endl;
        return -1;
    }
    fid << tle_string;
    fid.close();
    results.clear();

    // Read the TLE file
    satutils::TLEParser<double> tp(full_file);
    std::vector<satutils::SGP4Ephemeris<double>> sgp4_eph_vec;
    tp.get(sgp4_eph_vec);

    // results.push_back(test(sgp4_eph_vec[0].id, (std::string) "IRIDIUM 123"));
    // results.push_back(test(sgp4_eph_vec[0].N, navtools::TWO_PI<double> * 32238.0));
    // results.push_back(test(sgp4_eph_vec[0].set_number, 999));
    results.push_back(test(sgp4_eph_vec[0].catalog_id, 42804));
    results.push_back(test(sgp4_eph_vec[0].week, 2276));
    results.push_back(test(std::round(sgp4_eph_vec[0].t_oe * 1e5) / 1e5, 227926.43136));
    results.push_back(
            test(sgp4_eph_vec[0].n_dot, navtools::TWO_PI<double> / 2073600.0 * 0.00000091));
    results.push_back(test(sgp4_eph_vec[0].n_ddot, navtools::TWO_PI<double> / 2985984000.0 * 0.0));
    results.push_back(test(sgp4_eph_vec[0].B_star, 0.25259e-4));
    results.push_back(test(sgp4_eph_vec[0].i_0, navtools::DEG2RAD<double> * 86.3925));
    results.push_back(test(sgp4_eph_vec[0].OMEGA_0, navtools::DEG2RAD<double> * 63.8997));
    results.push_back(test(sgp4_eph_vec[0].e_0, 0.0001981));
    results.push_back(test(sgp4_eph_vec[0].omega, navtools::DEG2RAD<double> * 98.2724));
    results.push_back(test(sgp4_eph_vec[0].M_0, navtools::DEG2RAD<double> * 261.8696));
    results.push_back(test(sgp4_eph_vec[0].n_0, navtools::TWO_PI<double> / 1440.0 * 14.34218517));

    // check rinex results
    total_passed = 0;
    for (std::vector<bool>::iterator it = results.begin(); it != results.end(); ++it) {
        total_passed += *it;
    }
    if (total_passed == results.size()) {
        std::cout << BGRN << total_passed << "/" << results.size() << " TLE parser tests passed!"
                  << RST << std::endl;
    } else {
        std::cout << BRED << "TLE parser tests failed! " << total_passed << "/" << results.size()
                  << RST << " tests succeeded." << std::endl;
    }

    // Delete TLE file
    std::cout << std::endl;
    std::remove(full_file.c_str());

    return 0;
}