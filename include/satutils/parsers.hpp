/**
|========================================= parsers.hpp ============================================|
|                                                                                                  |
|   @file     include/satutils/parsers.hpp                                                         |
|   @brief    Tools to parse navigation RINEX files and TLE files into ephemeris parameters.       |
|   @date     July 2024                                                                            |
|                                                                                                  |
|==================================================================================================|
*/

// TODO: add function calls to only return requested satellites

#ifndef SATUTILS_PARSERS_HPP
#define SATUTILS_PARSERS_HPP

// #include <algorithm>
// #include <cmath>
// #include <iomanip>
// #include <sstream>

#include <chrono>
#include <fstream>
#include <iostream>
#include <iterator>
#include <map>
#include <memory>
#include <string>
#include <vector>

#include <satutils/common.hpp>
#include "navtools/constants.hpp"
#include "satutils/atmosphere.hpp"
#include "satutils/clock.hpp"
#include "satutils/ephemeris.hpp"

namespace satutils {

// enum ConstellationType = {GPS = 0, Galileo = 1, Glonass = 2};

//* ===== Rinex Tools ========================================================================== *//

template <typename Float>
class RinexParser {
  public:
    //! === RINEXPARSER ===
    RinexParser(std::string filename) {
        // safely open file
        fid_ = std::ifstream(filename);
        if (!fid_.is_open()) {
            std::cout << "navsat::RinexParser::RinexParser - Invalid file!" << std::endl;
            // return std::vector<std::unique_ptr<Ephemeris>>;
        } else if (fid_.bad()) {
            std::cout << "navsat::RinexParser::RinexParser - Failed to read file!" << std::endl;
            // return std::vector<std::unique_ptr<Ephemeris>>;
        }
        fid_.seekg(0);

        init();
    }

    //! === ~RINEXPARSER ===
    ~RinexParser() {
        fid_.close();
    }

    //! === GETGPS ===
    void getGps(
            time_point &tow,
            std::vector<KeplerEphemeris<Float>> &eph,
            std::vector<PolyClock<Float>> &clk,
            std::vector<Klobuchar<Float>> &atm) {
        // convert time ranges to GpsTime
        GpsTime t0;
        date2gpsTime(t0, tow);

        for (std::map<std::string, std::vector<std::pair<int, GpsTime>>>::iterator m =
                     gps_idx_.begin();
             m != gps_idx_.end();
             ++m) {
            std::cout << m->first << " ";

            KeplerEphemeris<Float> tmp_eph;
            PolyClock<Float> tmp_clk;
            Klobuchar<Float> tmp_atm;
            FindCorrectGps(m->second, t0, tmp_eph, tmp_clk, tmp_atm);
            eph.push_back(tmp_eph);
            clk.push_back(tmp_clk);
            atm.push_back(tmp_atm);

            std::cout << std::endl;
        }
    }

    //! === GETGALILEO ===
    void getGalileo(
            time_point &tow,
            std::vector<KeplerEphemeris<Float>> &eph,
            std::vector<PolyClock<Float>> &clk,
            std::vector<Klobuchar<Float>> &atm) {
        // convert time ranges to GpsTime
        GpsTime t0;
        date2gpsTime(t0, tow);

        for (std::map<std::string, std::vector<std::pair<int, GpsTime>>>::iterator m =
                     galileo_idx_.begin();
             m != galileo_idx_.end();
             ++m) {
            std::cout << m->first << " ";

            KeplerEphemeris<Float> tmp_eph;
            PolyClock<Float> tmp_clk;
            Klobuchar<Float> tmp_atm;
            FindCorrectGalileo(m->second, t0, tmp_eph, tmp_clk, tmp_atm);
            eph.push_back(tmp_eph);
            clk.push_back(tmp_clk);
            atm.push_back(tmp_atm);

            std::cout << std::endl;
        }
    }

  private:
    std::ifstream fid_;
    Float leap_seconds_;
    std::vector<Float> gps_iono_alpha_;
    std::vector<Float> gps_iono_beta_;
    std::vector<Float> galileo_iono_;
    std::map<std::string, std::vector<std::pair<int, GpsTime>>> gps_idx_;
    std::map<std::string, std::vector<std::pair<int, GpsTime>>> galileo_idx_;

    //! === INIT ===
    void init() {
        ParseHeader();

        // scan file to find file locations where satellites occur
        std::string line, sv_id;
        int idx = 0;
        GpsTime gt;
        while (!fid_.eof()) {
            std::getline(fid_, line, '\n');
            // --- GPS BLOCK ---
            if (line[0] == 'G') {
                idx = static_cast<int>(fid_.tellg()) - line.length() - 1;
                gt = ParseTimeBlock(line.substr(4, 20));
                sv_id = line.substr(0, 3);

                if (gps_idx_.find(sv_id) == gps_idx_.end()) {
                    // new prn
                    gps_idx_.insert(std::make_pair(sv_id, std::vector<std::pair<int, GpsTime>>(0)));
                }

                gps_idx_[sv_id].push_back(std::make_pair(idx, gt));
            }

            // --- GALILEO BLOCK ---
            else if (line[0] == 'E') {
                idx = static_cast<int>(fid_.tellg()) - line.length() - 1;
                gt = ParseTimeBlock(line.substr(4, 20));
                sv_id = line.substr(0, 3);

                if (galileo_idx_.find(sv_id) == galileo_idx_.end()) {
                    // new prn
                    galileo_idx_.insert(
                            std::make_pair(sv_id, std::vector<std::pair<int, GpsTime>>(0)));
                }

                galileo_idx_[sv_id].push_back(std::make_pair(idx, gt));
            }
        }
        fid_.clear(); // needed because end of file flag was likely set
    }

    //! === PARSEHEADER ===
    void ParseHeader() {
        const std::string gps_alpha_token = "GPSA";
        const std::string gps_beta_token = "GPSB";
        const std::string galileo_iono_token = "GAL";
        const std::string iono_token = "IONOSPHERIC CORR";
        const std::string time_token = "LEAP SECONDS";
        const std::string eoh_token = "END OF HEADER";
        const std::string comment_token = "COMMENT";

        // parse header
        std::string line;  // current line from file
        leap_seconds_ = 0;
        while (!fid_.eof()) {
            std::getline(fid_, line, '\n');

            // ignore comment lines
            if (line.find(comment_token) != std::string::npos) {
                continue;
            }

            // check for ionosphere parameters
            else if (line.find(iono_token) != std::string::npos) {
                std::istringstream buf(line.substr(5, 55));
                if (line.find(gps_alpha_token) != std::string::npos) {
                    // found gps alpha parameters
                    gps_iono_alpha_ = {
                            std::istream_iterator<Float>(buf), std::istream_iterator<Float>()};
                } else if (line.find(gps_beta_token) != std::string::npos) {
                    // found gps beta parameters
                    gps_iono_beta_ = {
                            std::istream_iterator<Float>{buf}, std::istream_iterator<Float>()};
                } else if (line.find(galileo_iono_token) != std::string::npos) {
                    // found gps beta parameters
                    galileo_iono_ = {
                            std::istream_iterator<Float>{buf}, std::istream_iterator<Float>()};
                }
            }

            // check for UTC leap second correction
            else if (line.find(time_token) != std::string::npos) {
                leap_seconds_ = std::stoi(line.substr(4, 6));
            }

            // check for end of header
            else if (line.find(eoh_token) != std::string::npos) {
                break;
            }
        }
    }

    //! === PARSETIMEBLOCK ===
    GpsTime ParseTimeBlock(const std::string &timeline) {
        // std::cout << "navsat::RinexParser::ParseTimeBlock - " << timeline << std::endl;
        std::istringstream buf(timeline);

        //! C++20 Function! - Requires GCC14 or Clang/libc++17
        time_point dt;
        std::chrono::from_stream(buf, "%Y %m %d %H %M %S", dt);
        // std::tm t;
        // buf >> std::get_time(&t, "%Y %m %d %H %M %S");
        // time_point dt = std::chrono::system_clock::from_time_t(std::mktime(&t));
        GpsTime gpsT;
        date2gpsTime(gpsT, dt);

        return gpsT;
    }

    //! === PARSENAVBLOCK ===
    std::vector<Float> ParseNavBlock(const std::string &navblock) {
        std::vector<Float> navdata;
        std::string navword;

        // split word every 19 characters
        for (size_t i = 0; i < navblock.length(); i += 19) {
            navword = navblock.substr(i, 19);
            if (navword.find_first_not_of(' ') == std::string::npos) {
                continue;
            }
            if (navword.find('D') != std::string::npos) {
                std::replace(navword.begin(), navword.end(), 'D', 'E');
            }
            // std::cout << "navsat::navDataSplitter - word = " << navword << std::endl;
            // std::cout << "navsat::navDataSplitter - result = " << std::stod(navword) <<
            // std::endl;
            navdata.push_back(std::stod(navword));
        }
        return navdata;
    }

    //! === MINPOSITIVEIDX ===
    size_t MinPositiveIdx(
            const std::vector<std::pair<int, GpsTime>> &arr, const Float &set_point = 0.0) {
        size_t result = 0;
        Float tmp;
        Float max_val = std::numeric_limits<Float>::max();
        for (size_t i = 0; i < arr.size(); i++) {
            tmp = arr[i].second.second - set_point;
            if ((0 < tmp) && (tmp < max_val)) {
                max_val = tmp;
                result = i;
            }
        }
        return result;
    }

    //! === FINDCORRECTGPS ===
    void FindCorrectGps(
            const std::vector<std::pair<int, GpsTime>> &v,
            GpsTime &t0,
            KeplerEphemeris<Float> &eph,
            PolyClock<Float> &clk,
            Klobuchar<Float> &atm) {
        // find correct index
        size_t i = MinPositiveIdx(v, t0.second);

        // parse first line
        std::string line, navline;
        fid_.seekg(v[i].first);
        std::getline(fid_, line, '\n');

        // parse entire block
        std::string sv_id = line.substr(0, 3);
        navline = line.substr(23, line.length());
        for (int i = 1; i < 8; i++) {
            std::getline(fid_, line, '\n');
            navline += line.substr(4, line.length());
        }
        std::vector<Float> param = ParseNavBlock(navline);
        // std::cout << navline << std::endl << std::endl;

        // eph.id = sv_id;
        // eph.leap_seconds = leap_seconds_;
        // eph.l2c_flag = static_cast<bool>(param[20]);   //
        // eph.l2p_flag = static_cast<bool>(param[22]);   //
        // eph.accuracy = param[23];                      // [m]
        // eph.health = static_cast<uint8_t>(param[24]);  //

        clk.t_oc = v[i].second.second;  // [s]
        clk.a_f0 = param[0];            // [s]
        clk.a_f1 = param[1];            // [s/s]
        clk.a_f2 = param[2];            // [s/s^2]
        clk.IODC = param[26];           // [s]
        clk.T_GD = param[25];           // [s]

        eph.IODE = param[3];            // [s]
        eph.C_rs = param[4];            // [m]
        eph.delta_n = param[5];         // [rad/s]
        eph.M_0 = param[6];             // [rad]
        eph.C_uc = param[7];            // [rad]
        eph.e = param[8];               //
        eph.C_us = param[9];            // [rad]
        eph.A = param[10] * param[10];  // [(m]
        eph.t_oe = param[11];           // [s]
        eph.C_ic = param[12];           // [rad]
        eph.OMEGA_0 = param[13];        // [rad]
        eph.C_is = param[14];           // [rad]
        eph.i_0 = param[15];            // [rad]
        eph.C_rc = param[16];           // [m]
        eph.omega = param[17];          // [rad]
        eph.OMEGA_DOT = param[18];      // [rad/s]
        eph.IDOT = param[19];           // [rad/s]
        eph.week = param[21];           //

        atm.alpha_0 = gps_iono_alpha_[0];  //
        atm.alpha_1 = gps_iono_alpha_[1];  //
        atm.alpha_2 = gps_iono_alpha_[2];  //
        atm.alpha_3 = gps_iono_alpha_[3];  //
        atm.beta_0 = gps_iono_beta_[0];    //
        atm.beta_1 = gps_iono_beta_[1];    //
        atm.beta_2 = gps_iono_beta_[2];    //
        atm.beta_3 = gps_iono_beta_[3];    //
    }

    //! === FINDCORRECTGALILEO ===
    void FindCorrectGalileo(
            const std::vector<std::pair<int, GpsTime>> &v,
            GpsTime &t0,
            KeplerEphemeris<Float> &eph,
            PolyClock<Float> &clk,
            Klobuchar<Float> &atm) {
        // find correct index
        size_t i = MinPositiveIdx(v, t0.second);

        // parse first line
        std::string line, navline;
        fid_.seekg(v[i].first);
        std::getline(fid_, line, '\n');

        // parse entire block
        std::string sv_id = line.substr(0, 3);
        navline = line.substr(23, line.length());
        for (int i = 1; i < 8; i++) {
            std::getline(fid_, line, '\n');
            navline += line.substr(4, line.length());
        }
        std::vector<Float> param = ParseNavBlock(navline);
        // std::cout << navline << std::endl << std::endl;

        // eph.id = sv_id;
        // eph.leap_seconds = h_.leap_seconds;
        // eph.accuracy = param[22];                                 // [m]
        // eph.health = static_cast<uint8_t>(param[23]);             //
        // eph.data_source_flag = static_cast<uint16_t>(param[20]);  //

        clk.t_oc = v[i].second.second;  // [s]
        clk.a_f0 = param[0];            // [s]
        clk.a_f1 = param[1];            // [s/s]
        clk.a_f2 = param[2];            // [s/s^2]
        clk.T_GD = param[24];           // [s]
        // clk.BGD_e5a_e1 = param[24];   // [s]
        // clk.BGD_e5b_e1 = param[25];   // [s]

        eph.IODE = param[3];            // [s]
        eph.C_rs = param[4];            // [m]
        eph.delta_n = param[5];         // [rad/s]
        eph.M_0 = param[6];             // [rad]
        eph.C_uc = param[7];            // [rad]
        eph.e = param[8];               //
        eph.C_us = param[9];            // [rad]
        eph.A = param[10] * param[10];  // [m]
        eph.t_oe = param[11];           // [s]
        eph.C_ic = param[12];           // [rad]
        eph.OMEGA_0 = param[13];        // [rad]
        eph.C_is = param[14];           // [rad]
        eph.i_0 = param[15];            // [rad]
        eph.C_rc = param[16];           // [m]
        eph.omega = param[17];          // [rad]
        eph.OMEGA_DOT = param[18];      // [rad/s]
        eph.IDOT = param[19];           // [rad/s]
        eph.week = param[21];           //

        atm.alpha_0 = galileo_iono_[0];  //
        atm.alpha_1 = galileo_iono_[1];  //
        atm.alpha_2 = galileo_iono_[2];  //
        atm.alpha_3 = galileo_iono_[3];  //
    }
};

//* ===== TLE Tools ============================================================================ *//

template <typename Float>
class TLEParser {
    // https://celestrak.org/publications/AIAA/2008-6770/AIAA-2008-6770.pdf
  public:
    //! === TLEPARSER ===
    TLEParser(std::string filename) {
        // safely open file
        fid_ = std::ifstream(filename);
        if (!fid_.is_open()) {
            std::cout << "navsat::TLEParser::TLEParser - Invalid file!" << std::endl;
        } else if (fid_.bad()) {
            std::cout << "navsat::TLEParser::TLEParser - Failed to read file!" << std::endl;
        }
        fid_.seekg(0);
    }

    //! === ~TLEPARSER ===
    ~TLEParser() {
        fid_.close();
    }

    //! === GET ===
    void get(std::vector<SGP4Ephemeris<Float>> &eph) {
        std::string line, sv_id;
        time_point dt;
        GpsTime gpsT;
        int catalog_id;
        Float n_dot, n_ddot, B_star, i_0, OMEGA_0, e, omega, M_0, n_0;

        while (!fid_.eof()) {
            // every 3 lines is a satellite
            for (int i = 0; i < 3; i++) {
                std::getline(fid_, line, '\n');
                switch (i) {
                    case (0): {
                        // satellite id
                        // std::cout << line << std::endl;
                        std::istringstream buf(line);
                        std::vector<std::string> tmp{
                                std::istream_iterator<std::string>(buf),
                                std::istream_iterator<std::string>()};
                        sv_id = tmp[0] + " " + tmp[1];
                        std::cout << sv_id << std::endl;
                        break;
                    }
                    case (1): {
                        // line 1
                        // std::cout << line << std::endl;
                        catalog_id = std::stoi(line.substr(2, 5));

                        std::istringstream buf(line.substr(18, 5));
                        std::chrono::from_stream(buf, "%y%j", dt);
                        uint64_t ns = static_cast<uint64_t>(
                                1e9 * 86400.0 * std::stod(line.substr(23, 9)));
                        // std::cout << "ns = " << ns << ", " << line.substr(23, 9) << std::endl;
                        dt += std::chrono::nanoseconds{ns};
                        date2gpsTime(gpsT, dt);

                        n_dot = std::stod(line.substr(33, 10));
                        n_ddot = assumedDecimalPoint(line.substr(44, 8));
                        B_star = assumedDecimalPoint(line.substr(53, 8));
                        // set_number = std::stoi(line.substr(64, 4));
                        // std::cout << "checksum line 1 = " << checksum(line.substr(0,
                        // line.length() - 1))
                        //           << std::endl;
                        if (checksum(line.substr(0, line.length() - 1)) !=
                            std::stoi(line.substr(68, 1))) {
                            std::cout << "navsat::TLEParser::getSatellites - ERROR " << sv_id
                                      << " line1 checksum incorrect" << std::endl;
                        }
                        break;
                    }
                    case (2): {
                        // line 2
                        // std::cout << line << std::endl;
                        i_0 = std::stod(line.substr(8, 8));
                        OMEGA_0 = std::stod(line.substr(17, 8));
                        e = assumedDecimalPoint(line.substr(26, 7));
                        omega = std::stod(line.substr(34, 8));
                        M_0 = std::stod(line.substr(43, 8));
                        n_0 = std::stod(line.substr(52, 11));
                        // N = std::stod(line.substr(63, 5));
                        // std::cout << "checksum line 2 = " << checksum(line.substr(0,
                        // line.length() - 1))
                        //           << std::endl;
                        if (checksum(line.substr(0, line.length() - 1)) !=
                            std::stoi(line.substr(68, 1))) {
                            std::cout << "navsat::TLEParser::getSatellites - ERROR " << sv_id
                                      << " line2 checksum incorrect" << std::endl;
                        }
                        break;
                    }
                }
            }
            // create satellite
            SGP4Ephemeris<Float> tmp;
            // tmp.id = sv_id;
            // tmp.set_number = set_number;
            // tmp.N = N;
            tmp.catalog_id = catalog_id;
            tmp.week = gpsT.week;
            tmp.t_oe = gpsT.second;
            tmp.n_dot = TWO_PI<Float> / 2073600.0 * n_dot;  // (XNDT2O) [rev/day^2] -> [rad/min^2]
            tmp.n_ddot =
                    TWO_PI<Float> / 2985984000.0 * n_ddot;  // (XNDD6O) [rev/day^3] -> [rad/min^3]
            tmp.B_star = B_star;                            // (BSTAR)
            tmp.i_0 = DEG2RAD<Float> * i_0;                 // (XINCL)  [deg] -> [rad]
            tmp.OMEGA_0 = DEG2RAD<Float> * OMEGA_0;         // (OMEGA0) [deg] -> [rad]
            tmp.e_0 = e;                                    // (E0)
            tmp.omega = DEG2RAD<Float> * omega;             // (XNODEO) [deg] -> [rad]
            tmp.M_0 = DEG2RAD<Float> * M_0;                 // (XMO)    [deg] -> [rad]
            tmp.n_0 = TWO_PI<Float> / 1440.0 * n_0;         // (XN0)    [rev/day] -> [rad/min]
            eph.push_back(tmp);
        }
    }

  private:
    std::ifstream fid_;

    //! === ASSUMEDDECIMALPOINT ===
    Float assumedDecimalPoint(const std::string &num) {
        std::string assumed_num{""};
        int loc = 0;
        if (num[0] == '-') {
            assumed_num += '-';
            loc = 1;
        }
        assumed_num += '.';

        for (char &c : num.substr(loc, num.length())) {
            switch (c) {
                case ' ':
                    break;
                case '-':
                    assumed_num += "E-";
                    break;
                case '+':
                    assumed_num += "E+";
                    break;
                default:
                    assumed_num += c;
            }
        }

        return std::stod(assumed_num);
    }

    //! === CHECKSUM ===
    int checksum(const std::string &navline) {
        int sum = 0;
        for (const char &c : navline) {
            switch (c) {
                case '-':
                    sum += 1;
                    break;
                case '1':
                    sum += 1;
                    break;
                case '2':
                    sum += 2;
                    break;
                case '3':
                    sum += 3;
                    break;
                case '4':
                    sum += 4;
                    break;
                case '5':
                    sum += 5;
                    break;
                case '6':
                    sum += 6;
                    break;
                case '7':
                    sum += 7;
                    break;
                case '8':
                    sum += 8;
                    break;
                case '9':
                    sum += 9;
                    break;
                default:
                    break;
            }
        }
        return sum % 10;
    }
};

}  // namespace satutils

#endif
