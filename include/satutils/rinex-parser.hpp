/**
 * *rinex-parser.hpp*
 *
 * =======  ========================================================================================
 * @file    satutils/rinex-parser.hpp
 * @brief   Tool to parse Rinex 3.0 file into ephemeris elements.
 * @date    January 2025
 * @author  Daniel Sturdivant <Auburn University GAVLAB>
 * @author  Blake Baker <Auburn University GAVLAB>
 * =======  ========================================================================================
 */

// TODO: add base ECEF based parser for Glonass constellation

#ifndef SATUTILS_RINEX_PARSER_HPP
#define SATUTILS_RINEX_PARSER_HPP

#include <chrono>
#include <fstream>
#include <iostream>
#include <map>
#include <sstream>
#include <string>

#include "satutils/atmosphere.hpp"
#include "satutils/ephemeris.hpp"
#include "satutils/time.hpp"

namespace satutils {

/**
 * *=== ParseTimeBlock ===*
 * @brief
 */
template <typename T>
void ParseTimeBlock(int &week, T &time, const std::string &timeline)
{
  int year = std::stoi(timeline.substr(0, 4));
  unsigned int month = std::stoi(timeline.substr(5, 2));
  unsigned int day = std::stoi(timeline.substr(8, 2));
  unsigned int hour = std::stoi(timeline.substr(11, 2));
  unsigned int min = std::stoi(timeline.substr(14, 2));
  unsigned int sec = std::stoi(timeline.substr(17, 2));
  std::chrono::system_clock::time_point tp = std::chrono::sys_days{
      std::chrono::month{month} / std::chrono::day{day} / std::chrono::year{year}};
  tp += std::chrono::hours{hour} + std::chrono::minutes{min} + std::chrono::seconds{sec};
  Date2GpsTime(week, time, tp);
}

/**
 * *=== ParseNavBlock ===*
 * @brief
 */
template <typename T>
std::vector<T> ParseNavBlock(const std::string &navblock)
{
  std::vector<T> navdata;
  std::string navword;
  std::istringstream iss;

  // split word every 19 characters
  for (size_t i = 0; i < navblock.length(); i += 19) {
    navword = navblock.substr(i, 19);
    if (navword.find_first_not_of(' ') == std::string::npos) continue;

    if (navword.find('D') != std::string::npos) {
      std::replace(navword.begin(), navword.end(), 'D', 'E');
    }

    navdata.push_back(std::stod(navword));
  }
  return navdata;
}

/**
 * *=== RinexParser ===*
 * @brief
 */
template <typename T>
std::map<std::string, std::pair<KlobucharElements<T>, KeplerElements<T>>> RinexParser(
    std::string filename)
{
  // safely open file
  std::ifstream fid = std::ifstream(filename);
  if (!fid.is_open()) {
    std::cout << "satutils::RinexParser - Invalid file!\n";
  } else if (fid.bad()) {
    std::cout << "satutils::RinexParser - Failed to read file!\n";
  }
  fid.seekg(0);

  // init parsing variables
  const std::string GPS_ALPHA_TOKEN = "GPSA";
  const std::string GPS_BETA_TOKEN = "GPSB";
  const std::string GALILEO_IONO_TOKEN = "GAL";
  const std::string IONO_TOKEN = "IONOSPHERIC CORR";
  const std::string TIME_TOKEN = "LEAP SECONDS";
  const std::string EOH_TOKEN = "END OF HEADER";
  const std::string COMMENT_TOKEN = "COMMENT";
  // int LEAP_SECONDS = 0;
  int week;
  T tmp_val, toc;
  KlobucharElements<T> GPS_IONO, GALILEO_IONO;
  KeplerElements<T> TMP_EPH;
  std::string line, sv_id, navline;
  std::map<std::string, std::pair<KlobucharElements<T>, KeplerElements<T>>> my_map;

  // --- Parse Header ---
  while (!fid.eof()) {
    std::getline(fid, line, '\n');

    // skip comment lines
    if (line.find(COMMENT_TOKEN) != std::string::npos) {
      continue;
    }

    // check for ionospheric parameters
    else if (line.find(IONO_TOKEN) != std::string::npos) {
      std::istringstream iss(line.substr(6, 47));
      if (line.find(GPS_ALPHA_TOKEN) != std::string::npos) {
        // found gps klobuchar alpha parameters
        int i = 0;
        //TODO exponent could be marked with 'D' instead of 'E'
        while (iss >> tmp_val) {
          switch (i) {
            case 0:
              GPS_IONO.a0 = tmp_val;
              break;
            case 1:
              GPS_IONO.a1 = tmp_val;
              break;
            case 2:
              GPS_IONO.a2 = tmp_val;
              break;
            case 3:
              GPS_IONO.a3 = tmp_val;
              break;
            default:
              break;
          }
          i++;
        }
      }
      else if (line.find(GPS_BETA_TOKEN) != std::string::npos) {
        // found gps klobuchar beta parameters
        int i = 0;
        while (iss >> tmp_val) {
          switch (i) {
            case 0:
              GPS_IONO.b0 = tmp_val;
              break;
            case 1:
              GPS_IONO.b1 = tmp_val;
              break;
            case 2:
              GPS_IONO.b2 = tmp_val;
              break;
            case 3:
              GPS_IONO.b3 = tmp_val;
              break;
            default:
              break;
          }
          i++;
        }
      }
      else if (line.find(GALILEO_IONO_TOKEN) != std::string::npos) {
        // found galileo ionosphere alpha parameters
        int i = 0;
        while (iss >> tmp_val) {
          switch (i) {
            case 0:
              GALILEO_IONO.a0 = tmp_val;
              break;
            case 1:
              GALILEO_IONO.a1 = tmp_val;
              break;
            case 2:
              GALILEO_IONO.a2 = tmp_val;
              break;
            case 3:
              GALILEO_IONO.a3 = tmp_val;
              break;
            default:
              break;
          }
          i++;
        }
      }
      continue;
    }

    // check for UTC leap second correction
    // if (line.find(TIME_TOKEN) != std::string::npos) {
    //   LEAP_SECONDS = std::stoi(line.substr(4, 6));
    //   continue;
    // }

    // check for end of header
    else if (line.find(EOH_TOKEN) != std::string::npos) {
      break;
    }
  }

  // --- Parse Rinex Entries ---
  while (!fid.eof()) {
    std::getline(fid, line, '\n');
    if (line.size() < 4) continue;
    sv_id = line.substr(0, 3);
    ParseTimeBlock<T>(week, toc, line.substr(4, 20));
    navline = line.substr(23, line.length());
    //TODO this assumes the same number of lines for each constellation type
    for (int i = 1; i < 8; i++) {
      std::getline(fid, line, '\n');
      navline += line.substr(4, line.length());
    }
    std::vector<T> nav = ParseNavBlock<T>(navline);

    if (sv_id[0] == 'G') {
      // GPS satellite data
      // TMP_EPH.week = week;                              // nav[21]
      // TMP_EPH.l2c_flag = static_cast<bool>(nav[20]);  //
      // TMP_EPH.l2p_flag = static_cast<bool>(nav[22]);  //
      TMP_EPH.iode = nav[3];       // [s]
      TMP_EPH.iodc = nav[26];      // [s]
      TMP_EPH.toe = nav[11];       // [s]
      TMP_EPH.toc = toc;           // [s]
      TMP_EPH.tgd = nav[25];       // [s]
      TMP_EPH.af2 = nav[2];        // [s/s^2]
      TMP_EPH.af1 = nav[1];        // [s/s]
      TMP_EPH.af0 = nav[0];        // [s]
      TMP_EPH.e = nav[8];          //
      TMP_EPH.sqrtA = nav[10];     // [m]
      TMP_EPH.deltan = nav[5];     // [rad/s]
      TMP_EPH.m0 = nav[6];         // [rad]
      TMP_EPH.omega0 = nav[13];    // [rad]
      TMP_EPH.omega = nav[17];     // [rad]
      TMP_EPH.omegaDot = nav[18];  // [rad/s]
      TMP_EPH.i0 = nav[15];        // [rad]
      TMP_EPH.iDot = nav[19];      // [rad/s]
      TMP_EPH.cuc = nav[7];        // [rad]
      TMP_EPH.cus = nav[9];        // [rad]
      TMP_EPH.cic = nav[12];       // [rad]
      TMP_EPH.cis = nav[14];       // [rad]
      TMP_EPH.crc = nav[16];       // [m]
      TMP_EPH.crs = nav[4];        // [m]
      TMP_EPH.ura = nav[23];       // [m]
      TMP_EPH.health = nav[24];
      my_map.insert({sv_id, std::make_pair(GPS_IONO, TMP_EPH)});
    }
    else if (sv_id[0] == 'E') {
      // Galileo satellite data

      // eph.id = sv_id;
      // TMP_EPH.leap_seconds = h_.leap_seconds;
      // TMP_EPH.data_source_flag = static_cast<uint16_t>(nav[20]);  //
      // TMP_EPH.BGD_e5a_e1 = nav[24];   // [s]
      // TMP_EPH.BGD_e5b_e1 = nav[25];   // [s]
      // TMP_EPH.week = nav[21];           //

      TMP_EPH.iode = nav[3];       // [s]
      TMP_EPH.iodc = nav[3];       // [s]
      TMP_EPH.toe = nav[11];       // [s]
      TMP_EPH.toc = toc;           // [s]
      TMP_EPH.tgd = nav[24];       // [s]
      TMP_EPH.af2 = nav[2];        // [s/s^2]
      TMP_EPH.af1 = nav[1];        // [s/s]
      TMP_EPH.af0 = nav[0];        // [s]
      TMP_EPH.e = nav[8];          //
      TMP_EPH.sqrtA = nav[10];     // [m]
      TMP_EPH.deltan = nav[5];     // [rad/s]
      TMP_EPH.m0 = nav[6];         // [rad]
      TMP_EPH.omega0 = nav[13];    // [rad]
      TMP_EPH.omega = nav[17];     // [rad]
      TMP_EPH.omegaDot = nav[18];  // [rad/s]
      TMP_EPH.i0 = nav[15];        // [rad]
      TMP_EPH.iDot = nav[19];      // [rad/s]
      TMP_EPH.cuc = nav[7];        // [rad]
      TMP_EPH.cus = nav[9];        // [rad]
      TMP_EPH.cic = nav[12];       // [rad]
      TMP_EPH.cis = nav[14];       // [rad]
      TMP_EPH.crc = nav[16];       // [m]
      TMP_EPH.crs = nav[4];        // [m]
      TMP_EPH.ura = nav[22];       // [m]
      TMP_EPH.health = nav[23];    //
      my_map.insert({sv_id, std::make_pair(GALILEO_IONO, TMP_EPH)});
    }
  }

  return my_map;
};


}  // namespace satutils

#endif
