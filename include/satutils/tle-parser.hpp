/**
 * *tle-parser.hpp*
 *
 * =======  ========================================================================================
 * @file    satutils/tle-parser.hpp
 * @brief   Tool to parse tle file into SGP4 elements.
 * @date    January 2025
 * @author  Daniel Sturdivant <Auburn University GAVLAB>
 * @author  Blake Baker <Auburn University GAVLAB>
 * =======  ========================================================================================
 */

// TODO: figure out a nice-easy way to map all satellite norad id to strings (celestrak?)

#ifndef SATUTILS_TLE_PARSER_HPP
#define SATUTILS_TLE_PARSER_HPP

#include <chrono>
#include <fstream>
#include <iostream>
#include <map>
#include <sstream>

#include "satutils/ephemeris.hpp"
#include "satutils/time.hpp"

namespace satutils {

/**
 * *=== CheckSum ===*
 * @brief calculate the checksum of the line
 */
int CheckSum(const std::string &navline) {
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
};

/**
 * *=== AssumedDecimalPoint ===*
 * @brief places an assumed decimal point at the begining of a number
 */
template <typename T>
void AssumedDecimalPoint(T &assumed_num, const std::string &num) {
  std::stringstream ss;
  int loc = 0;
  if (num[0] == '-') {
    ss << '-';
    loc = 1;
  }
  ss << '.';

  for (char &c : num.substr(loc, num.length())) {
    switch (c) {
      case ' ':
        break;
      case '-':
        ss << "E-";
        break;
      case '+':
        ss << "E+";
        break;
      default:
        ss << c;
    }
  }

  ss >> assumed_num;
};

/**
 * *=== TleParser ===*
 * @brief Parse a TLE file and return its contents
 */
template <typename T>
std::map<std::string, Sgp4Elements<T>> TleParser(std::string filename) {
  // safely open file
  std::ifstream fid = std::ifstream(filename);
  if (!fid.is_open()) {
    std::cout << "satutils::TleParser - Invalid file!\n";
  } else if (fid.bad()) {
    std::cout << "satutils::TleParser - Failed to read file!\n";
  }
  fid.seekg(0);

  // initialize parsing variables
  std::string line, sv_id, word;
  std::chrono::system_clock::time_point tp;
  T catalog_id, week, toe, nDot, nDDot, Bstar, i0, omega0, e0, omega, m0, n0;
  int year, day;
  uint64_t ns;
  std::stringstream ss;
  std::map<std::string, Sgp4Elements<T>> my_map;
  Sgp4Elements<T> my_eph;

  // read file
  while (!fid.eof()) {
    // every 3 lines is a satellite
    for (int i = 0; i < 3; i++) {
      std::getline(fid, line, '\n');
      switch (i) {
        case 0:
          // satellite id (2 words)
          ss << line;
          while (std::getline(ss, word, ' ')) {
            if (word.length() > 0) {
              sv_id += word;
              sv_id += " ";
            }
          }
          sv_id = sv_id.substr(0, sv_id.length() - 1);
          break;
        case 1:
          catalog_id = std::stoi(line.substr(2, 5));
          year = std::stoi(line.substr(18, 2)) + 2000;
          day = std::stoi(line.substr(20, 3));
          tp = std::chrono::sys_days{std::chrono::January / 0 / std::chrono::year{year}};
          tp += std::chrono::days{day};
          ns = static_cast<uint64_t>(1e9 * 86400.0 * std::stod(line.substr(23, 9)));
          tp += std::chrono::nanoseconds{ns};
          Date2GpsTime(week, toe, tp);
          nDot = std::stod(line.substr(33, 10));
          AssumedDecimalPoint<T>(nDDot, line.substr(44, 8));
          AssumedDecimalPoint<T>(Bstar, line.substr(53, 8));
          if (CheckSum(line.substr(0, line.length() - 1)) != std::stoi(line.substr(68, 1))) {
            std::cout << "satutils::TleParser - ERROR " << sv_id << " line1 checksum incorrect\n";
          }
          break;
        case 2:
          i0 = std::stod(line.substr(8, 8));
          omega0 = std::stod(line.substr(17, 8));
          AssumedDecimalPoint<T>(e0, line.substr(26, 7));
          omega = std::stod(line.substr(34, 8));
          m0 = std::stod(line.substr(43, 8));
          n0 = std::stod(line.substr(52, 11));
          if (CheckSum(line.substr(0, line.length() - 1)) != std::stoi(line.substr(68, 1))) {
            std::cout << "satutils::TleParser - ERROR " << sv_id << " line2 checksum incorrect\n";
          }
          break;
      }
      ss.str("");
      ss.clear();
    }

    // std::cout << "sv_id: " << sv_id << "\ncatalog_id: " << catalog_id << "\nweek: " << week
    //           << "\ntoe: " << toe << "\nBstar: " << Bstar << "\ne0: " << e0
    //           << "\nomega0: " << omega0 << "\noemga: " << omega << "\ni0: " << i0 << "\nn0: " <<
    //           n0
    //           << "\nnDot: " << nDot << "\nnDDot: " << nDDot << "\nm0: " << m0 << "\n";

    // add item to map
    my_eph = {
        catalog_id,
        week,
        toe,
        Bstar,
        e0,
        navtools::DEG2RAD<T> * omega0,
        navtools::DEG2RAD<T> * omega,
        navtools::DEG2RAD<T> * i0,
        navtools::TWO_PI<T> / 1440.0 * n0,
        navtools::TWO_PI<T> / 2073600.0 * nDot,
        navtools::TWO_PI<T> / 2985984000.0 * nDDot,
        navtools::DEG2RAD<T> * m0};
    my_map.insert({sv_id, my_eph});
  }

  // done
  fid.close();
  return my_map;
};

}  // namespace satutils

#endif