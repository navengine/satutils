/**
 * *satutils-python.cpp*
 *
 * =======  ========================================================================================
 * @file    src/satutils-python.cpp
 * @brief   PyBind11 wrapper for using satutils in python!
 * @date    March 2025
 * =======  ========================================================================================
 */

#include <pybind11/chrono.h>
#include <pybind11/eigen.h>
#include <pybind11/pybind11.h>
#include <pybind11/stl.h>

#include "satutils/atmosphere.hpp"
#include "satutils/code-gen.hpp"
#include "satutils/ephemeris.hpp"
#include "satutils/gnss-constants.hpp"
#include "satutils/gps-lnav.hpp"
#include "satutils/rinex-parser.hpp"
#include "satutils/time.hpp"
#include "satutils/tle-parser.hpp"

namespace py = pybind11;
using namespace satutils;

PYBIND11_MODULE(_satutils_core, h) {
  h.doc() = R"pbdoc(
      SatUtils
      ========

      A set of utilities common for both the simulation and processing of signals from various 
      satellite systems

      Contains the following submodules:

        1. `atmosphere`
        2. `codegen`
        3. `ephemeris`
        4. `gpslnav`
        5. `rinexparser`
        6. `time`
        6. `tleparser`
      )pbdoc";

  h.attr("__version__") = "1.0.0";

  //! === GNSS Constants ===========================================================================
  h.attr("GPS_PI") = GPS_PI<double>;
  h.attr("TWO_GPS_PI") = TWO_GPS_PI<double>;
  h.attr("GPS_L1_FREQUENCY") = GPS_L1_FREQUENCY<double>;
  h.attr("GPS_L2_FREQUENCY") = GPS_L2_FREQUENCY<double>;
  h.attr("GPS_L5_FREQUENCY") = GPS_L5_FREQUENCY<double>;
  h.attr("GPS_CA_CODE_RATE") = GPS_CA_CODE_RATE<double>;
  h.attr("GPS_L2_CODE_RATE") = GPS_L2_CODE_RATE<double>;
  h.attr("GPS_L5_CODE_RATE") = GPS_L5_CODE_RATE<double>;
  h.attr("GPS_CA_CODE_LENGTH") = GPS_CA_CODE_LENGTH;
  h.attr("GPS_L2CM_CODE_LENGTH") = GPS_L2CM_CODE_LENGTH;
  h.attr("GPS_L2CL_CODE_LENGTH") = GPS_L2CL_CODE_LENGTH;
  h.attr("GPS_L5_CODE_LENGTH") = GPS_L5_CODE_LENGTH;
  h.attr("LNAV_SUBFRAME_SIZE") = LNAV_SUBFRAME_SIZE;
  h.attr("LNAV_WORD_SIZE") = LNAV_WORD_SIZE;
  h.attr("LNAV_PREAMBLE_BITS") = LNAV_PREAMBLE_BITS;
  h.attr("LNAV_INV_PREAMBLE_BITS") = LNAV_INV_PREAMBLE_BITS;

  h.attr("GALILEO_E5_FREQUENCY") = GALILEO_E5_FREQUENCY<double>;
  h.attr("GALILEO_E5A_FREQUENCY") = GALILEO_E5A_FREQUENCY<double>;
  h.attr("GALILEO_E5B_FREQUENCY") = GALILEO_E5B_FREQUENCY<double>;
  h.attr("GALILEO_E6_FREQUENCY") = GALILEO_E6_FREQUENCY<double>;
  h.attr("GALILEO_E6_CODE_RATE") = GALILEO_E6_CODE_RATE<double>;
  h.attr("GALILEO_E1_DATA_RATE") = GALILEO_E1_DATA_RATE<double>;
  h.attr("GALILEO_E5A_DATA_RATE") = GALILEO_E5A_DATA_RATE<double>;
  h.attr("GALILEO_E5B_DATA_RATE") = GALILEO_E5B_DATA_RATE<double>;
  h.attr("GALILEO_E6_DATA_RATE") = GALILEO_E6_DATA_RATE<double>;
  h.attr("GALILEO_E1_CODE_LENGTH") = GALILEO_E1_CODE_LENGTH;
  h.attr("GALILEO_E5_CODE_LENGTH") = GALILEO_E5_CODE_LENGTH;
  h.attr("GALILEO_E6_CODE_LENGTH") = GALILEO_E6_CODE_LENGTH;

  h.attr("SGP_AE") = SGP_AE<double>;
  h.attr("SGP_XKMPER") = SGP_XKMPER<double>;
  h.attr("SGP_S") = SGP_S<double>;
  h.attr("SGP_QOMS2T") = SGP_QOMS2T<double>;
  h.attr("SGP_XKE") = SGP_XKE<double>;
  h.attr("SGP_CK2") = SGP_CK2<double>;
  h.attr("SGP_CK4") = SGP_CK4<double>;
  h.attr("SGP_A3OVK2") = SGP_A3OVK2<double>;

  //! === Atmosphere ===============================================================================
  py::module_ atm = h.def_submodule(
      "atmosphere",
      R"pbdoc(
      Atmosphere
      ==========
      
      GNSS atmospheric corrections.)pbdoc");

  // KlobucharElements
  py::class_<KlobucharElements<double>>(atm, "KlobucharElements")
      .def(py::init<>())
      .def_readwrite("a0", &KlobucharElements<double>::a0)
      .def_readwrite("a1", &KlobucharElements<double>::a1)
      .def_readwrite("a2", &KlobucharElements<double>::a2)
      .def_readwrite("a3", &KlobucharElements<double>::a3)
      .def_readwrite("b0", &KlobucharElements<double>::b0)
      .def_readwrite("b1", &KlobucharElements<double>::b1)
      .def_readwrite("b2", &KlobucharElements<double>::b2)
      .def_readwrite("b3", &KlobucharElements<double>::b3)
      .doc() = R"pbdoc(
               KlobucharElements
               =================

               Struct containing polynomial coefficients for ionospheric corrections
               )pbdoc";

  // IonoModel
  py::class_<IonoModel<double>>(atm, "IonoModel")
      .def(py::init<>())
      .def(py::init<const KlobucharElements<double> &>(), py::arg("klob"))
      .def(
          "SetKlobuchar",
          &IonoModel<double>::SetKlobuchar,
          py::arg("klob"),
          R"pbdoc(
          SetKlobuchar
          ============

          Set the Klobuchar elements

          Parameters
          ----------

          klob : KlobucharElements

              Klobuchar elements struct
          )pbdoc")
      .def(
          "GetKlobuchar",
          &IonoModel<double>::GetKlobuchar,
          R"pbdoc(
          GetKlobuchar
          ============

          Get the Klobuchar elements

          Returns
          -------

          klob : KlobucharElements

              Klobuchar elements struct
          )pbdoc")
      .def(
          "CalcIonoDelay",
          &IonoModel<double>::CalcIonoDelay,
          py::arg("Tow"),
          py::arg("lat"),
          py::arg("lon"),
          py::arg("az"),
          py::arg("el"),
          py::arg("gamma"),
          R"pbdoc(
          CalcIonoDelay
          ============

          Estimates the ionospheric delay based on the Klobuchar model

          Parameters
          ----------

          ToW : double

              GPS time of week in seconds [s]

          lat : double

              geodetic latitude [rad]

          lon : double

              geodetic longitude [rad]

          az : double

              azimuth angle to satellite [rad]

          el : double

              elevation angle to satellite [rad]

          gamma : double

              frequency dependant scaling factor (defaults to 1.0)

          Returns
          -------
          
          Iono : double

              ionospheric time delay [s]
          )pbdoc")
      .doc() = R"pbdoc(
               IonoModel
               =========
               
               Klobuchar based ionospheric delay model
               )pbdoc";

  // TropoModel
  py::class_<TropoModel<double>>(atm, "TropoModel")
      .def(py::init<>())
      .def(
          "CalcTropoDelay",
          &TropoModel<double>::CalcTropoDelay,
          py::arg("DoY"),
          py::arg("lat"),
          py::arg("h"),
          py::arg("el"),
          R"pbdoc(
          CalcTropoDelay
          ==============

          Estimates the tropospheric delay based on dry and wet delays

          Parameters
          ----------

          DoY : double

              current day of the year (Jan 1 = 0, Dec 31  = 365)

          lat : double

              geodetic latitude [rad]

          h : double

              geodetic altitude [m]

          el : double

              elevation angle to satellite [rad]

          Returns
          -------
          
          Tropo : double

              tropospheric time delay [s]
          )pbdoc")
      .doc() = R"pbdoc(
               TropoModel
               ==========
               
               Simple tropospheric model that does not use meterological data.
               )pbdoc";

  //! === Code-Gen =================================================================================
  py::module_ cgn = h.def_submodule(
      "codegen",
      R"pbdoc(
      Code-Gen
      ========
      
      Constellation gold code generators.)pbdoc");

  // CodeGenCA
  cgn.def(
      "CodeGenCA",
      py::overload_cast<bool[1023], uint8_t>(&CodeGenCA),
      py::arg("ca_code"),
      py::arg("prn_id"),
      R"pbdoc(
      CodeGenCA
      =========
      
      Generates the GPS L1 C/A gold code.

      Parameters
      ----------

      ca_code : np.ndarray

          size 1023 array of booleans

      prn_id : int

          Desired prn code to generate.

      Returns
      -------

      e : np.ndarray

          size 3 RPY euler angles [rad]
      )pbdoc");

  //! === Ephemeris ================================================================================
  py::module_ eph = h.def_submodule(
      "ephemeris",
      R"pbdoc(
      Ephemeris
      =========
      
      Satellite ephemeris structures.)pbdoc");

  // KeplerElements
  py::class_<KeplerElements<double>>(atm, "KeplerElements")
      .def(py::init<>())
      .def_readwrite("iode", &KeplerElements<double>::iode)
      .def_readwrite("iodc", &KeplerElements<double>::iodc)
      .def_readwrite("toe", &KeplerElements<double>::toe)
      .def_readwrite("tgd", &KeplerElements<double>::tgd)
      .def_readwrite("af2", &KeplerElements<double>::af2)
      .def_readwrite("af1", &KeplerElements<double>::af1)
      .def_readwrite("af0", &KeplerElements<double>::af0)
      .def_readwrite("e", &KeplerElements<double>::e)
      .def_readwrite("sqrtA", &KeplerElements<double>::sqrtA)
      .def_readwrite("deltan", &KeplerElements<double>::deltan)
      .def_readwrite("m0", &KeplerElements<double>::m0)
      .def_readwrite("omega0", &KeplerElements<double>::omega0)
      .def_readwrite("omega", &KeplerElements<double>::omega)
      .def_readwrite("omegaDot", &KeplerElements<double>::omegaDot)
      .def_readwrite("i0", &KeplerElements<double>::i0)
      .def_readwrite("iDot", &KeplerElements<double>::iDot)
      .def_readwrite("cuc", &KeplerElements<double>::cuc)
      .def_readwrite("cus", &KeplerElements<double>::cus)
      .def_readwrite("cic", &KeplerElements<double>::cic)
      .def_readwrite("cis", &KeplerElements<double>::cis)
      .def_readwrite("crc", &KeplerElements<double>::crc)
      .def_readwrite("crs", &KeplerElements<double>::crs)
      .def_readwrite("ura", &KeplerElements<double>::ura)
      .def_readwrite("health", &KeplerElements<double>::health)
      .doc() = R"pbdoc(
               KeplerElements
               ==============

               Ephemerides based on Keplerian orbital elements
               )pbdoc";

  // Sgp4Elements
  py::class_<Sgp4Elements<double>>(atm, "Sgp4Elements")
      .def(py::init<>())
      .def_readwrite("catalog_id", &Sgp4Elements<double>::catalog_id)
      .def_readwrite("week", &Sgp4Elements<double>::week)
      .def_readwrite("toe", &Sgp4Elements<double>::toe)
      .def_readwrite("Bstar", &Sgp4Elements<double>::Bstar)
      .def_readwrite("e0", &Sgp4Elements<double>::e0)
      .def_readwrite("omega0", &Sgp4Elements<double>::omega0)
      .def_readwrite("omega", &Sgp4Elements<double>::omega)
      .def_readwrite("i0", &Sgp4Elements<double>::i0)
      .def_readwrite("n0", &Sgp4Elements<double>::n0)
      .def_readwrite("nDot", &Sgp4Elements<double>::nDot)
      .def_readwrite("nDDot", &Sgp4Elements<double>::nDDot)
      .def_readwrite("m0", &Sgp4Elements<double>::m0)
      .doc() = R"pbdoc(
               Sgp4Elements
               ============

               Ephemerides based on SGP4 orbital elements
               )pbdoc";

  // KeplerEphem
  py::class_<KeplerEphem<double>>(atm, "KeplerEphem")
      .def(py::init<>())
      .def(py::init<const KeplerEphem<double> &>(), py::arg("eph"))
      .def(
          "SetEphemerides",
          &KeplerEphem<double>::SetEphemerides,
          py::arg("eph"),
          R"pbdoc(
          SetEphemerides
          ==============

          Set the Keplerian ephemeris elements

          Parameters
          ----------

          eph : KeplerElements

              Keplerian ephemeris elements struct
          )pbdoc")
      .def(
          "GetEphemerides",
          &KeplerEphem<double>::GetEphemerides,
          R"pbdoc(
          GetEphemerides
          ==============

          Get the Keplerian ephemeris elements

          Returns
          -------

          eph : KeplerElements

              Keplerian ephemeris elements struct
          )pbdoc")
      .def(
          "init",
          &KeplerEphem<double>::init,
          R"pbdoc(
          init
          ====

          Initialize additional ephemeris constants
          )pbdoc")
      .def(
          "CalcNavStates",
          &KeplerEphem<double>::CalcNavStates<true>,
          py::arg("clk"),
          py::arg("pos"),
          py::arg("vel"),
          py::arg("acc"),
          py::arg("transmit_time"),
          R"pbdoc(
          CalcIonoDelay
          ============

          Calculates satellite position, velocity, and acceleration using ephemeris

          Parameters
          ----------

          clk : np.ndarray

              Satellite clock corrections vector [s, s/s, s/s^2]

          pos : np.ndarray

              Satellite ECEF position vector [m]

          vel : np.ndarray

              Satellite ECEF velocity vector [m/s]

          acc : np.ndarray

              Satellite ECEF acceleration vector [m/s^2]

          transmit_time : np.ndarray

              GPS system/transmitter time (TOW) of the satellite (accounting for transit time from satellite to receiver) [gps seconds]
          )pbdoc")
      .doc() = R"pbdoc(
               KeplerEphem
               ===========
               
               Satellite navigation state estimator using Keplerian ephemeris elements
               )pbdoc";

  // Sgp4Ephem
  py::class_<Sgp4Ephem<double>>(atm, "Sgp4Ephem")
      .def(py::init<>())
      .def(py::init<const Sgp4Ephem<double> &>(), py::arg("eph"))
      .def(
          "SetEphemerides",
          &Sgp4Ephem<double>::SetEphemerides,
          py::arg("eph"),
          R"pbdoc(
          SetEphemerides
          ==============

          Set the SGP4 ephemeris elements

          Parameters
          ----------

          eph : Sgp4Elements

              SGP4 ephemeris elements struct
          )pbdoc")
      .def(
          "GetEphemerides",
          &Sgp4Ephem<double>::GetEphemerides,
          R"pbdoc(
          GetEphemerides
          ==============

          Get the SGP4 ephemeris elements

          Returns
          -------

          eph : Sgp4Elements

              SGP4 ephemeris elements struct
          )pbdoc")
      .def(
          "init",
          &Sgp4Ephem<double>::init,
          R"pbdoc(
          init
          ====

          Initialize additional ephemeris constants
          )pbdoc")
      .def(
          "CalcNavStates",
          &Sgp4Ephem<double>::CalcNavStates,
          py::arg("pos"),
          py::arg("vel"),
          py::arg("transmit_time"),
          R"pbdoc(
          CalcIonoDelay
          ============

          Calculates satellite position and velocity using SGP4 ephemeris

          Parameters
          ----------

          pos : np.ndarray

              Satellite ECEF position vector [m]

          vel : np.ndarray

              Satellite ECEF velocity vector [m/s]

          transmit_time : np.ndarray

              GPS system/transmitter time (TOW) of the satellite (accounting for transit time from satellite to receiver) [gps seconds]
          )pbdoc")
      .doc() = R"pbdoc(
               Sgp4Ephem
               ===========
               
               Satellite navigation state estimator using SGP4 ephemeris elements
               )pbdoc";

  //! === GPS-LNAV =================================================================================
  py::module_ lnav = h.def_submodule(
      "gpslnav",
      R"pbdoc(
      GPS-LNAV
      ========
      
      Implementation of GPS L1 C/A navigation message utils.)pbdoc");

  // GpsLnav
  py::class_<GpsLnav<double>>(atm, "GpsLnav")
      .def(py::init<>())
      .def(py::init<const KeplerElements<double> &>(), py::arg("eph"))
      .def(py::init<const KlobucharElements<double> &>(), py::arg("klob"))
      .def(
          py::init<const KeplerElements<double> &, const KlobucharElements<double> &>(),
          py::arg("eph"),
          py::arg("klob"))
      .def(
          "SetEphemerides",
          &GpsLnav<double>::SetEphemerides,
          py::arg("eph"),
          R"pbdoc(
          SetEphemerides
          ==============

          Set the Keplerian ephemeris elements

          Parameters
          ----------

          eph : KeplerElements

              Keplerian ephemeris elements struct
          )pbdoc")
      .def(
          "GetEphemerides",
          &GpsLnav<double>::GetEphemerides,
          R"pbdoc(
          GetEphemerides
          ==============

          Get the Keplerian ephemeris elements

          Returns
          -------

          eph : KeplerElements

              Keplerian ephemeris elements struct
          )pbdoc")
      .def(
          "SetKlobuchar",
          &GpsLnav<double>::SetKlobuchar,
          py::arg("klob"),
          R"pbdoc(
          SetKlobuchar
          ============

          Set the Klobuchar elements

          Parameters
          ----------

          klob : KlobucharElements

              Klobuchar elements struct
          )pbdoc")
      .def(
          "GetKlobuchar",
          &GpsLnav<double>::GetKlobuchar,
          R"pbdoc(
          GetKlobuchar
          ============

          Get the Klobuchar elements

          Returns
          -------

          klob : KlobucharElements

              Klobuchar elements struct
          )pbdoc")
      .def(
          "GetWeekNumber",
          &GpsLnav<double>::GetWeekNumber,
          R"pbdoc(
          GetWeekNumber
          =============

          Gets current GPS week number

          Returns
          -------

          week : int

              GPS Week number
          )pbdoc")
      .def(
          "GetTimeOfWeek",
          &GpsLnav<double>::GetTimeOfWeek,
          R"pbdoc(
          GetTimeOfWeek
          =============

          Gets current GPS time of week

          Returns
          -------

          ToW : double

              GPS time of week [s]
          )pbdoc")
      .def(
          "AreEphemeridesParsed",
          &GpsLnav<double>::AreEphemeridesParsed,
          R"pbdoc(
          AreEphemeridesParsed
          ====================

          Check if subframes 1, 2, and 3 have been parsed

          Returns
          -------

          status : bool

              True|False based on if subframe 1,2 and 3 have been parsed
          )pbdoc")
      .def(
          "SetNextBit",
          &GpsLnav<double>::SetNextBit,
          py::arg("bit"),
          R"pbdoc(
          SetNextBit
          ==========

          Read the next navigation data bit and triggers ephemeris parsing attempts

          Parameters
          ----------

          bit : bool

              Navigation data bit parsed by tracking loop

          Returns
          -------

          status : bool

              True|False based on if a subframe was successfully parsed
          )pbdoc")
      .def(
          "ParseSubframe",
          &GpsLnav<double>::ParseSubframe,
          R"pbdoc(
          ParseSubframe
          =============

          Attempts to parse subframe

          Returns
          -------

          status : bool

              True|False based on if a subframe was successfully parsed
          )pbdoc")
      .def(
          "ParityCheck",
          &GpsLnav<double>::ParityCheck,
          py::arg("gpsword"),
          py::arg("D29star"),
          py::arg("D30star"),
          R"pbdoc(
          ParityCheck
          ===========

          GPS data bit parity check (IS-GPS-200N pg. 139)

          Parameters
          ----------

          gpsword : uint32

              GPS word to evaluate

          D29star : bool

              29th bit of previous gps word

          D30star : bool

              30th bit of previous gps word

          Returns
          -------

          status : bool

              Parity success or failure
          )pbdoc")
      .def(
          "LoadPreamble",
          &GpsLnav<double>::LoadPreamble,
          R"pbdoc(
          LoadPreamble
          ============

          Reads words 1 and 2 of each subframe (IS-GPS-200N pg. 92)
          )pbdoc")
      .def(
          "LoadSubframe1",
          &GpsLnav<double>::LoadSubframe1,
          R"pbdoc(
          LoadSubframe1
          =============

          Reads GPS LNAV subframe 1 (IS-GPS-200N pg. 80 & 97)
          )pbdoc")
      .def(
          "LoadSubframe2",
          &GpsLnav<double>::LoadSubframe2,
          R"pbdoc(
          LoadSubframe2
          =============

          Reads GPS LNAV subframe 2 (IS-GPS-200N pg. 81 & 105)
          )pbdoc")
      .def(
          "LoadSubframe3",
          &GpsLnav<double>::LoadSubframe3,
          R"pbdoc(
          LoadSubframe3
          =============

          Reads GPS LNAV subframe 3 (IS-GPS-200N pg. 82 & 105)
          )pbdoc");
  // .def(
  //     "LoadSubframe4",
  //     &GpsLnav<double>::LoadSubframe4,
  //     R"pbdoc(
  //     LoadSubframe4
  //     =============

  //     Reads GPS LNAV subframe 4 (IS-GPS-200N pg. )
  //     )pbdoc")
  // .def(
  //     "LoadSubframe5",
  //     &GpsLnav<double>::LoadSubframe5,
  //     R"pbdoc(
  //     LoadSubframe5
  //     =============

  //     Reads GPS LNAV subframe 5 (IS-GPS-200N pg. )
  //     )pbdoc");

  //! === Rinex Parser =============================================================================
  py::module_ rnx = h.def_submodule(
      "rinexparser",
      R"pbdoc(
      Rinex Parser
      ============
      
      Tool to parse Rinex 3.0 file into ephemeris elements.)pbdoc");

  // ParseTimeBlock
  rnx.def(
      "ParseTimeBlock",
      &ParseTimeBlock<double>,
      py::arg("week"),
      py::arg("sec"),
      py::arg("line"),
      R"pbdoc(
      ParseTimeBlock
      ==============

      Parse the time of a rinex block

      Parameters
      ----------

      week : int

          GPS week number

      sec : double

          GPS second

      line : string

          String containing the time information from the rinex block
      )pbdoc");

  // ParseNavBlock
  rnx.def(
      "ParseNavBlock",
      &ParseNavBlock<double>,
      py::arg("line"),
      R"pbdoc(
      ParseNavBlock
      =============

      Parse the navigation data of a rinex block to a vector

      Parameters
      ----------

      line : string

          String containing the navigation information from the rinex block

      Returns
      -------

      blk : np.ndarray

          Vector of rinex navigation data
      )pbdoc");

  // RinexParser
  rnx.def(
      "RinexParser",
      &RinexParser<double>,
      py::arg("filename"),
      R"pbdoc(
        RinexParser
        ===========

        Parse a Rinex file and return its contents

        Parameters
        ----------

        filename : string

            Filepath to rinex.

        Returns
        -------

        x : dict

            Parsed rinex information.
        )pbdoc");

  //! === Time =====================================================================================
  py::module_ time = h.def_submodule(
      "time",
      R"pbdoc(
      Time
      ====
      
      GNSS time conversions.)pbdoc");

  time.attr("HALF_WEEK") = HALF_WEEK<double>;
  time.attr("WEEK") = WEEK<double>;
  time.attr("S_PER_DAY") = S_PER_DAY<double>;

  // Date2GpsTime
  time.def(
      "Date2GpsTime",
      &Date2GpsTime<int, double>,
      py::arg("week"),
      py::arg("sec"),
      py::arg("dt"),
      R"pbdoc(
      Date2GpsTime
      ============

      Convert a date-time into a GPS week and second of week

      Parameters
      ----------

      week : int

          GPS week number

      sec : double

          GPS second

      dt : datetime

          System clock time point
      )pbdoc");

  // CheckGpsSecond
  time.def(
      "CheckGpsSecond",
      &CheckGpsSecond<double>,
      py::arg("t"),
      R"pbdoc(
      CheckGpsSecond
      ==============

      Ensure that the time is within a GPS week

      Parameters
      ----------

      t : double

          Time to check [gps seconds]

      Returns
      -------

      t : double

          Corrected time [gps seconds]
      )pbdoc");

  //! === TLE Parser ===============================================================================
  py::module_ tle = h.def_submodule(
      "tleparser",
      R"pbdoc(
      TLE Parser
      ==========
      
      Tool to parse tle file into SGP4 elements.)pbdoc");

  // CheckSum
  tle.def(
      "CheckSum",
      &CheckSum,
      py::arg("line"),
      R"pbdoc(
      CheckSum
      =======

      Calculate the checksum of the TLE line

      Parameters
      ----------

      line : string

          An information line from a TLE

      Returns
      -------

      x : int

          Calculated checksum
      )pbdoc");

  // AssumedDecimalPoint
  tle.def(
      "AssumedDecimalPoint",
      &AssumedDecimalPoint<double>,
      py::arg("assumed_num"),
      py::arg("in_str"),
      R"pbdoc(
        AssumedDecimalPoint
        ===================
  
        Places an assumed decimal point at the begining of a number
  
        Parameters
        ----------

        assumed_num : double

            Resulting floating point number
  
        in_str : string
  
            String containing the number
        )pbdoc");

  // TleParser
  tle.def(
      "TleParser",
      &TleParser<double>,
      py::arg("filename"),
      R"pbdoc(
        TleParser
        =========

        Parse a TLE file and return its contents

        Parameters
        ----------

        filename : string

            Filepath to TLE.

        Returns
        -------

        x : dict

            Parsed TLE information.
        )pbdoc");
}