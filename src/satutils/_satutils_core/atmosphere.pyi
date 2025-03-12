"""

      Atmosphere
      ==========
      
      GNSS atmospheric corrections.
"""
from __future__ import annotations
import numpy
import typing
__all__ = ['GpsLnav', 'IonoModel', 'KeplerElements', 'KeplerEphem', 'KlobucharElements', 'Sgp4Elements', 'Sgp4Ephem', 'TropoModel']
class GpsLnav:
    def AreEphemeridesParsed(self) -> bool:
        """
                  AreEphemeridesParsed
                  ====================
        
                  Check if subframes 1, 2, and 3 have been parsed
        
                  Returns
                  -------
        
                  status : bool
        
                      True|False based on if subframe 1,2 and 3 have been parsed
        """
    def GetEphemerides(self) -> KeplerElements:
        """
                  GetEphemerides
                  ==============
        
                  Get the Keplerian ephemeris elements
        
                  Returns
                  -------
        
                  eph : KeplerElements
        
                      Keplerian ephemeris elements struct
        """
    def GetKlobuchar(self) -> KlobucharElements:
        """
                  GetKlobuchar
                  ============
        
                  Get the Klobuchar elements
        
                  Returns
                  -------
        
                  klob : KlobucharElements
        
                      Klobuchar elements struct
        """
    def GetTimeOfWeek(self) -> float:
        """
                  GetTimeOfWeek
                  =============
        
                  Gets current GPS time of week
        
                  Returns
                  -------
        
                  ToW : double
        
                      GPS time of week [s]
        """
    def GetWeekNumber(self) -> int:
        """
                  GetWeekNumber
                  =============
        
                  Gets current GPS week number
        
                  Returns
                  -------
        
                  week : int
        
                      GPS Week number
        """
    def LoadPreamble(self) -> None:
        """
                  LoadPreamble
                  ============
        
                  Reads words 1 and 2 of each subframe (IS-GPS-200N pg. 92)
        """
    def LoadSubframe1(self) -> None:
        """
                  LoadSubframe1
                  =============
        
                  Reads GPS LNAV subframe 1 (IS-GPS-200N pg. 80 & 97)
        """
    def LoadSubframe2(self) -> None:
        """
                  LoadSubframe2
                  =============
        
                  Reads GPS LNAV subframe 2 (IS-GPS-200N pg. 81 & 105)
        """
    def LoadSubframe3(self) -> None:
        """
                  LoadSubframe3
                  =============
        
                  Reads GPS LNAV subframe 3 (IS-GPS-200N pg. 82 & 105)
        """
    def ParityCheck(self, gpsword: int, D29star: bool, D30star: bool) -> bool:
        """
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
        """
    def ParseSubframe(self) -> bool:
        """
                  ParseSubframe
                  =============
        
                  Attempts to parse subframe
        
                  Returns
                  -------
        
                  status : bool
        
                      True|False based on if a subframe was successfully parsed
        """
    def SetEphemerides(self, eph: KeplerElements) -> None:
        """
                  SetEphemerides
                  ==============
        
                  Set the Keplerian ephemeris elements
        
                  Parameters
                  ----------
        
                  eph : KeplerElements
        
                      Keplerian ephemeris elements struct
        """
    def SetKlobuchar(self, klob: KlobucharElements) -> None:
        """
                  SetKlobuchar
                  ============
        
                  Set the Klobuchar elements
        
                  Parameters
                  ----------
        
                  klob : KlobucharElements
        
                      Klobuchar elements struct
        """
    def SetNextBit(self, bit: bool) -> bool:
        """
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
        """
    @typing.overload
    def __init__(self) -> None:
        ...
    @typing.overload
    def __init__(self, eph: KeplerElements) -> None:
        ...
    @typing.overload
    def __init__(self, klob: KlobucharElements) -> None:
        ...
    @typing.overload
    def __init__(self, eph: KeplerElements, klob: KlobucharElements) -> None:
        ...
class IonoModel:
    """
    
                   IonoModel
                   =========
                   
                   Klobuchar based ionospheric delay model
                   
    """
    def CalcIonoDelay(self, Tow: float, lat: float, lon: float, az: float, el: float, gamma: float) -> float:
        """
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
        """
    def GetKlobuchar(self) -> KlobucharElements:
        """
                  GetKlobuchar
                  ============
        
                  Get the Klobuchar elements
        
                  Returns
                  -------
        
                  klob : KlobucharElements
        
                      Klobuchar elements struct
        """
    def SetKlobuchar(self, klob: KlobucharElements) -> None:
        """
                  SetKlobuchar
                  ============
        
                  Set the Klobuchar elements
        
                  Parameters
                  ----------
        
                  klob : KlobucharElements
        
                      Klobuchar elements struct
        """
    @typing.overload
    def __init__(self) -> None:
        ...
    @typing.overload
    def __init__(self, klob: KlobucharElements) -> None:
        ...
class KeplerElements:
    """
    
                   KeplerElements
                   ==============
    
                   Ephemerides based on Keplerian orbital elements
                   
    """
    af0: float
    af1: float
    af2: float
    cic: float
    cis: float
    crc: float
    crs: float
    cuc: float
    cus: float
    deltan: float
    e: float
    health: float
    i0: float
    iDot: float
    iodc: float
    iode: float
    m0: float
    omega: float
    omega0: float
    omegaDot: float
    sqrtA: float
    tgd: float
    toe: float
    ura: float
    def __init__(self) -> None:
        ...
class KeplerEphem:
    """
    
                   KeplerEphem
                   ===========
                   
                   Satellite navigation state estimator using Keplerian ephemeris elements
                   
    """
    def CalcNavStates(self, clk: numpy.ndarray[numpy.float64[3, 1], numpy.ndarray.flags.writeable], pos: numpy.ndarray[numpy.float64[3, 1], numpy.ndarray.flags.writeable], vel: numpy.ndarray[numpy.float64[3, 1], numpy.ndarray.flags.writeable], acc: numpy.ndarray[numpy.float64[3, 1], numpy.ndarray.flags.writeable], transmit_time: float) -> None:
        """
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
        """
    def GetEphemerides(self) -> KeplerElements:
        """
                  GetEphemerides
                  ==============
        
                  Get the Keplerian ephemeris elements
        
                  Returns
                  -------
        
                  eph : KeplerElements
        
                      Keplerian ephemeris elements struct
        """
    def SetEphemerides(self, eph: KeplerElements) -> None:
        """
                  SetEphemerides
                  ==============
        
                  Set the Keplerian ephemeris elements
        
                  Parameters
                  ----------
        
                  eph : KeplerElements
        
                      Keplerian ephemeris elements struct
        """
    @typing.overload
    def __init__(self) -> None:
        ...
    @typing.overload
    def __init__(self, eph: KeplerEphem) -> None:
        ...
    def init(self) -> None:
        """
                  init
                  ====
        
                  Initialize additional ephemeris constants
        """
class KlobucharElements:
    """
    
                   KlobucharElements
                   =================
    
                   Struct containing polynomial coefficients for ionospheric corrections
                   
    """
    a0: float
    a1: float
    a2: float
    a3: float
    b0: float
    b1: float
    b2: float
    b3: float
    def __init__(self) -> None:
        ...
class Sgp4Elements:
    """
    
                   Sgp4Elements
                   ============
    
                   Ephemerides based on SGP4 orbital elements
                   
    """
    Bstar: float
    catalog_id: float
    e0: float
    i0: float
    m0: float
    n0: float
    nDDot: float
    nDot: float
    omega: float
    omega0: float
    toe: float
    week: float
    def __init__(self) -> None:
        ...
class Sgp4Ephem:
    """
    
                   Sgp4Ephem
                   ===========
                   
                   Satellite navigation state estimator using SGP4 ephemeris elements
                   
    """
    def CalcNavStates(self, pos: numpy.ndarray[numpy.float64[3, 1]], vel: numpy.ndarray[numpy.float64[3, 1]], transmit_time: float) -> None:
        """
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
        """
    def GetEphemerides(self) -> Sgp4Elements:
        """
                  GetEphemerides
                  ==============
        
                  Get the SGP4 ephemeris elements
        
                  Returns
                  -------
        
                  eph : Sgp4Elements
        
                      SGP4 ephemeris elements struct
        """
    def SetEphemerides(self, eph: Sgp4Elements) -> None:
        """
                  SetEphemerides
                  ==============
        
                  Set the SGP4 ephemeris elements
        
                  Parameters
                  ----------
        
                  eph : Sgp4Elements
        
                      SGP4 ephemeris elements struct
        """
    @typing.overload
    def __init__(self) -> None:
        ...
    @typing.overload
    def __init__(self, eph: Sgp4Ephem) -> None:
        ...
    def init(self) -> None:
        """
                  init
                  ====
        
                  Initialize additional ephemeris constants
        """
class TropoModel:
    """
    
                   TropoModel
                   ==========
                   
                   Simple tropospheric model that does not use meterological data.
                   
    """
    def CalcTropoDelay(self, DoY: float, lat: float, h: float, el: float) -> float:
        """
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
        """
    def __init__(self) -> None:
        ...
