"""

Ephemeris
=========

Satellite ephemeris structures.
"""

from __future__ import annotations
import numpy
import typing

__all__ = ["KeplerElements", "KeplerEphem", "Sgp4Elements", "Sgp4Ephem"]

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
    toc: float
    toe: float
    ura: float
    def __init__(self) -> None: ...

class KeplerEphem:
    """

    KeplerEphem
    ===========

    Satellite navigation state estimator using Keplerian ephemeris elements

    """

    def CalcNavStates(
        self,
        clk: numpy.ndarray[numpy.float64[3, 1], numpy.ndarray.flags.writeable],
        pos: numpy.ndarray[numpy.float64[3, 1], numpy.ndarray.flags.writeable],
        vel: numpy.ndarray[numpy.float64[3, 1], numpy.ndarray.flags.writeable],
        acc: numpy.ndarray[numpy.float64[3, 1], numpy.ndarray.flags.writeable],
        transmit_time: float,
    ) -> None:
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
    def __init__(self) -> None: ...
    @typing.overload
    def __init__(self, eph: KeplerEphem) -> None: ...
    def init(self) -> None:
        """
        init
        ====

        Initialize additional ephemeris constants
        """

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
    def __init__(self) -> None: ...

class Sgp4Ephem:
    """

    Sgp4Ephem
    ===========

    Satellite navigation state estimator using SGP4 ephemeris elements

    """

    def CalcNavStates(
        self,
        pos: numpy.ndarray[numpy.float64[3, 1]],
        vel: numpy.ndarray[numpy.float64[3, 1]],
        transmit_time: float,
    ) -> None:
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
    def __init__(self) -> None: ...
    @typing.overload
    def __init__(self, eph: Sgp4Ephem) -> None: ...
    def init(self) -> None:
        """
        init
        ====

        Initialize additional ephemeris constants
        """
