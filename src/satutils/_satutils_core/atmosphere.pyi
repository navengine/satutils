"""

Atmosphere
==========

GNSS atmospheric corrections.
"""

from __future__ import annotations
import typing

__all__ = ["IonoModel", "KlobucharElements", "TropoModel"]

class IonoModel:
    """

    IonoModel
    =========

    Klobuchar based ionospheric delay model

    """

    def CalcIonoDelay(
        self, Tow: float, lat: float, lon: float, az: float, el: float, gamma: float
    ) -> float:
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
    def __init__(self) -> None: ...
    @typing.overload
    def __init__(self, klob: KlobucharElements) -> None: ...

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
    def __init__(self) -> None: ...

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

    def __init__(self) -> None: ...
