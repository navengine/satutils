"""

GPS-LNAV
========

Implementation of GPS L1 C/A navigation message utils.
"""

from __future__ import annotations
import satutils._satutils_core.atmosphere
import satutils._satutils_core.ephemeris
import typing

__all__ = ["GpsLnav"]

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

    def GetEphemerides(self) -> satutils._satutils_core.ephemeris.KeplerElements:
        """
        GetEphemerides
        ==============

        Get the Keplerian ephemeris elements

        Returns
        -------

        eph : KeplerElements

            Keplerian ephemeris elements struct
        """

    def GetKlobuchar(self) -> satutils._satutils_core.atmosphere.KlobucharElements:
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

    def SetEphemerides(self, eph: satutils._satutils_core.ephemeris.KeplerElements) -> None:
        """
        SetEphemerides
        ==============

        Set the Keplerian ephemeris elements

        Parameters
        ----------

        eph : KeplerElements

            Keplerian ephemeris elements struct
        """

    def SetKlobuchar(self, klob: satutils._satutils_core.atmosphere.KlobucharElements) -> None:
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
    def __init__(self) -> None: ...
    @typing.overload
    def __init__(self, eph: satutils._satutils_core.ephemeris.KeplerElements) -> None: ...
    @typing.overload
    def __init__(self, klob: satutils._satutils_core.atmosphere.KlobucharElements) -> None: ...
    @typing.overload
    def __init__(
        self,
        eph: satutils._satutils_core.ephemeris.KeplerElements,
        klob: satutils._satutils_core.atmosphere.KlobucharElements,
    ) -> None: ...
