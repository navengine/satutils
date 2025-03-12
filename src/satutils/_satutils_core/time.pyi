"""

      Time
      ====
      
      GNSS time conversions.
"""
from __future__ import annotations
import datetime
__all__ = ['CheckGpsSecond', 'Date2GpsTime', 'HALF_WEEK', 'S_PER_DAY', 'WEEK']
def CheckGpsSecond(t: float) -> float:
    """
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
    """
def Date2GpsTime(week: int, sec: float, dt: datetime.datetime) -> None:
    """
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
    """
HALF_WEEK: float = 302400.0
S_PER_DAY: float = 86400.0
WEEK: float = 604800.0
