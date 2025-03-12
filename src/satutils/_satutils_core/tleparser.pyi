"""

      TLE Parser
      ==========
      
      Tool to parse tle file into SGP4 elements.
"""
from __future__ import annotations
import satutils._satutils_core.atmosphere
__all__ = ['AssumedDecimalPoint', 'CheckSum', 'TleParser']
def AssumedDecimalPoint(assumed_num: float, in_str: str) -> None:
    """
            AssumedDecimalPoint
            ===================
      
            Places an assumed decimal point at the begining of a number
      
            Parameters
            ----------
    
            assumed_num : double
    
                Resulting floating point number
      
            in_str : string
      
                String containing the number
    """
def CheckSum(line: str) -> int:
    """
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
    """
def TleParser(filename: str) -> dict[str, satutils._satutils_core.atmosphere.Sgp4Elements]:
    """
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
    """
