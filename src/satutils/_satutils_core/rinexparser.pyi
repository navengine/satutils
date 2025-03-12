"""

      Rinex Parser
      ============
      
      Tool to parse Rinex 3.0 file into ephemeris elements.
"""
from __future__ import annotations
import satutils._satutils_core.atmosphere
__all__ = ['ParseNavBlock', 'ParseTimeBlock', 'RinexParser']
def ParseNavBlock(line: str) -> list[float]:
    """
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
    """
def ParseTimeBlock(week: int, sec: float, line: str) -> None:
    """
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
    """
def RinexParser(filename: str) -> dict[str, tuple[satutils._satutils_core.atmosphere.KlobucharElements, satutils._satutils_core.atmosphere.KeplerElements]]:
    """
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
    """
