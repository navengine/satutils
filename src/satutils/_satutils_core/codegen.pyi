"""

      Code-Gen
      ========
      
      Constellation gold code generators.
"""
from __future__ import annotations
__all__ = ['CodeGenCA']
def CodeGenCA(ca_code: bool, prn_id: int) -> None:
    """
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
    """
