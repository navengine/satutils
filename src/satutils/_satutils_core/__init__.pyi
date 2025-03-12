"""

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
      
"""
from __future__ import annotations
from . import atmosphere
from . import codegen
from . import ephemeris
from . import gpslnav
from . import rinexparser
from . import time
from . import tleparser
__all__ = ['GALILEO_E1_CODE_LENGTH', 'GALILEO_E1_DATA_RATE', 'GALILEO_E5A_DATA_RATE', 'GALILEO_E5A_FREQUENCY', 'GALILEO_E5B_DATA_RATE', 'GALILEO_E5B_FREQUENCY', 'GALILEO_E5_CODE_LENGTH', 'GALILEO_E5_FREQUENCY', 'GALILEO_E6_CODE_LENGTH', 'GALILEO_E6_CODE_RATE', 'GALILEO_E6_DATA_RATE', 'GALILEO_E6_FREQUENCY', 'GPS_CA_CODE_LENGTH', 'GPS_CA_CODE_RATE', 'GPS_L1_FREQUENCY', 'GPS_L2CL_CODE_LENGTH', 'GPS_L2CM_CODE_LENGTH', 'GPS_L2_CODE_RATE', 'GPS_L2_FREQUENCY', 'GPS_L5_CODE_LENGTH', 'GPS_L5_CODE_RATE', 'GPS_L5_FREQUENCY', 'GPS_PI', 'LNAV_INV_PREAMBLE_BITS', 'LNAV_PREAMBLE_BITS', 'LNAV_SUBFRAME_SIZE', 'LNAV_WORD_SIZE', 'SGP_A3OVK2', 'SGP_AE', 'SGP_CK2', 'SGP_CK4', 'SGP_QOMS2T', 'SGP_S', 'SGP_XKE', 'SGP_XKMPER', 'TWO_GPS_PI', 'atmosphere', 'codegen', 'ephemeris', 'gpslnav', 'rinexparser', 'time', 'tleparser']
GALILEO_E1_CODE_LENGTH: int = 4092
GALILEO_E1_DATA_RATE: float = 250.0
GALILEO_E5A_DATA_RATE: float = 50.0
GALILEO_E5A_FREQUENCY: float = 1176450000.0
GALILEO_E5B_DATA_RATE: float = 250.0
GALILEO_E5B_FREQUENCY: float = 1207140000.0
GALILEO_E5_CODE_LENGTH: int = 10230
GALILEO_E5_FREQUENCY: float = 1191795000.0
GALILEO_E6_CODE_LENGTH: int = 5115
GALILEO_E6_CODE_RATE: float = 5115000.0
GALILEO_E6_DATA_RATE: float = 1000.0
GALILEO_E6_FREQUENCY: float = 1278750000.0
GPS_CA_CODE_LENGTH: int = 1023
GPS_CA_CODE_RATE: float = 1023000.0
GPS_L1_FREQUENCY: float = 1575420000.0
GPS_L2CL_CODE_LENGTH: int = 767250
GPS_L2CM_CODE_LENGTH: int = 10230
GPS_L2_CODE_RATE: float = 511500.0
GPS_L2_FREQUENCY: float = 1227600000.0
GPS_L5_CODE_LENGTH: int = 10230
GPS_L5_CODE_RATE: float = 10230000.0
GPS_L5_FREQUENCY: float = 1176450000.0
GPS_PI: float = 3.1415926535898
LNAV_INV_PREAMBLE_BITS: int = 116
LNAV_PREAMBLE_BITS: int = 139
LNAV_SUBFRAME_SIZE: int = 300
LNAV_WORD_SIZE: int = 30
SGP_A3OVK2: float = 0.0046780659154137035
SGP_AE: float = 1.0
SGP_CK2: float = 0.00054131345
SGP_CK4: float = 6.0765e-07
SGP_QOMS2T: float = 1.88027916e-09
SGP_S: float = 1.01222928
SGP_XKE: float = 0.07436691613317341
SGP_XKMPER: float = 6378.135
TWO_GPS_PI: float = 6.2831853071796
__version__: str = '1.0.0'
