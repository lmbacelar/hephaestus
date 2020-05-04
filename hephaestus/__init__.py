"""
Python module containing calibration data and lookup functions for standard
thermocouples of types B, C, D, E, G, J, K, N, R, S, T, and a few more
"""

__author__ = "lmbacelar@gmail.com"
__copyright__ = "public domain"
__version__ = "0.0.1"

from . import transfer_functions
from . import nist
from . import astm
from . import omega

# assemble thermocouple list.
# note the order: type G from OMEGA source is replaced by the ASTM one.
thermocouples = dict()
thermocouples.update(omega.thermocouples)
thermocouples.update(nist.thermocouples)
thermocouples.update(astm.thermocouples)
