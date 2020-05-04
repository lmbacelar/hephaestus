"""    
This module contains thermocouple reference functions for types C,D,G.

Thermocouple reference functions are polynomials taken from
    "Tungsten-Rhenium Thermocouples Calibration Equivalents"
    http://www.omega.com/temperature/z/pdf/z202.pdf
"""

__copyright__ = "public domain"

import numpy as np
from .transfer_functions import Thermocouple, Polynomial_Gaussian

_source = "OMEGA Inc. z202.pdf"
_cal = "IPTS-68"

thermocouples = {
    # The coefficients of this polynomial have been converted from
    # a Fahrenheit polynomial with 7 significant figures, and have
    # been stored here with longer precision to avoid further inaccuracy.
    "G": Thermocouple(
        Polynomial_Gaussian(
            table=[
                [
                    0.0,
                    2315.0,
                    np.array(
                        [
                            -2.2222283359680003e-16,
                            2.2112702944358399e-12,
                            -1.0316119658501839e-08,
                            2.1425207201941232e-05,
                            1.2905824431600024e-03,
                            0.0000000000000000e00,
                        ]
                    ),
                    None,
                ]
            ],
            invtable=None,
            calibration=_cal,
            source=_source + ", type G",
        ),
        ttype="Type G",
    ),
    # The coefficients of this polynomial have been converted from
    # a Fahrenheit polynomial with 7 significant figures, and have
    # been stored here with longer precision to avoid further inaccuracy.
    "C": Thermocouple(
        Polynomial_Gaussian(
            table=[
                [
                    0.0,
                    2315.0,
                    np.array(
                        [
                            -4.9446064258560002e-16,
                            3.6006582486412798e-12,
                            -1.0489145155399069e-08,
                            1.2252598548103214e-05,
                            1.3387722982319094e-02,
                            0.0000000000000000e00,
                        ]
                    ),
                    None,
                ]
            ],
            invtable=None,
            calibration=_cal,
            source=_source + ", type C",
        ),
        ttype="Type C",
    ),
    "D": Thermocouple(
        Polynomial_Gaussian(
            table=[
                [
                    0.0,
                    783.0,
                    np.array(
                        [
                            -1.4240735e-15,
                            7.9498033e-12,
                            -1.8464573e-8,
                            2.0592621e-5,
                            9.5685256e-3,
                            0.0,
                        ]
                    ),
                    None,
                ],
                [
                    783.0,
                    2320.0,
                    np.array(
                        [
                            -7.9026726e-16,
                            5.3743821e-12,
                            -1.4935266e-8,
                            1.8666488e-5,
                            9.9109462e-3,
                            0.0,
                        ]
                    ),
                    None,
                ],
            ],
            invtable=None,
            calibration=_cal,
            source=_source + ", type D",
        ),
        ttype="Type D",
    ),
}

del _source
del _cal

# end of module
