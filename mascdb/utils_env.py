# -----------------------------------------------------------------------------.
# Copyright (c) 2021-2025 MASCDB developers
#
# This program is free software: you can redistribute it and/or modify
# it under the terms of the MIT License.
#
# This program is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.
#
# You should have received a copy of the MIT License
# along with this program.  If not, see <https://opensource.org/license/mit/>.
# -----------------------------------------------------------------------------.
"""MASCDB auxiliary environmental functions."""


import numpy as np


def wet_bulb_t(t, rh):
    """
    Returns Wet bulb temperature estimated from T and RH.

    Parameters
    ----------
    t : float, int, list or numpy.ndarray
        Temperature in degree Celsius
    rh : float, int, list or numpy.ndarray
        Relative humidity in percentage

    Returns
    -------
    array-like
        Wet bulb temperature in °C, same data type as as t/rh

    """
    tt = np.asarray(t)
    rr = np.asarray(rh)

    t_wetbulb = (
        tt * np.arctan(0.152 * (rr + 8.3136) ** (1 / 2))
        + np.arctan(tt + rr)
        - np.arctan(rr - 1.6763)
        + 0.00391838 * (rr) ** (3 / 2) * np.arctan(0.0231 * rr)
        - 4.686
    )

    return t_wetbulb
