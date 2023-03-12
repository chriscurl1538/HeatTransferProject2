"""
Calculations
"""

from __init__ import ureg, Q_
from math import log, pi
import numpy as np


def calc_q_temp_oh(constants=None, props=None):
    mdot_c = constants["mdot_c"]
    temp_oc = constants["temp_oc"]
    temp_ic = constants["temp_ic"]
    cp_c = props["cp_c"]

    q = mdot_c * cp_c * (temp_oc - temp_ic)

    mdot_h = constants["mdot_h"]
    T_ih = constants["T_ih"]
    cp_h = props["cp_h"]

    T_oh = (q / (mdot_h * cp_h) - T_ih) * -1
    return q, T_oh


def calc_cr_cmin_cmax(constants=None, props=None):
    mdot_h = constants["mdot_h"]
    mdot_c = constants["mdot_c"]
    cp_h = props["cp_h"]
    cp_c = props["cp_c"]

    c_h = mdot_h * cp_h
    c_c = mdot_c * cp_c

    if c_h < c_c:
        c_r = c_h / c_c
        c_min = c_h
        c_max = c_c
    else:
        c_r = c_c / c_h
        c_min = c_c
        c_max = c_h
    return c_r, c_min, c_max


def calc_eta(constants=None, c_min=None, q=None):
    temp_ih = constants["temp_ih"]
    temp_ic = constants["temp_ic"]

    q_max = c_min * (temp_ih - temp_ic)
    eta = q / q_max
    return eta


def calc_ntu(eta=None, c_r=None, c_min=None):
    F = (eta * c_r - 1)/(eta - 1)
    eta_1 = (F - 1)/(F - c_r)
    E = (2/eta_1 - (1 + c_r)) * (1 + c_r**2)**(-1/2)
    ntu = -1 * (1 + c_r**2)**(-1/2) * log((E-1)/(E+1))
    UA = ntu * c_min
    return UA


"""
Parametric calculations
"""


# TODO: Do parametric calc in main.py
def calc_hc(constants=None, props=None, tube_diameter_outer=None, tube_thickness=None):
    mdot_c = constants["mdot_c"]
    mu_c = props["mu_c"]
    pr_c = props["pr_c"]
    k_c = props["k_c"]

    inner_diam = tube_diameter_outer - 2 * tube_thickness
    re_d = (4 * mdot_c)/(pi * inner_diam * mu_c)

    if re_d > 3500:
        # Flow is turbulent
        nu_d = 0.023 * re_d**(4/5) * pr_c**0.4
    else:
        # Flow is laminar
        nu_d = 4.36     # Uniform heat flux

    h_c = (nu_d * k_c)/inner_diam
    return h_c


def calc_hh():
    return None


# TODO: Remember to include number of tubes in area calculations
def calc_length_parametric():
    return None
