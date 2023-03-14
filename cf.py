"""
Calculations
"""

from __init__ import ureg, Q_
from math import log, pi


def calc_q_temp_oh(constants=None, props=None):
    mdot_c = constants["mdot_c"]
    assert mdot_c.units == ureg.kg / ureg.seconds
    temp_oc = constants["temp_oc"]
    temp_ic = constants["temp_ic"]
    cp_c = props["cp_c"]
    assert cp_c.units == ureg.joules / (ureg.kg * ureg.degK)

    q = (mdot_c * cp_c * (temp_oc - temp_ic)).to('watts')

    mdot_h = constants["mdot_h"]
    assert mdot_h.units == ureg.kg / ureg.seconds
    temp_ih = constants["temp_ih"]
    cp_h = props["cp_h"]
    assert cp_h.units == ureg.joules / (ureg.kg * ureg.degK)

    temp_oh = ((q / (mdot_h * cp_h) - temp_ih) * -1).to('degK')
    return q, temp_oh


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

    h_c = ((nu_d * k_c)/inner_diam).to(ureg.watts / (ureg.meters**2 * ureg.degK))
    return h_c


def calc_re_d_max(constants=None, props=None, tube_diameter_outer=None):
    rho_h = props["rho_h"]
    mu_h = props["mu_h"]
    nu_h = mu_h / rho_h
    velocity = constants['air_vel']
    s_t = constants["s_t"]
    s_d = constants["s_d"]

    # Step 1 - calculate max velocity of air through tube bank
    vel_max = s_t * velocity / (2 * (s_d - tube_diameter_outer))

    # Step 2 - calculate reynolds number using vmax
    re_d_max = vel_max * tube_diameter_outer / nu_h

    return re_d_max


def calc_nu_hh_20_plus_NL(constants=None, props=None, re_d_max=None, tube_diameter_outer=None):
    # Use different c1 coefficients based on Re_D_max values
    # Assumed staggered array

    s_t = constants['s_t']
    s_l = constants['s_l']
    pr_h = props['pr_h']
    pr_hs = props['pr_hs']
    k_h = props['k_h']

    if 10 <= re_d_max <= 10**2:
        c1 = 0.9
        m = 0.4
        nu_d = c1 * re_d_max ** m * pr_h ** 0.36 * (pr_h / pr_hs) ** 0.25
    elif 10**2 < re_d_max <= 10**3:
        c1 = 0.683
        m = 0.466
        nu_d = c1 * re_d_max**m * pr_h**(1/3)
    elif 10**3 < re_d_max <= 2 * 10**5 and s_t/s_l < 2:
        c1 = 0.35 * (s_t/s_l)**(1/5)
        m = 0.6
        nu_d = c1 * re_d_max ** m * pr_h ** 0.36 * (pr_h / pr_hs) ** 0.25
    elif 10**3 < re_d_max <= 2 * 10**5 and 2 < s_t/s_l:
        c1 = 0.4
        m = 0.6
        nu_d = c1 * re_d_max ** m * pr_h ** 0.36 * (pr_h / pr_hs) ** 0.25
    elif 2 * 10**5 < re_d_max < 2 * 10**6:
        c1 = 0.022
        m = 0.84
        nu_d = c1 * re_d_max ** m * pr_h ** 0.36 * (pr_h / pr_hs) ** 0.25
    else:
        Exception("Re_D_max is out of range")

    h_h = ((nu_d * k_h) / tube_diameter_outer).to(ureg.watts / (ureg.meters**2 * ureg.degK))

    return nu_d, h_h


def calc_temp_tube(tube_diameter_outer=None, tube_thickness=None, h_h=None, temp_film_c=None,
                   temp_film_h=None, k_copper_guess=None):
    # Compare to estimate from above function for Pr_s
    # Use guess temp value to eval k_copper

    inner_diam = tube_diameter_outer - 2 * tube_thickness
    r1 = 0.5 * inner_diam.to(ureg.meters)
    r2 = 0.5 * tube_diameter_outer.to(ureg.meters)
    h_h_value = Q_(h_h, ureg.watts / (ureg.meters**2 * ureg.degK))

    temp = (k_copper_guess * temp_film_c + h_h_value * r2 * temp_film_h * log(r2/r1)) /\
                (k_copper_guess + h_h_value * r2 * log(r2/r1))
    temp_tube = temp.to('degK')

    return temp_tube


def calc_hh_less_20_NL(props=None, nu_20_plus_NL=None, N_L=None, tube_diameter_outer=None):
    k_h = props['k_h']

    c2_values = {
        '1': 0.64,
        '2': 0.76,
        '3': 0.84,
        '4': 0.89,
        '5': 0.92,
        '6': 0.92,
        '7': 0.95,
        '8': 0.95,
        '9': 0.95,
        '10': 0.97,
        '11': 0.97,
        '12': 0.97,
        '13': 0.98,
        '14': 0.98,
        '15': 0.98,
        '16': 0.99,
        '17': 0.99,
        '18': 0.99,
        '19': 0.99
    }
    nu_d = nu_20_plus_NL * c2_values[str(N_L)]
    h_h = (nu_d * k_h / tube_diameter_outer).to(ureg.watts / (ureg.meters**2 * ureg.degK))
    return h_h


# TODO: Remember to include number of tubes in area calculations
def calc_length(copper_props=None, UA=None, h_c=None, h_h=None, tube_diameter_outer=None, tube_thickness=None):
    # Uses overall HT coefficient equations
    inner_diam = tube_diameter_outer - 2 * tube_thickness

    k_copper = copper_props["k_copper"]
    r2 = tube_diameter_outer.to(ureg.meters) * 0.5
    r1 = inner_diam.to(ureg.meters) * 0.5
    h_c = Q_(h_c, ureg.watts / (ureg.meters**2 * ureg.degK))
    h_h = Q_(h_h, ureg.watts / (ureg.meters**2 * ureg.degK))


    l = UA * (1/(h_c * 2 * pi * r1) + log(r2/r1)/(2 * pi * k_copper) + 1/(h_h * 2 * pi * r2))
    length = l.to(ureg.meters)
    return length
