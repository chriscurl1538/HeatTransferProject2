"""
Project 2
"""

from __init__ import ureg, Q_
import cf
import numpy as np


def main():
    # Define Constants
    constants = {
        "temp_amb": 273 * ureg.degK,
        "temp_ic": 293 * ureg.degK,
        "temp_oc": 313 * ureg.degK,
        "temp_ih": 747 * ureg.degK,
        "air_vel": 3000 * ureg.meters/ureg.seconds,   # TODO: Modify if necessary
        "mdot_h": 1736.4 * (ureg.kg / ureg.seconds),
        "mdot_c": (1700 * ureg.kg / ureg.seconds),  # TODO: Modify if necessary
        "s_d": 50 * ureg.mm,
        "s_l": 40 * ureg.mm,
        "s_t": 40 * ureg.mm
                 }

    # Define Parametric variables
    tube_diameter_outer = [Q_(5, 'mm'), Q_(10, 'mm'), Q_(15, 'mm')]
    tube_thickness = [Q_(0.5, 'mm'), Q_(0.75, 'mm'), Q_(1, 'mm')]

    # Properties
    air_props_temp_inlet_r440a = {
        # Air properties @ inlet temp, 747K
        "rho_h": 0.4643 * (ureg.kg / ureg.meters ** 3),
        "cp_h": 1087 * (ureg.joules / (ureg.kg * ureg.degK)),
        "mu_h": 354.6 * 10 ** (-7) * (ureg.pascals * ureg.seconds),
        "nu_h": 76.37 * 10 ** (-6) * (ureg.meters ** 2 / ureg.seconds),
        "pr_h": 0.702,

        # Fluid (liquid) 1 => R440A @ 303K
        "temp_sat": 226.5 * ureg.degK,
        "rho_c": 897.6 * (ureg.kg / ureg.meters ** 3),
        "cp_c": 1801.1 * (ureg.joules / (ureg.kg * ureg.degK)),
        "mu_c": 1.6173 * 10 ** (-4) * (ureg.pascals * ureg.seconds),
        "nu_c": 1.8018 * 10 ** (-7) * (ureg.meters ** 2 / ureg.seconds),
        "pr_c": 2.986,
        "k_c": 0.0975 * (ureg.watts / (ureg.meters * ureg.degK))
    }

    air_props_temp_inlet_r410a = {
        # Air properties @ inlet temp, 747K
        "rho_h": 0.4643 * (ureg.kg / ureg.meters ** 3),
        "cp_h": 1087 * (ureg.joules / (ureg.kg * ureg.degK)),
        "mu_h": 354.6 * 10 ** (-7) * (ureg.pascals * ureg.seconds),
        "nu_h": 76.37 * 10 ** (-6) * (ureg.meters ** 2 / ureg.seconds),
        "k_h": 53.65 * 10 ** (-3) * (ureg.watts / (ureg.meters * ureg.degK)),
        "pr_h": 0.702,

        # Fluid (liquid) 2 => R410A @ 303K
        "temp_sat": 221.8 * ureg.degK,
        "rho_c": 1058 * (ureg.kg / ureg.meters ** 3),
        "cp_c": 1708.1 * (ureg.joules / (ureg.kg * ureg.degK)),
        "mu_c": 1.1798 * 10 ** (-4) * (ureg.pascals * ureg.seconds),
        "nu_c": 1.1144 * 10 ** (-7) * (ureg.meters ** 2 / ureg.seconds),
        "pr_c": 2.264,
        "k_c": 0.089 * (ureg.watts / (ureg.meters * ureg.degK))
    }

    """
    Step 1: Find temp_film using air props at inlet temp
    """

    temp_oh_inlet_props_r440a = cf.calc_q_temp_oh(constants=constants, props=air_props_temp_inlet_r440a)[1]
    temp_oh_inlet_props_r410a = cf.calc_q_temp_oh(constants=constants, props=air_props_temp_inlet_r410a)[1]

    # Film temps
    temp_film_h = 730.5 * ureg.degK
    temp_film_c_r410a = 0.5 * (constants['temp_ic'] + constants['temp_oc'])
    temp_film_c_r440a = 0.5 * (constants['temp_ic'] + constants['temp_oc'])

    # Properties updated (temp_film = 730.5K)
    r440a_props_temp_film = {
        # Air properties @ film temp
        "rho_h": 0.4809 * (ureg.kg / ureg.meters ** 3),
        "cp_h": 1081 * (ureg.joules / (ureg.kg * ureg.degK)),
        "mu_h": 346.7 * 10 ** (-7) * (ureg.pascals * ureg.seconds),
        "nu_h": 72.24 * 10 ** (-6) * (ureg.meters ** 2 / ureg.seconds),
        "k_h": 53.65 * 10**(-3) * (ureg.watts/(ureg.meters * ureg.degK)),
        "pr_h": 0.699,
        "pr_hs": 0.707,     # Note that Pr_s is evaluated at 305K

        # Fluid (liquid) 1 => R440A @ 303K
        "temp_sat": 226.5 * ureg.degK,
        "rho_c": 897.6 * (ureg.kg / ureg.meters ** 3),
        "cp_c": 1801.1 * (ureg.joules / (ureg.kg * ureg.degK)),
        "mu_c": 1.6173 * 10 ** (-4) * (ureg.pascals * ureg.seconds),
        "nu_c": 1.8018 * 10 ** (-7) * (ureg.meters ** 2 / ureg.seconds),
        "pr_c": 2.986,
        "k_c": 0.0975 * (ureg.watts / (ureg.meters * ureg.degK))
    }

    # Properties updated (temp_film = 731.6K)
    r410a_props_temp_film = {
        # Air properties @ film temp
        "rho_h": 0.4809 * (ureg.kg / ureg.meters ** 3),
        "cp_h": 1081 * (ureg.joules / (ureg.kg * ureg.degK)),
        "mu_h": 346.7 * 10 ** (-7) * (ureg.pascals * ureg.seconds),
        "nu_h": 72.24 * 10 ** (-6) * (ureg.meters ** 2 / ureg.seconds),
        "k_h": 53.65 * 10 ** (-3) * (ureg.watts / (ureg.meters * ureg.degK)),
        "pr_h": 0.699,
        "pr_hs": 0.707,  # Note that Pr_s is evaluated at 730K

        # Fluid (liquid) 2 => R410A @ 303K
        "temp_sat": 221.8 * ureg.degK,
        "rho_c": 1058 * (ureg.kg / ureg.meters ** 3),
        "cp_c": 1708.1 * (ureg.joules / (ureg.kg * ureg.degK)),
        "mu_c": 1.1798 * 10 ** (-4) * (ureg.pascals * ureg.seconds),
        "nu_c": 1.1144 * 10 ** (-7) * (ureg.meters ** 2 / ureg.seconds),
        "pr_c": 2.264,
        "k_c": 0.089 * (ureg.watts / (ureg.meters * ureg.degK))
    }

    """
    Step 2: Non-parametric calcs - Find UA value
    """

    # R440A
    q_r440a, temp_oh_r440a = cf.calc_q_temp_oh(constants=constants, props=r440a_props_temp_film)
    c_r_r440a, c_min_r440a, c_max_r440a = cf.calc_cr_cmin_cmax(constants=constants, props=r440a_props_temp_film)
    eta_r440a = cf.calc_eta(constants=constants, c_min=c_min_r440a, q=q_r440a)
    UA_r440a = cf.calc_ntu(eta=eta_r440a, c_r=c_r_r440a, c_min=c_min_r440a)

    # R410A
    q_r410a, temp_oh_r410a = cf.calc_q_temp_oh(constants=constants, props=r410a_props_temp_film)
    c_r_r410a, c_min_r410a, c_max_r410a = cf.calc_cr_cmin_cmax(constants=constants, props=r410a_props_temp_film)
    eta_r410a = cf.calc_eta(constants=constants, c_min=c_min_r410a, q=q_r410a)
    UA_r410a = cf.calc_ntu(eta=eta_r410a, c_r=c_r_r410a, c_min=c_min_r410a)

    """
    Step 3: Parametric calcs - find shell-side and tube-side convection coefficients
    """

    dimensions_mm = np.full(shape=[3, 3], dtype=list, fill_value=0)

    h_c_r440a = np.full(shape=[3, 3], fill_value=0)
    re_d_max_r440a = np.full(shape=[3, 3], fill_value=0)
    nu_20_plus_tubes_r440a = np.full(shape=[3, 3], fill_value=0)
    h_h_20_plus_tubes_r440a = np.full(shape=[3, 3], fill_value=0)
    h_h_less_20_tubes_r440a = np.full(shape=[3, 3], fill_value=0)

    h_c_r410a = np.full(shape=[3, 3], fill_value=0)
    re_d_max_r410a = np.full(shape=[3, 3], fill_value=0)
    nu_20_plus_tubes_r410a = np.full(shape=[3, 3], fill_value=0)
    h_h_20_plus_tubes_r410a = np.full(shape=[3, 3], fill_value=0)
    h_h_less_20_tubes_r410a = np.full(shape=[3, 3], fill_value=0)

    temp_copper_tube_20_plus_NL_r410a = np.full(shape=[3, 3], fill_value=0)
    temp_copper_tube_less_20_NL_r410a = np.full(shape=[3, 3], fill_value=0)
    temp_copper_tube_20_plus_NL_r440a = np.full(shape=[3, 3], fill_value=0)
    temp_copper_tube_less_20_NL_r440a = np.full(shape=[3, 3], fill_value=0)

    length_total_20_plus_tubes_r410a = np.full(shape=[3, 3], fill_value=0)
    length_total_less_20_tubes_r410a = np.full(shape=[3, 3], fill_value=0)
    length_total_20_plus_tubes_r440a = np.full(shape=[3, 3], fill_value=0)
    length_total_less_20_tubes_r440a = np.full(shape=[3, 3], fill_value=0)

    for i, outer_diam in enumerate(tube_diameter_outer):
        for j, thickness in enumerate(tube_thickness):
            st = [outer_diam.magnitude, thickness.magnitude]
            dimensions_mm[i, j] = st

            # R440A
            h_c_r440a[i, j] = cf.calc_hc(constants=constants, props=r440a_props_temp_film, tube_diameter_outer=outer_diam,
                                   tube_thickness=thickness).magnitude
            re_d_max_r440a[i, j] = cf.calc_re_d_max(constants=constants, props=r440a_props_temp_film, tube_diameter_outer=outer_diam)

            nu_20_plus_tubes_r440a[i, j] = cf.calc_nu_hh_20_plus_NL(constants=constants, props=r440a_props_temp_film,
                                                                                       re_d_max=re_d_max_r440a[i, j],
                                                                                       tube_diameter_outer=outer_diam)[0]
            h_h_20_plus_tubes_r440a[i, j] = cf.calc_nu_hh_20_plus_NL(constants=constants, props=r440a_props_temp_film,
                                                                                       re_d_max=re_d_max_r440a[i, j],
                                                                                       tube_diameter_outer=outer_diam)[1].magnitude
            h_h_less_20_tubes_r440a[i, j] = cf.calc_hh_less_20_NL(props=r440a_props_temp_film, nu_20_plus_NL=nu_20_plus_tubes_r440a[i, j],
                                                            tube_diameter_outer=outer_diam, N_L=15).magnitude   # TODO: N_L is placeholder

            # R410A
            h_c_r410a[i, j] = cf.calc_hc(constants=constants, props=r410a_props_temp_film, tube_diameter_outer=outer_diam,
                                   tube_thickness=thickness).magnitude

            re_d_max_r410a[i, j] = cf.calc_re_d_max(constants=constants, props=r410a_props_temp_film, tube_diameter_outer=outer_diam)

            nu_20_plus_tubes_r410a[i, j] = cf.calc_nu_hh_20_plus_NL(constants=constants, props=r410a_props_temp_film,
                                                                    re_d_max=re_d_max_r410a[i, j],
                                                                    tube_diameter_outer=outer_diam)[0]
            h_h_20_plus_tubes_r410a[i, j] = cf.calc_nu_hh_20_plus_NL(constants=constants, props=r410a_props_temp_film,
                                                                     re_d_max=re_d_max_r410a[i, j],
                                                                     tube_diameter_outer=outer_diam)[1].magnitude

            h_h_less_20_tubes_r410a[i, j] = cf.calc_hh_less_20_NL(props=r410a_props_temp_film, nu_20_plus_NL=nu_20_plus_tubes_r410a[i, j],
                                                            tube_diameter_outer=outer_diam, N_L=15).magnitude  # TODO: N_L is placeholder

            # Calc T_s of tube for copper props
            temp_copper_guess = 305 * ureg.degK
            k_copper_guess = 401 * (ureg.watts / (ureg.meters * ureg.degK))

            temp_copper_tube_20_plus_NL_r410a[i, j] = cf.calc_temp_tube(tube_diameter_outer=outer_diam, tube_thickness=thickness,
                                                                  h_h=h_h_20_plus_tubes_r410a[i, j], temp_film_c=temp_film_c_r410a,
                                                                  temp_film_h=temp_film_h, k_copper_guess=k_copper_guess).magnitude
            temp_copper_tube_less_20_NL_r410a[i, j] = cf.calc_temp_tube(tube_diameter_outer=outer_diam, tube_thickness=thickness,
                                                                  h_h=h_h_less_20_tubes_r410a[i, j], temp_film_c=temp_film_c_r410a,
                                                                  temp_film_h=temp_film_h, k_copper_guess=k_copper_guess).magnitude

            temp_copper_tube_20_plus_NL_r440a[i, j] = cf.calc_temp_tube(tube_diameter_outer=outer_diam, tube_thickness=thickness,
                                                                  h_h=h_h_20_plus_tubes_r440a[i, j], temp_film_c=temp_film_c_r440a,
                                                                  temp_film_h=temp_film_h, k_copper_guess=k_copper_guess).magnitude
            temp_copper_tube_less_20_NL_r440a[i, j] = cf.calc_temp_tube(tube_diameter_outer=outer_diam, tube_thickness=thickness,
                                                                  h_h=h_h_less_20_tubes_r440a[i, j], temp_film_c=temp_film_c_r440a,
                                                                  temp_film_h=temp_film_h, k_copper_guess=k_copper_guess).magnitude

            copper_props_305k = {
                "rho_copper": 8933 * ureg.kg / ureg.meters**3,
                "cp_copper": 385 * ureg.joules / (ureg.kg * ureg.degK),
                "k_copper": 401 * (ureg.watts / (ureg.meters * ureg.degK))
            }

            # Calc HX length
            length_total_20_plus_tubes_r410a[i, j] = cf.calc_length(copper_props=copper_props_305k, UA=UA_r410a, h_c=h_c_r410a[i, j],
                                                              h_h=h_h_20_plus_tubes_r410a[i, j], tube_diameter_outer=outer_diam,
                                                              tube_thickness=thickness).magnitude
            length_total_less_20_tubes_r410a[i, j] = cf.calc_length(copper_props=copper_props_305k, UA=UA_r410a, h_c=h_c_r410a[i, j],
                                                              h_h=h_h_less_20_tubes_r410a[i, j], tube_diameter_outer=outer_diam,
                                                              tube_thickness=thickness).magnitude
            length_total_20_plus_tubes_r440a[i, j] = cf.calc_length(copper_props=copper_props_305k, UA=UA_r440a, h_c=h_c_r440a[i, j],
                                                              h_h=h_h_20_plus_tubes_r440a[i, j], tube_diameter_outer=outer_diam,
                                                              tube_thickness=thickness).magnitude
            length_total_less_20_tubes_r440a[i, j] = cf.calc_length(copper_props=copper_props_305k, UA=UA_r440a, h_c=h_c_r440a[i, j],
                                                              h_h=h_h_less_20_tubes_r440a[i, j], tube_diameter_outer=outer_diam,
                                                              tube_thickness=thickness).magnitude

    print("Outer diameter, Thickness (mm)")
    print("...")
    print(dimensions_mm)

    print("R440A - Convection Coefficient (W/m2K)")
    print("...")
    print(h_c_r440a)

    print("Exhaust - Convection Coefficient, 15 tubes (W/m2K)")
    print("...")
    print(h_h_less_20_tubes_r440a)

    print("Exhaust - Convection Coefficient, >20 tubes (W/m2K)")
    print("...")
    print(h_h_20_plus_tubes_r440a)

    print("R440A - Copper Tube Temperature (K)")
    print("...")
    print(temp_copper_tube_less_20_NL_r440a)

    print("R440A - HX length (Nt = 15)")
    print("...")
    print(length_total_less_20_tubes_r440a)

    print("R440A - HX length (Nt > 20)")
    print("...")
    print(length_total_20_plus_tubes_r440a)

    print("R410A - Convection Coefficient (W/m2K)")
    print("...")
    print(h_c_r410a)

    print("Exhaust - Convection Coefficient, 15 tubes (W/m2K)")
    print("...")
    print(h_h_less_20_tubes_r410a)

    print("Exhaust - Convection Coefficient, >20 tubes (W/m2K)")
    print("...")
    print(h_h_20_plus_tubes_r410a)

    print("R410A - Copper Tube Temperature (K)")
    print("...")
    print(temp_copper_tube_less_20_NL_r410a)

    print("R410A - HX length (Nt = 15)")
    print("...")
    print(length_total_less_20_tubes_r410a)

    print("R410A - HX length (Nt > 20)")
    print("...")
    print(length_total_20_plus_tubes_r410a)


if __name__ == "__main__":
    main()
