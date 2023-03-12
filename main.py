"""
Project 2
"""

from __init__ import ureg, Q_


def main():
    # Constants
    constants = {
        "temp_amb": 273 * ureg.degK,
        "temp_ic": 293 * ureg.degK,
        "temp_oc": 313 * ureg.degK,
        "temp_ih": (474 + 273) * ureg.degK,
        "mdot_h": 1736.4 * (ureg.meters**3 / ureg.seconds),
        "mdot_c": (1 * ureg.kg / ureg.seconds)  # TODO: Modify if necessary
                 }

    # Parametric variables TODO: select values
    tube_diameter_outer = [5, 10, 15] * ureg.mm
    tube_thickness = [0.5, 0.75, 1] * ureg.mm

    # Properties TODO: Add properties here
    props = {}
