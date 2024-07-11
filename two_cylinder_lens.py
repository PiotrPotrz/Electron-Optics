import numpy as np
import matplotlib.pyplot as plt
from scipy import integrate

import paraxial_ray_equation

#print("TO JEST PRZYKŁAD DLA SOCZEWKI Z DWOMA CYLINDRAMI")

omega = 1.318


def axial_potential_narrow_gap(v1, v2, radius, z):
    """
    Function calculates axial potential distribution for double cylinder lens with infinitely
    narrow distance between electrodes s.
    Equation (3 - 132)
    :param v1: potential of the first electrode
    :param v2: potential of the second electrode
    :param radius: radius of the lens
    :param z: the axial coordinate
    :return: axial potential distribution
    """
    return (v1+v2)/2 + (v2-v1)/2 * np.tanh(omega * z / radius)

def axial_potential_general(v1, v2, radius, z, s):
    """
    Function calculates axial potential distribution for double cylinder lens with
    distance between electrodes equal to s.
    Equation (3 - 131)
    :param v1: potential of the first electrode
    :param v2: potential of the second electrode
    :param radius: radius of the lens
    :param s: distance between electrodes
    :return:
    """
    inlog_term = np.cosh(omega * (z + 0.5 * s)/radius) / np.cosh(omega * (z - 0.5 * s)/radius)
    return (v1 + v2)/2 + radius *(v2 - v1)/(2 * omega * s) * np.log(inlog_term)

def focal_length(v1, v2, radius, z, u):
    """
    Calculating position of the focal lengths of the lens.

    Zaniża wyniki w porównaniu z Hartringiem, Readem dla dużych v2/v1
    :param v1: potential of the first electrode
    :param v2: potential of the second electrode
    :param radius: radius of the lens
    :param z: axial coordinate
    :param u: axial potential distribution
    :return: [object side focal length, image side focal length]
    """
    # Determining u0 parameter - potential on the z axis where the particle velocity is zero
    u0 = u[0]
    # Checking what happens if u0 == 0
    u0 = 0
    # Calculating ratio of voltages
    ratio_v = (v2 - u0)/(v1 - u0)
    #print(f"ratio_v = {ratio_v}, u0 ={u0}")


    # Calculating the first focal length (object side f1)

    Rpf1 = 0.495 * (np.abs(ratio_v)**0.25) * ((ratio_v + 1)/(ratio_v - 1)*np.log(np.abs(ratio_v)) - 2)

    f1 = radius/Rpf1

    # Calculating image side focal length

    Rpf2 = 0.495 * (np.abs(1/ratio_v) ** 0.25) * ((ratio_v + 1) / (ratio_v - 1) * np.log(np.abs(ratio_v)) - 2)
    f2 = radius/Rpf2

    return np.array([f1,f2])


