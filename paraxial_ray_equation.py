import numpy as np
import scipy as sci
import matplotlib.pyplot as plt
import findiff as fd


def paraxial_ray_equation(z, u, r0, rp0):
    """
    :param z: axial coordinate
    :param u: axial potential distribution
    :param r0: initial ray position
    :param rp0: initial ray position derivative
    :return: ray position
    """
    # sprawdzam czy u jest gdzies zerowe
    u = np.where(u == 0, np.max(u)*10**(-19),u)

    # obliczam pochodne

    """
    Tak to wygladalo wczesniej.
    du = np.gradient(u)/np.gradient(z)
    d2u = np.gradient(du)/np.gradient(z)
    """
    du = np.gradient(u, z)
    d2u = np.gradient(du, z) # zmieniono z z[1]-z[0]
    def f(z, r, rp, u, du, d2u):
        return -du*rp/(2*u) - d2u*r/(4*u)

    r_calc = np.zeros(len(z))
    rp_calc = np.zeros(len(z))
    h = z[1] - z[0]

    r_calc[0] = r0
    rp_calc[0] = rp0

    for i in range(len(z) - 1):
        k1 = f(z[i], r_calc[i], rp_calc[i],u[i],du[i],d2u[i]) * h * h / 2
        k2 = f(z[i] + h / 2, r_calc[i] + h * rp_calc[i] / 2 + k1 / 4, rp_calc[i] + k1 / h,u[i],du[i],d2u[i]) * h * h / 2
        k3 = f(z[i] + h / 2, r_calc[i] + h * rp_calc[i] / 2 + k1 / 4, rp_calc[i] + k2 / h,u[i],du[i],d2u[i]) * h * h / 2
        k4 = f(z[i] + h, r_calc[i] + h * rp_calc[i] + k3, rp_calc[i] + 2 * k3 / h,u[i],du[i],d2u[i]) * h * h / 2
        r_calc[i + 1] = r_calc[i] + h * rp_calc[i] + (k1 + k2 + k3) / 3
        rp_calc[i + 1] = rp_calc[i] + (k1 + 2 * (k2 + k3) + k4) / (3 * h)
        if (k2 - k3) / (k1 - k2) > 0.03:
            print(f'ratio:{(k2-k3)/(k1-k2)}')


    return r_calc


def add_ray_after(z,r, z_max):
    """
    Function adds straight line showing ray after leaving the lens.
    :param z: space parameter
    :param r: ray path
    :return: path of ray after leaving lens
    """
    a = (r[len(r) - 1] - r[len(r) - 2]) / (z[len(r) - 1] - z[len(r) - 2])

    b = r[len(r) - 1] - a * z[len(r) - 1]

    z_after = np.linspace(z[len(z)-1],z_max,len(z))
    r_after = a*z_after + b

    z_fin = np.append(z,z_after)
    r_fin = np.append(r,r_after)

    return [z_fin, r_fin]

def add_ray_after_reverse(z, r, z_min):
    rr = r
    zr = z[::-1]

    a = (rr[-1]-rr[-2])/(zr[-1]-zr[-2])
    b = rr[-1] - a * zr[-1]

    z_min = np.abs(z_min)

    zr_after = np.linspace(zr[-1],-z_min)
    rr_after = a * zr_after + b

    z_fin = np.append(zr,zr_after)
    r_fin = np.append(rr,rr_after)

    return [z_fin, r_fin]



def paraxial_ray_equation_RK4(z, u, r0, rp0):
    """
    Calculated from RK4
    :param z: axial coordinate
    :param u: axial potential distribution
    :param r0: initial ray position
    :param rp0: initial ray position derivative
    :return: ray position
    """
    # sprawdzam czy u jest gdzies zerowe
    u = np.where(u == 0, np.max(u)*10**(-19),u)

    h = z[1] - z[0]

    # obliczam pochodne

    """du = np.gradient(u, z)
    d2u = np.gradient(du, z)
    """

    d_dz = fd.FinDiff(0,h,1)
    d2_dz2 = fd.FinDiff(0,h,2)

    du = d_dz(u)    # 11.11.2023 r. zmiana na FINDIFF
    d2u = d2_dz2(u) # 11.11.2023 r. zmiana na FINDIFF


    def f(z, r, rp, u, du, d2u):
        return rp

    def g(z, r, rp, u, du, d2u):
        return -du*rp/(2*u) - d2u*r/(4*u)



    r_calc = np.zeros(len(z))
    rp_calc = np.zeros(len(z))


    r_calc[0] = r0
    rp_calc[0] = rp0

    for i in range(len(z) - 1):
        k0 = h * f(z[i],r_calc[i],rp_calc[i],u[i],du[i],d2u[i])
        l0 = h * g(z[i],r_calc[i],rp_calc[i],u[i],du[i],d2u[i])

        k1 = h * f(z[i] + 0.5*h, r_calc[i] + 0.5*k0, rp_calc[i] + 0.5*l0, u[i], du[i], d2u[i])
        l1 = h * g(z[i] + 0.5*h, r_calc[i] + 0.5*k0, rp_calc[i] + 0.5*l0, u[i], du[i], d2u[i])

        k2 = h * f(z[i] + 0.5*h, r_calc[i] + 0.5*k1, rp_calc[i] + 0.5*l1, u[i], du[i], d2u[i])
        l2 = h * g(z[i] + 0.5*h, r_calc[i] + 0.5*k1, rp_calc[i] + 0.5*l1, u[i], du[i], d2u[i])

        k3 = h * f(z[i] + h, r_calc[i] + k2, rp_calc[i] + l2, u[i], du[i], d2u[i])
        l3 = h * g(z[i] + h, r_calc[i] + k2, rp_calc[i] + l2, u[i], du[i], d2u[i])

        r_calc[i + 1] = r_calc[i] + (k0 + 2*k1 + 2*k2 + k3 )/6
        rp_calc[i + 1] = rp_calc[i] + (l0 + 2*l1 + 2*l2 + l3 )/6        #dodaje podzielenie na 6 w odu przypadkach 09.11.2023 r.




    return r_calc


def paraxial_ray_equation_RK4_derivative(z, u, r0, rp0):
    """
    Calculated from RK4
    :param z: axial coordinate
    :param u: axial potential distribution
    :param r0: initial ray position
    :param rp0: initial ray position derivative
    :return: ray position
    """
    # sprawdzam czy u jest gdzies zerowe
    u = np.where(u == 0, np.max(u)*10**(-19),u)

    h = z[1] - z[0]

    # obliczam pochodne

    """du = np.gradient(u, z)
    d2u = np.gradient(du, z)
    """

    d_dz = fd.FinDiff(0,h,1)
    d2_dz2 = fd.FinDiff(0,h,2)

    du = d_dz(u)    # 11.11.2023 r. zmiana na FINDIFF
    d2u = d2_dz2(u) # 11.11.2023 r. zmiana na FINDIFF


    def f(z, r, rp, u, du, d2u):
        return rp

    def g(z, r, rp, u, du, d2u):
        return -du*rp/(2*u) - d2u*r/(4*u)



    r_calc = np.zeros(len(z))
    rp_calc = np.zeros(len(z))


    r_calc[0] = r0
    rp_calc[0] = rp0

    for i in range(len(z) - 1):
        k0 = h * f(z[i],r_calc[i],rp_calc[i],u[i],du[i],d2u[i])
        l0 = h * g(z[i],r_calc[i],rp_calc[i],u[i],du[i],d2u[i])

        k1 = h * f(z[i] + 0.5*h, r_calc[i] + 0.5*k0, rp_calc[i] + 0.5*l0, u[i], du[i], d2u[i])
        l1 = h * g(z[i] + 0.5*h, r_calc[i] + 0.5*k0, rp_calc[i] + 0.5*l0, u[i], du[i], d2u[i])

        k2 = h * f(z[i] + 0.5*h, r_calc[i] + 0.5*k1, rp_calc[i] + 0.5*l1, u[i], du[i], d2u[i])
        l2 = h * g(z[i] + 0.5*h, r_calc[i] + 0.5*k1, rp_calc[i] + 0.5*l1, u[i], du[i], d2u[i])

        k3 = h * f(z[i] + h, r_calc[i] + k2, rp_calc[i] + l2, u[i], du[i], d2u[i])
        l3 = h * g(z[i] + h, r_calc[i] + k2, rp_calc[i] + l2, u[i], du[i], d2u[i])

        r_calc[i + 1] = r_calc[i] + (k0 + 2*k1 + 2*k2 + k3 )/6
        rp_calc[i + 1] = rp_calc[i] + (l0 + 2*l1 + 2*l2 + l3 )/6        #dodaje podzielenie na 6 w odu przypadkach 09.11.2023 r.




    return rp_calc





def characteristic_function(z,u)->np.array:
    """
    Equation (7-3) Electron and Ion Optics
    u0 = 0

    :param z: axial coordinate
    :param u: axial potential
    :return: characteristic function value list
    """
    du = np.gradient(u,z)
    u0 = 0
    return 3*((du/(u-u0))**2)/16


def focal_length(z: np.array, u: np.array):
    """
    Numerically aproximated focal lengths for thin electrostatic lenses.
    Equations (4 - 117) and (4 - 118) - Electron and Ion Optics
    :param z: axial coordinate
    :param u: axial potential
    :return: [Object side focal length, Image side focal length]
    """
    a = z[0]
    b = z[-1]
    u0 = 0 # potential where the speed of the particle is equal to zero - means at the infinity distance from the lens

    one_per_f1 = np.sqrt(np.sqrt((u[-1]-u0)/(u[0]-u0))) * np.trapz(characteristic_function(z,u), x=z)

    one_per_f2 = np.sqrt(np.sqrt((u[0]-u0)/(u[-1]-u0))) * np.trapz(characteristic_function(z,u), x=z)

    return [1/one_per_f1, 1/one_per_f2]


def many_rays_angles(
        z: np.array, u: np.array, r0, number_of_rays,
        max_angle = 0.05, min_angle = - 0.05,
        text = 'PROMIEŃ ELEKTRONOWY', d = None
    ):

    step = (max_angle - min_angle)/number_of_rays

    starting_derivtive_list = np.arange(max_angle,min_angle,-step)

    rays_list = []

    for i in range(len(starting_derivtive_list)):
        rays_list.append(paraxial_ray_equation_RK4(z,u,r0,starting_derivtive_list[i]))

    return rays_list

def focal_point(z :np.array, u :np.array):
    ztdz = np.trapz(z*characteristic_function(z,u), x=z)
    tdz = np.trapz(characteristic_function(z,u), x=z)
    focal_point1 = (ztdz-1)/tdz
    focal_point2 = (1+ztdz)/tdz
    return [focal_point1, focal_point2]

def center_of_field(z,u):
    ztdz = np.trapz(z*characteristic_function(z,u), x=z)
    tdz = np.trapz(characteristic_function(z,u), x=z)
    return ztdz/tdz

def reverse_paraxial_ray_equation_RK4(z_standard,u_standard,r0,rp0):
    """
    Funkcja wykonuje obliczenia promienia osiowego wpadającego z drugiej strony soczewki.
    Funkcja wykorzystuje funkcję paraxial_ray_eauation_RK4 dla odwróconych wartości.
    :param z_standard: Współrzędna osiowa [m]
    :param u_standard:Potencjał na osi [V]
    :param r0: Promień początkowy [m].
    :param rp0: Kąt pod którym wpada promień (początkowa wartość pochodnej).
    :return: Trajektoria promienia wpadającego od końca soczewki.
    """
    #z_reversed = z_standard[::-1]
    u_reversed = u_standard[::-1]
    return paraxial_ray_equation_RK4(z_standard,u_reversed,r0,rp0)

def reverse_paraxial_ray_equation_RK4_2(z_standard,u_standard,r0,rp0):
    """
    Funkcja wykonuje obliczenia promienia osiowego wpadającego z drugiej strony soczewki.
    Funkcja wykorzystuje funkcję paraxial_ray_eauation_RK4 dla odwróconych wartości.
    :param z_standard: Współrzędna osiowa [m]
    :param u_standard:Potencjał na osi [V]
    :param r0: Promień początkowy [m].
    :param rp0: Kąt pod którym wpada promień (początkowa wartość pochodnej).
    :return: Trajektoria promienia wpadającego od końca soczewki.
    """
    u_reversed = u_standard[::-1]
    return paraxial_ray_equation_RK4(z_standard,u_reversed,r0,rp0)

def many_rays_angles_reverse(zr: np.array, ur: np.array, r0, number_of_rays,
        max_angle = 0.05, min_angle = - 0.05,
        text = 'PROMIEŃ ELEKTRONOWY', d = None):

    z = zr
    u = ur

    step = (max_angle - min_angle)/number_of_rays

    starting_derivtive_list = np.arange(max_angle,min_angle,-step)

    rays_list = []

    for i in range(len(starting_derivtive_list)):
        rays_list.append(reverse_paraxial_ray_equation_RK4(z,u,r0,starting_derivtive_list[i]))

    return rays_list