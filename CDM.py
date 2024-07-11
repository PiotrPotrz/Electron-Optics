import numpy as np
import matplotlib.pyplot as plt
import scipy
import boundary_conditions
import paraxial_ray_equation
import time

epsilon_zero = 8.854187817 * 10**(-12)


def k_function(zi,zj,Ri,Rj): # tez tau
    return np.sqrt((4*Ri*Rj/((Ri+Rj)**2 + (zj-zi)**2)))

def K_eliptic_integral(k:float):
    # Aproksymacja dla małych k
    eta = 1 - k*k
    return 1.386294361 + 0.097932891*eta + 0.054544409*eta**2 + 0.032024666 * eta**3 - (0.5 + 0.124750742 * eta +
                                                                                        0.060118519*eta*2 + 0.010944912
                                                                                        * eta**3)*np.log(eta)

def A_function(si,Ri,Rj,zi,zj):
    return si*k_function(zi,zj,Ri,Rj) * K_eliptic_integral(k_function(zi,zj,Ri,Rj))/(4*np.pi**2*epsilon_zero*np.sqrt(Ri*Rj))


def G_cylindrical(tau,si,ri,rj,zi,zj):
    """
    Function used to create G (A) matrix (3.371) (Electron and Ion Optics).
    Used while evaluating density distribution and axial potential distribution.
    :param tau: Tau/k function given by Eq (3.372)
    :param si: The surface of a cylindrical strip
    :param ri: Radial position of the element affecting the evaluated element
    :param rj: Radial position of evaluated element
    :param zi: Axial position of the element affecting the evaluated element
    :param zj: Axial position of evaluated element
    :return: The matrix element
    """
    return (scipy.special.ellipk(tau**2)*si)/(2*np.pi**2*epsilon_zero*np.sqrt((rj+ri)**2 + (zj-zi)**2)) # zmieniam na tau**2 - wynika to z różnicy definicji funkcji w numpy i EaIO


def cylindrical_strip_surface(r,l):
    return 2*np.pi*r*np.abs(l)

def b_surface(R,r,l):
    return 2*np.pi*r*np.abs(l) + np.pi * (R ** 2 - r ** 2)


def CDM(z: np.array, boundary_conditions: tuple):
    start_time = time.time()

    z_total = boundary_conditions[0]
    r_total = boundary_conditions[1]
    v_boundary = boundary_conditions[2]
    dz_length = boundary_conditions[3]

    # zamieniamy na macierze
    v_boundary = np.matrix(v_boundary)
    density_total = np.matrix(np.zeros(len(z_total)))

    g_matrix = np.zeros((v_boundary.shape[1], v_boundary.shape[1]))

    start_matrix_G_TIME = time.time()

    g1timez = 0

    for j in range(v_boundary.shape[1]):
        for i in range(v_boundary.shape[1]):
            if i != j:
                g_matrix[i][j] = G_cylindrical(k_function(z_total[i], z_total[j], r_total[i], r_total[j]) #z_total[1] - z_total[0]
                                               , cylindrical_strip_surface(r_total[i],dz_length[int(i//(len(z_total)/len(dz_length)))]), # dodajem dz length zamiast z_total[1] - z_total[0] dz_length[int(i//(len(z_total)/len(dz_length)))]
                                               r_total[i], r_total[j], z_total[i], z_total[j])
                g1timez += 1
                print(
                    f"calculating G matrix for density calculation : {(((j - 1) * v_boundary.shape[1] + i) / ((v_boundary.shape[1] - 1) * v_boundary.shape[1])) * 100}%")

    density_total = v_boundary @ np.linalg.inv(g_matrix)

    end_matrix_G_TIME = time.time()

    G_time = end_matrix_G_TIME - start_matrix_G_TIME

    print(f"G {g_matrix.shape}, V_bound = {v_boundary.shape}, density total = {density_total.shape}")

    dens_list = np.ndarray.flatten(np.array(density_total))

    # Obliczanie potencjału na osi

    # obliczanie drugiej macierzy G - G_axial

    G_axial = np.zeros((v_boundary.shape[1], len(z)))
    print(G_axial.shape)

    start_V_time = time.time()

    g2timez = 0

    for j in range(len(z)):
        for i in range(v_boundary.shape[1]):
            if z[j] != z_total[i]:
                G_axial[i][j] = G_cylindrical(
                    k_function(z_total[i], z[j], r_total[i], 0),
                    cylindrical_strip_surface(r_total[i], dz_length[int(i//(len(z_total)/len(dz_length)))]),
                    r_total[i],
                    0,
                    z_total[i],
                    z[j]
                )
                print(
                    f"calculating G matrix for potential calculation : {(((j) * G_axial.shape[0] + i) / (G_axial.shape[0] * G_axial.shape[1])) * 100}%")
                g2timez += 1

    v_axial = density_total @ G_axial

    end_V_time = time.time()

    V_time = end_V_time - start_V_time

    v_axial_list = np.ndarray.flatten(np.array(v_axial))

    """
    v_axial_list = v_axial_list[1:]
    v_axial_list = v_axial_list[:-1]

    z = z[1:]
    z = z[:-1]
    """

    end_time = time.time()

    full_time = end_time - start_time

    print(f"Obliczanie gęstości = {G_time}s\n"
          f"OBliczanie potencjału = {V_time}s\n"
          f"Całkowity czas = {full_time}s\n")

    return v_axial_list


def CDM_thick(z: np.array, boundary_conditions: tuple):
    start_time = time.time()

    z_total = boundary_conditions[0]
    r_total = boundary_conditions[1]
    v_boundary = boundary_conditions[2]
    dz_length = boundary_conditions[3]
    R = boundary_conditions[4]

    e_z_points = int(len(z_total)//len(dz_length)) # liczba punktów na elektrodzie

    point_list = []

    # tworzę listę indeksów brzegowych

    for i in range(1,len(dz_length)+1):
        point_list.append([(i - 1) * e_z_points, i * e_z_points -1])

    point_list = np.ndarray.flatten(np.array(point_list))


    # zamieniamy na macierze
    v_boundary = np.matrix(v_boundary)
    density_total = np.matrix(np.zeros(len(z_total)))

    g_matrix = np.zeros((v_boundary.shape[1], v_boundary.shape[1]))

    start_matrix_G_TIME = time.time()

    g1timez = 0

    for j in range(v_boundary.shape[1]):
        for i in range(v_boundary.shape[1]):
            if i != j:
                if np.isin(i,point_list) == False:
                    g_matrix[i][j] = G_cylindrical(k_function(z_total[i], z_total[j], r_total[i], r_total[j]) #z_total[1] - z_total[0]
                                               , cylindrical_strip_surface(r_total[i],dz_length[int(i//(len(z_total)/len(dz_length)))]), # dodajem dz length zamiast z_total[1] - z_total[0] dz_length[int(i//(len(z_total)/len(dz_length)))]
                                               r_total[i], r_total[j], z_total[i], z_total[j])
                else:
                    g_matrix[i][j] = G_cylindrical(k_function(z_total[i], z_total[j], r_total[i], r_total[j])
                                                   , b_surface(R,r_total[i], dz_length[int(i // (len(z_total) / len(dz_length)))]),
                                                   r_total[i], r_total[j], z_total[i], z_total[j])
                g1timez += 1
                print(
                    f"calculating G matrix for density calculation : {(((j - 1) * v_boundary.shape[1] + i) / ((v_boundary.shape[1] - 1) * v_boundary.shape[1])) * 100}%")

    density_total = v_boundary @ np.linalg.inv(g_matrix)

    end_matrix_G_TIME = time.time()

    G_time = end_matrix_G_TIME - start_matrix_G_TIME

    print(f"G {g_matrix.shape}, V_bound = {v_boundary.shape}, density total = {density_total.shape}")

    dens_list = np.ndarray.flatten(np.array(density_total))

    # Obliczanie potencjału na osi

    # obliczanie drugiej macierzy G - G_axial

    G_axial = np.zeros((v_boundary.shape[1], len(z)))
    print(G_axial.shape)

    start_V_time = time.time()

    g2timez = 0

    for j in range(len(z)):
        for i in range(v_boundary.shape[1]):
            if z[j] != z_total[i]:
                if np.isin(i,point_list) == False:
                    G_axial[i][j] = G_cylindrical(
                    k_function(z_total[i], z[j], r_total[i], 0),
                    cylindrical_strip_surface(r_total[i], dz_length[int(i//(len(z_total)/len(dz_length)))]),
                    r_total[i],
                    0,
                    z_total[i],
                    z[j]
                    )
                else:
                    G_axial[i][j] = G_cylindrical(
                        k_function(z_total[i], z[j], r_total[i], 0),
                        b_surface(R,r_total[i], dz_length[int(i // (len(z_total) / len(dz_length)))]),

                        r_total[i],
                        0,
                        z_total[i],
                        z[j]
                    )
                print(
                    f"calculating G matrix for potential calculation : {(((j) * G_axial.shape[0] + i) / (G_axial.shape[0] * G_axial.shape[1])) * 100}%")
                g2timez += 1

    v_axial = density_total @ G_axial

    end_V_time = time.time()

    V_time = end_V_time - start_V_time

    v_axial_list = np.ndarray.flatten(np.array(v_axial))

    """
    v_axial_list = v_axial_list[1:]
    v_axial_list = v_axial_list[:-1]

    z = z[1:]
    z = z[:-1]
    """

    end_time = time.time()

    full_time = end_time - start_time

    print(f"Obliczanie gęstości = {G_time}s\n"
          f"OBliczanie potencjału = {V_time}s\n"
          f"Całkowity czas = {full_time}s\n")

    return v_axial_list

