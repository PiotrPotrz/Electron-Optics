import numpy as np
import matplotlib.pyplot as plt


def boundary_conditions(electrodes_length: list, slit_length: list, radii: list, voltages: list, point_number: int, outer_radius = 0.0) -> tuple:
    """
    Funkcja zadająca warunki brzegowe do metody gęstości ładunku (CDM).
    :param electrodes_length: lista długości elektrod
    :param slit_length: lista długości przerw
    :param radii: lista długości promieni
    :param voltages: lista potencjałów na elektrodach
    :param point_number: ilośc punktów - wpływa na szybkość i dokładność obliczeń
    :return: (lista z elektrod, lista promieni elektrod, lista potencjałów)
    """
    geometry = ([item for pair in zip(electrodes_length, slit_length) for item in pair]
                + electrodes_length[len(slit_length):]
                + slit_length[len(electrodes_length):])

    for i in range(1, len(geometry)):
        geometry[i] = geometry[i] + geometry[i - 1]

    L = geometry[-1]

    for k in range(len(geometry)):
        geometry[k] = geometry[k] - L / 2

    geometry = np.insert(geometry, 0, -L / 2)

    z_boundary = []

    z_step_array = []

    for i in range(len(electrodes_length)):
        z_boundary.append(np.linspace(geometry[int(2 * i)], geometry[int(2 * i + 1)], point_number))
        z_step_array.append(z_boundary[i][1] - z_boundary[i][0])

    z_boundary_total = []

    for i in range(len(electrodes_length)):
        z_boundary_total = np.append(z_boundary_total, z_boundary[i])

    r = []

    for i in range(len(electrodes_length)):
        r.append(np.ones(point_number) * radii[i])

    total_radii = []

    for i in range(len(electrodes_length)):
        total_radii = np.append(total_radii, r[i])

    v = []

    for i in range(len(electrodes_length)):
        v.append(np.ones(point_number) * voltages[i])

    boundary_voltage = []

    for i in range(len(electrodes_length)):
        boundary_voltage = np.append(boundary_voltage, v[i])

    return (z_boundary_total, total_radii, boundary_voltage, z_step_array, outer_radius)






def boundary_conditions_inclined_at_ends(electrodes_length: list, slit_length: list, radii: list, voltages: list, point_number: int, outer_radius = 0.0) -> tuple:
    electrodes_length = [0.5 * radii[0] * np.sqrt(3)] + electrodes_length + [0.5 * radii[-1] * np.sqrt(3)]
    slit_length = [0] + slit_length + [0]
    voltages = [voltages[0]] + voltages + [voltages[-1]]

    geometry = ([item for pair in zip(electrodes_length, slit_length) for item in pair]
                + electrodes_length[len(slit_length):]
                + slit_length[len(electrodes_length):])

    for i in range(1, len(geometry)):
        geometry[i] = geometry[i] + geometry[i - 1]

    L = geometry[-1]

    for k in range(len(geometry)):
        geometry[k] = geometry[k] - L / 2

    geometry = np.insert(geometry, 0, -L / 2)

    z_boundary = []

    z_step_array = []

    for i in range(len(electrodes_length)):
        z_boundary.append(np.linspace(geometry[int(2 * i)], geometry[int(2 * i + 1)], point_number))
        z_step_array.append(z_boundary[i][1] - z_boundary[i][0])

    z_boundary_total = []

    for i in range(len(electrodes_length)):
        z_boundary_total = np.append(z_boundary_total, z_boundary[i])

    def first_r(z):
        return radii[0]/(z_boundary[0][-1]-z_boundary[0][0])*(z -z_boundary[0][-1]) + radii[0] + 0.0001 * radii[0]

    def last_r(z):
        return -radii[-1]/(z_boundary[-1][-1]-z_boundary[-1][0])*(z -z_boundary[-1][-1]) + 0.0001 * radii[-1]

    r = []

    for i in range(len(electrodes_length)):
        if i!= 0 and i!=len(electrodes_length)-1:
            r.append(np.ones(point_number) * radii[i-1])
            print(i)
        elif i==0:
            r.append(first_r(z_boundary[0][:]))
        elif i==len(electrodes_length)-1:
            r.append(last_r(z_boundary[-1][:]))


    total_radii = []

    for i in range(len(electrodes_length)):
        total_radii = np.append(total_radii, r[i])

    v = []

    for i in range(len(electrodes_length)):
        v.append(np.ones(point_number) * voltages[i])

    boundary_voltage = []

    for i in range(len(electrodes_length)):
        boundary_voltage = np.append(boundary_voltage, v[i])

    return (z_boundary_total, total_radii, boundary_voltage, z_step_array, outer_radius)
