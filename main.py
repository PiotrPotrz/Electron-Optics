import matplotlib.pyplot as plt
import numpy
import matplotlib.pyplot
import findiff
import os
import sys

import numpy as np

import two_cylinder_lens
import CDM
import boundary_conditions
import kontrola
import paraxial_ray_equation
import potential_fix

voltages = []
electrode_length = []
slit_length = []
radii = []
end = 'k'
i = 0
number_of_electrodes = 0

d = 0.01

number = 0
number_e = 0
directory_path = 'C:\\Users\Piotr\PycharmProjects\CIRCULAR APERTURE\wyniki'



if __name__ == '__main__':
    print('PROGRAM OPTYKA ELEKTRONOWA WERSJA 0.0')
    print('podaj ilość elektrod')
    number_of_electrodes = input()
    if number_of_electrodes != end:
        number_of_electrodes = int(number_of_electrodes)
    else:
        exit()

    print('podaj długości elektrod - długość skalowana względem stałej d = 0.01 m')
    for i in range(number_of_electrodes):
        e = input(f'elektroda {i + 1} =')
        if e != end:
            electrode_length.append(float(e) * d)
        else:
            exit()

    print('podaj długość przerw - długość skalowana względem stałej d = 0.01 m')
    for i in range(number_of_electrodes-1):
        e = input(f'przerwa {i+1} =')
        if e != end:
            slit_length.append(float(e) * d)
        else:
            exit()
    if number_of_electrodes == 1:
        slit_length = [0]

    print('podaj potencjały na elektrodach')
    for i in range(number_of_electrodes):
        e = input(f'potencjał na eletrodzie {i+1} [V]')
        if e!= end:
            try:
                voltages.append(int(e))
            except ValueError:
                voltages.append(float(e))
        else:
            exit()


    print('podaj promienie eletrod - skalowane względem stałej d = 0.01 m')
    for i in range(number_of_electrodes):
        e = input(f'promień elektrody {i+1}: ')
        if e!= end:
            radii.append(float(e) * d)
        else:
            exit()



    print('wybierz dokładność obliczeń')
    print('1: niska (krótki czas obliczeń)')
    print('2: średnia')
    print('3: wysoka (długi czas obliczeń)')
    print('4: spersonalizowana')
    print('5: koniec programu')
    e = int(input('Dokładność:  '))

    match e:
        case 1:
            number = 200
            number_e = 2 * number
        case 2:
            number = 300
            number_e = 2 * number
        case 3:
            number = 400
            number_e = 2 * number
        case 4:
            print('Wybrano ustawienia spersonalizowane')
            print('Domyślny model ustawienia siatki: węzły przestrzeni = n, węzły elektrod = 2n')
            print('Dla niektórych przypadków, takich jak \"poszarpane wyniki\" zwiększenie liczby węzłow elektrod znacznie poprawia wyniki.')
            number = int(input('podaj liczbę węzłów przestrzenii: '))
            number_e = int(input('podaj liczbę węzłów na każdej eletrodzie: '))
        case 5:
            exit()

    print(f'n = {number}')
    print(f'n_e = {number_e}')

    ele_l = sum(electrode_length) + sum(slit_length)

    z = np.linspace(-0.5 * ele_l, 0.5 * ele_l, number)

    print('WYBIERZ MODEL OBLICZENIOWY DO METODY GĘSTOŚCI ŁADUNKU')
    print('1: NORMALNY - dla niektórych przypadków zakłamuje rozkład potencjału na brzegach - w celu dobrych obliczeń,\n'
          'szczególnie dla soczewek cienkich, konieczne przycinanie.')
    print('2: ULEPSZONY - czas obliczeń z jego wykorzystaniem jest znacznie dłuższy, nie potrzebne jest tu zwykle przycinanie potencjału.\n'
          'Model wykorzystuje zakończenie za pomocą cylindrów na osi, w celu poprawy potencjału na osi.')
    print('3: Brak metody gęstości ładunku (CDM) - możliwe tylko dla soczewki złożonej z dwóch takich samych współosiowych cylindrów')
    print('4: oba modele obliczeniowe - WIELEOKROTNIE DŁUŻSZY CZAS OBLICZEŃ')
    print('5: KONIEC')

    e = int(input('WYBRANY MODEL: '))



    match e:
        case 1:
            print('model normalny')
            type1 = 'normal'
            type2 = None
        case 2:
            print('model ulepszony')
            type2 = 'inclined'
            type1 = None
        case 3:
            print('brak CDM')
            type1 = None
            type2 = None
        case 4:
            print('oba modele obliczeniowe')
            type1 = 'normal'
            type2 = 'inclined'
        case 5:
            exit()



# Tworzenie ścieżki
    if type1 == 'normal':

        file_name1 = f'{electrode_length},{slit_length},{radii},{voltages},{number},{number_e},{z[0]},{z[-1]},{type1}.txt'
        file_path1 = os.path.join(directory_path, file_name1)


    if type2 == 'inclined':
        file_name2 = f'{electrode_length},{slit_length},{radii},{voltages},{number},{number_e},{z[0]},{z[-1]},{type2}.txt'
        file_path2 = os.path.join(directory_path, file_name2)


    if type1 == 'normal' and type2 == None:
        print('normal')
        bound1 = boundary_conditions.boundary_conditions(electrode_length, slit_length, radii, voltages, number_e)
        # sprawdzanie czy plik istnieje w folderze
        if os.path.isfile(file_path1) == True:
            avn = np.genfromtxt(file_path1,delimiter=',')
        else:
            avn = CDM.CDM(z,bound1)
            np.savetxt(file_path1, avn, fmt='%.15g', delimiter='\t')
        plt.plot(z, avn,label = 'potencjał osiowy')
        plt.scatter(bound1[0],bound1[2],c='r',s=0.25, label = 'V na elektrodach')
        plt.grid()
        plt.legend()
        plt.title(f'POTENCJAŁY = {voltages} V,\nDŁUGOŚCI ELETROD = {electrode_length} m,\nODSTĘPY = {slit_length} m,\nPROMIENIE = {radii} m')
        plt.show()

        kontrola.kontrola_1(z,avn,radii, bound1, voltages, slit_length)



    elif type1 == None and type2 == 'inclined':
        print('inclined')
        bound2 = boundary_conditions.boundary_conditions_inclined_at_ends(electrode_length,slit_length,radii,voltages,number_e)
        bound1 = boundary_conditions.boundary_conditions(electrode_length, slit_length, radii, voltages, number_e)
        if os.path.isfile(file_path2) == True:
            avi = np.genfromtxt(file_path2, delimiter=',')
        else:
            avi = CDM.CDM(z, bound2)
            np.savetxt(file_path2, avi, fmt='%.15g', delimiter='\t')
        plt.plot(z, avi, label='potencjał osiowy')
        plt.scatter(bound1[0], bound1[2], c='r', s=0.25, label='V na elektrodach')
        plt.grid()
        plt.legend()
        plt.title(
            f'POTENCJAŁY = {voltages} V,\nDŁUGOŚCI ELETROD = {electrode_length} m,\nODSTĘPY = {slit_length} m,\nPROMIENIE = {radii} m')
        plt.show()

        kontrola.kontrola_1(z,avi,radii, bound1, voltages, slit_length)
    elif type1 == 'normal' and type2 == 'inclined':
        print('both')
        bound1 = boundary_conditions.boundary_conditions(electrode_length, slit_length, radii, voltages, number_e)
        bound2 = boundary_conditions.boundary_conditions_inclined_at_ends(electrode_length, slit_length, radii,
                                                                          voltages, number_e)

        if os.path.isfile(file_path1) == True and os.path.isfile(file_path2) == True:
            avn = np.genfromtxt(file_path1, delimiter=',')
            avi = np.genfromtxt(file_path2, delimiter=',')
        elif os.path.isfile(file_path1) == True and os.path.isfile(file_path2) == False:
            avn = np.genfromtxt(file_path1, delimiter=',')
            avi = CDM.CDM(z, bound2)
            np.savetxt(file_path2, avi, fmt='%.15g', delimiter='\t')
        elif os.path.isfile(file_path1) == False and os.path.isfile(file_path2) == True:
            avi = np.genfromtxt(file_path2, delimiter=',')
            avn = CDM.CDM(z, bound1)
            np.savetxt(file_path1, avn, fmt='%.15g', delimiter='\t')
        elif os.path.isfile(file_path1) == False and os.path.isfile(file_path2) == False:
            avi = CDM.CDM(z,bound2)
            avn = CDM.CDM(z, bound1)
            np.savetxt(file_path1, avn, fmt='%.15g', delimiter='\t')
            np.savetxt(file_path2, avi, fmt='%.15g', delimiter='\t')
        plt.plot(z, avi, label='potencjał osiowy ulepszony', c = 'k')
        plt.plot(z, avn, label='potencjał osiowy normalny', c = 'b')
        plt.scatter(bound2[0], bound2[2], c='r', s=0.25, label='V na elektrodach')
        plt.grid()
        plt.legend()
        plt.title(
            f'POTENCJAŁY = {voltages} V,\nDŁUGOŚCI ELETROD = {electrode_length} m,\nODSTĘPY = {slit_length} m,\nPROMIENIE = {radii} m')
        plt.show()

        kontrola.kontrola_2(z,avn=avn,avi=avi,radii=radii)

    elif type1 == None and type2 == None:
        if len(voltages) == 2 and np.all(np.unique(np.array(radii) == radii[0])) == True and np.all(np.unique(np.array(electrode_length) == electrode_length[0])) == True:
            print('dwie eleektrody')
            radii = list(radii)
            print(f'promienie = {radii} m')
            electrode_length = list(electrode_length)
            print(f'długości eletrod = {electrode_length} m')

            ava = two_cylinder_lens.axial_potential_general(voltages[0], voltages[1],radii[0],z,slit_length[0])
            plt.plot(z, ava)
            plt.grid()
            plt.title(f'v2/v1 = {voltages[1]/voltages[0]}, przerwa/średnica = {slit_length[0]/(2*radii[0])}')
            plt.show()
            print(type1)
            print(type2)
            kontrola.kontrola_0(z,ava,radii)

        else:
            print('bląd - soczewka nie jest dwuelektrodowa lub nie ma takich samych promieni')
            print(f'promienie = {radii}')
            print(f'ilość eletrod: {number_of_electrodes}')
            print(f'eletrody = {electrode_length}')
            print(f'slitlength = {slit_length}')
            print(f'v = {voltages}')
            print(f'lenvolt = {len(voltages)}')
            print(np.all(np.unique(radii == radii[0])))