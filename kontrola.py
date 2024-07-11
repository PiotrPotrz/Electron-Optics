import matplotlib.pyplot as plt
import numpy
import matplotlib.pyplot
import findiff
import os

import numpy as np

import two_cylinder_lens
import CDM
import boundary_conditions
import paraxial_ray_equation
import potential_fix


def kontrola_1(z,av,radii, bound = [], voltages = [], slits = []):
    end = 'k'
    e = None
    d = 0.01

    d2_dz2 = findiff.FinDiff(0, z[1]-z[0], 2)
    while e != end:
        print('\n\nWYBIERZ CO ROBIĆ DALEJ:\n\n')
        print('1: Oblicz trajektrodię promieni elektronowych lecacą z lewej strony')
        print('2: Oblicz trajektrodię promieni elektronowych lecacą z prawej strony')
        print('3: Oblicz wydłużoną trajektorię promieni z lewej strony')
        print('4: Oblicz wydłużoną trajektorię promieni z prawej strony')
        print('5: Wyświetl druga pochodną potencjału')
        print('6: Przytnij potencjał i zastąp stałymi wartościami - pomaga usunąć błędy na brzegach')
        print('k: Zakończ')
        e = input('wybór: ')
        if e=='k':
            exit()
        else:
            e = int(e)

        if e == 1:
            print('trajektoria z lewej')
            print('a: Poglądowa ilość promieni (do wybrania r0, kąt_max, kąt_min, ilość)')
            print('b: Konktretny promień (r0, kąt)')
            print('c: Poglądowa ilość promieni (do wybrania r0, kąt_max, kąt_min, ilość) i rozkład potencjału')
            print('d: Poglądowa ilość promieni (do wybrania r0, kąt_max, kąt_min, ilość) i rozkład potencjału z potencjałem analitycznym')
            choice = input('wybór:')
            if choice == 'a':
                print('a')
                r0 = float(input('r0 [d] = ')) * d
                rp0max = float(input('kąt maksymalny [rad] = '))
                rp0min = float(input('kąt minimalny [rad] = '))
                n = int(input('liczba promieni = '))
                ray_list = paraxial_ray_equation.many_rays_angles(z,av,r0,n,rp0max,rp0min)
                for i in range(len(ray_list)):
                    plt.plot(z/((2*radii[0])),ray_list[i]/(2*radii[0]))
                plt.xlabel('z/D',fontsize=24)
                plt.ylabel('r/D',fontsize=24)
                plt.xticks(fontsize=18)
                plt.yticks(fontsize=18)
                plt.grid()
                plt.show()


            elif choice == 'b':
                print('b')
                r0 = float(input('r0 [d] = ')) * d
                rp0 = float(input('kąt [rad] = '))
                plt.plot(z/(2*radii[0]),paraxial_ray_equation.paraxial_ray_equation_RK4(z,av,r0,rp0)/(2*radii[0]))
                plt.xlabel('z/D', fontsize=24)
                plt.ylabel('r/D', fontsize=24)
                plt.xticks(fontsize=18)
                plt.yticks(fontsize=18)
                plt.grid()
                plt.show()
            elif choice == 'c':
                print('c')
                r0 = float(input('r0 [d] = ')) * d
                rp0max = float(input('kąt maksymalny [rad] = '))
                rp0min = float(input('kąt minimalny  [rad] = '))
                n = int(input('liczba promieni = '))
                fig, ax = plt.subplots(2, 1, figsize=(12, 8))
                ax[0].plot(z / (2 * radii[0]), av, c='blue')
                ax[0].set_title(f'Rozkład potencjału dla {voltages} V')
                ax[0].grid()
                ax[0].set_ylabel('V [V]')
                ax[0].set_xlabel('z/D')
                # ax[0].set_xticks(fontsize=18)
            # ax[0].set_yticks(fontsize=18)

                ax[0].scatter(bound[0] / (2 * radii[0]), bound[2], s=0.5, c='r')
            # ax[0].legend()

                ray_list = paraxial_ray_equation.many_rays_angles(z, av, r0, n, max_angle=rp0max,
                                                              min_angle=rp0min)

                for i in range(len(ray_list)):
                    ax[1].plot(z / (2 * radii[0]), ray_list[i] / (2 * radii[0]))

                ax[1].grid()
                ax[1].set_title('PROMIENIE ELEKTRONOWE')
                ax[1].set_xlabel('z/D')
                ax[1].set_ylabel('r/D')
            # ax[1].set_xticks(fontsize=18)
            # ax[1].set_yticks(fontsize=18)

                plt.show()
            else:
                print('d')
                ava = two_cylinder_lens.axial_potential_general(voltages[0], voltages[1], bound[1][0], z, slits[0])
                r0 = float(input('r0 [d] = ')) * d
                rp0max = float(input('kąt maksymalny [rad] = '))
                rp0min = float(input('kąt minimalny [rad] = '))
                n = int(input('liczba promieni = '))
                fig, ax = plt.subplots(2, 1, figsize=(12, 8))
                ax[0].plot(z / (2 * radii[0]), av, c='blue', label='CDM')
                ax[0].plot(z / (2 * radii[0]), ava, c='green', label='analityczny')
                ax[0].set_title(f'Rozkład potencjału dla {voltages} V')
                ax[0].grid()
                ax[0].set_ylabel('V [V]')
                ax[0].set_xlabel('z/D')
            # ax[0].tick_params(axis='x', labelsize=18)
            # ax[0].tick_params(axis='y', labelsize=18)
                ax[0].legend()
                ax[0].scatter(bound[0] / (2 * radii[0]), bound[2], s=0.5, c='r')

                ray_list = paraxial_ray_equation.many_rays_angles(z, av, r0, n, max_angle=rp0max,
                                                              min_angle=rp0min)

                for i in range(len(ray_list)):
                    ax[1].plot(z / (2 * radii[0]), ray_list[i] / (2 * radii[0]))

                ax[1].grid()
                ax[1].set_title('PROMIENIE ELEKTRONOWE')
                ax[1].set_xlabel('z/D')
                ax[1].set_ylabel('r/D')
            # ax[1].tick_params(axis='x', labelsize=18)
            # ax[1].tick_params(axis='y', labelsize=18)

                plt.show()
        elif e == 2:
            print('trajektroria z prawej')
            print('a: Poglądowa ilość promieni (do wybrania r0, kąt_max, kąt_min, ilość)')
            print('b: Konktretny promień (r0, kąt)')
            choice = input('wybór:')
            if choice == 'a':
                print('a')
                r0 = float(input('r0 [d] = ')) * d
                rp0max = float(input('kąt_max [rad] = '))
                rp0min = float(input('kąt_min [rad] = '))
                n = int(input('liczba promieni = '))
                rray_list = paraxial_ray_equation.many_rays_angles_reverse(z, av, r0, n, rp0max, rp0min)
                for i in range(len(rray_list)):
                    plt.plot(z[::-1]/((2*radii[0])), rray_list[i]/(2*radii[0]))
                plt.xlabel('z/D', fontsize=24)
                plt.ylabel('r/D', fontsize=24)
                plt.xticks(fontsize=18)
                plt.yticks(fontsize=18)
                plt.grid()
                plt.show()
            else:
                print('b')
                r0 = float(input('r0 [d] = ')) * d
                rp0 = float(input('kąt [rad] = '))
                plt.plot(z[::-1] / (2*radii[0]), paraxial_ray_equation.reverse_paraxial_ray_equation_RK4(z, av, r0, rp0)/(2*radii[0]))
                plt.xlabel('z/D', fontsize=24)
                plt.ylabel('r/D', fontsize=24)
                plt.xticks(fontsize=18)
                plt.yticks(fontsize=18)
                plt.grid()
                plt.show()
        elif e == 3:
            print('trajektrodia z lewej wydłużona')
            print('a: Poglądowa ilość promieni (do wybrania r0, kąt_max, kąt_min, ilość)')
            print('b: Konktretny promień (r0, kąt)')
            choice = input('wybór: ')
            if choice == 'a':
                print('a')
                r0 = float(input('r0 [d] = ')) * d
                rp0max = float(input('kąt maksymalny [rad] = '))
                rp0min = float(input('kąt minimalny [rad] = '))
                n = int(input('liczba promieni = '))
                add_len = float(input('wydłużenie [d] =')) * d
                ray_list = paraxial_ray_equation.many_rays_angles(z, av, r0, n, rp0max, rp0min)
                raya_list = []
                for i in range(np.shape(ray_list)[0]):
                    raya_list.append(paraxial_ray_equation.add_ray_after(z, ray_list[i], add_len))
                    plt.plot(raya_list[i][0]/((2*radii[0])), raya_list[i][1]/(2*radii[0]))
                plt.xlabel('z/D', fontsize=24)
                plt.ylabel('r/D', fontsize=24)
                plt.xticks(fontsize=18)
                plt.yticks(fontsize=18)
                plt.grid()
                plt.show()

            else:
                print('b')
                r0 = float(input('r0 [d] = ')) * d
                rp0 = float(input('kąt [rad] = '))
                add_len = float(input('wydłużenie [d] =')) * d
                ray = paraxial_ray_equation.paraxial_ray_equation_RK4(z,av,r0,rp0)
                aray = paraxial_ray_equation.add_ray_after(z,ray,add_len)
                plt.plot(aray[0]/(2*radii[0]),aray[1]/(2*radii[0]))
                plt.xlabel('z/D', fontsize=24)
                plt.ylabel('r/D', fontsize=24)
                plt.xticks(fontsize=18)
                plt.yticks(fontsize=18)
                plt.grid()
                plt.show()
        elif e == 4:
            print('trajektoria z prawej wydłużona')
            print('a: Poglądowa ilość promieni (do wybrania r0, kąt_max, kąt_min, ilość)')
            print('b: Konktretny promień (r0, kąt)')
            choice = input('wybór: ')
            if choice == 'a':
                print('a')
                r0 = float(input('r0 [d] = ')) * d
                rp0max = float(input('kąt maksymalny [rad] = '))
                rp0min = float(input('kąt minimalny [rad] = '))
                n = int(input('liczba promieni = '))
                add_len = float(input('wydłużenie [d] =')) * d
                rray_list = paraxial_ray_equation.many_rays_angles_reverse(z, av, r0, n, rp0max, rp0min)
                rraya_list = []
                for i in range(np.shape(rray_list)[0]):
                    rraya_list.append(paraxial_ray_equation.add_ray_after_reverse(z, rray_list[i], add_len))
                    plt.plot(rraya_list[i][0]/(2*radii[0]),rraya_list[i][1]/(2*radii[0]))
                plt.xlabel('z/D', fontsize=24)
                plt.ylabel('r/D', fontsize=24)
                plt.xticks(fontsize=18)
                plt.yticks(fontsize=18)
                plt.grid()
                plt.show()
            else:
                print('b')
                r0 = float(input('r0 [d] = ')) * d
                rp0 = float(input('kąt [rad] = '))
                add_len = float(input('wydłużenie [d] =')) * d
                rray = paraxial_ray_equation.reverse_paraxial_ray_equation_RK4(z, av, r0, rp0)
                rraya = paraxial_ray_equation.add_ray_after_reverse(z, rray, add_len)
                plt.plot(rraya[0]/(2*radii[0]), rraya[1]/(2*radii[0]))
                plt.xlabel('z/D', fontsize=24)
                plt.ylabel('r/D', fontsize=24)
                plt.xticks(fontsize=18)
                plt.yticks(fontsize=18)
                plt.grid()
                plt.show()
        elif e == 5:
            print('druga pochodna potencjału')
            d2avdz2 = d2_dz2(av)
            plt.xlabel('z/D', fontsize=24)
            plt.plot(z/(2*radii[0]),d2avdz2)
            plt.show()
        elif e == 6:
            print('Przycięcie potencjału')
            finish = 0
            while finish == 0:
                print('wybierz o ile punktów przyciąć potencjał')
                cut_num = int(input('ilość punktów do ucięcia'))
                new_av = potential_fix.fix_replace(av,cut_num)
                d2avdz2 = d2_dz2(av)
                d2_new_av_dz2 = d2_dz2(new_av)
                plt.plot(d2avdz2, c = 'k', label = 'druga pochodna potencjału pierwotnego')
                plt.plot(d2_new_av_dz2, c = 'r', label = f'druga pochodna przycięta {cut_num}')
                plt.xlabel('numer punktu')
                plt.legend()
                plt.grid()
                plt.show()
                print('CZY ZASTOSOWAĆ ZMIANĘ?')
                print('1: TAK')
                print('2: NIE, zmieniam dalej')
                print('3: KONIEC PROGRAMU')
                e = input('wybór: ')
                e = int(e)

                if e == 1:
                    av = new_av
                    finish = 1
                elif e == 2:
                    finish = 0
                elif e == 3:
                    exit()

    return 0


def kontrola_0(z,av,radii):
    end = 'k'
    e = None
    d = 0.01

    d2_dz2 = findiff.FinDiff(0, z[1]-z[0], 2)
    while e != end:
        print('\n\nWYBIERZ CO ROBIĆ DALEJ:\n\n')
        print('1: Oblicz trajektrodię promieni elektronowych lecacą z lewej strony')
        print('2: Oblicz trajektrodię promieni elektronowych lecacą z prawej strony')
        print('3: Oblicz wydłużoną trajektorię promieni z lewej strony')
        print('4: Oblicz wydłużoną trajektorię promieni z prawej strony')
        print('5: Wyświetl druga pochodną potencjału')
        print('k: Zakończ')
        e = input('wybór: ')
        if e=='k':
            exit()
        else:
            e = int(e)

        if e == 1:
            print('trajektoria z lewej')
            print('a: Poglądowa ilość promieni (do wybrania r0, kąt_max, kąt_min, ilość)')
            print('b: Konktretny promień (r0, kąt)')
            choice = input('wybór:')
            if choice == 'a':
                print('a')
                r0 = float(input('r0 [d] = ')) * d
                rp0max = float(input('kąt maksymalny [rad] = '))
                rp0min = float(input('kąt minimalny [rad] = '))
                n = int(input('liczba promieni = '))
                ray_list = paraxial_ray_equation.many_rays_angles(z,av,r0,n,rp0max,rp0min)
                for i in range(len(ray_list)):
                    plt.plot(z/(2*radii[0]),ray_list[i]/(2*radii[0]))
                plt.xlabel('z/D', fontsize=24)
                plt.ylabel('r/D', fontsize=24)
                plt.grid()
                plt.show()


            else:
                print('b')
                r0 = float(input('r0 [d] = ')) * d
                rp0 = float(input('kąt [rad] = '))
                plt.plot(z/(2*radii[0]),paraxial_ray_equation.paraxial_ray_equation_RK4(z,av,r0,rp0)/(2*radii[0]))
                plt.title(f'TRAJETORIA ELEKTRONU r0 = {r0} d, kąt = {rp0} rad')
                plt.xlabel('z/D', fontsize=24)
                plt.ylabel('r/D', fontsize=24)
                plt.grid()
                plt.show()
        elif e == 2:
            print('trajektroria z prawej')
            print('a: Poglądowa ilość promieni (do wybrania r0, kąt_max, kąt_min, ilość)')
            print('b: Konktretny promień (r0, kąt)')
            choice = input('wybór:')
            if choice == 'a':
                print('a')
                r0 = float(input('r0 [d] = ')) * d
                rp0max = float(input('kąt maksymalny [rad] = '))
                rp0min = float(input('kąt minimalny [rad] = '))
                n = int(input('liczba promieni = '))
                rray_list = paraxial_ray_equation.many_rays_angles_reverse(z, av, r0, n, rp0max, rp0min)
                for i in range(len(rray_list)):
                    plt.plot(z[::-1]/(2*radii[0]), rray_list[i]/(2*radii[0]))
                plt.grid()
                plt.xlabel('z/D', fontsize=24)
                plt.ylabel('r/D', fontsize=24)
                plt.show()
            else:
                print('b')
                r0 = float(input('r0 [d] = ')) * d
                rp0 = float(input('kąt [rad] = '))
                plt.plot(z[::-1] / 2*radii[0], paraxial_ray_equation.reverse_paraxial_ray_equation_RK4(z, av, r0, rp0)/(2*radii[0]))
                plt.title(f'TRAJETORIA ELEKTRONU r0 = {r0} d, kąt = {rp0} rad')
                plt.xlabel('z/D', fontsize=24)
                plt.ylabel('r/D', fontsize=24)
                plt.grid()
                plt.show()
        elif e == 3:
            print('trajektrodia z lewej wydłużona')
            print('a: Poglądowa ilość promieni (do wybrania r0, kąt_max, kąt_min, ilość)')
            print('b: Konktretny promień (r0, kąt)')
            choice = input('wybór: ')
            if choice == 'a':
                print('a')
                r0 = float(input('r0 [d] = ')) * d
                rp0max = float(input('kąt maksymalny [rad] = '))
                rp0min = float(input('kąt minimalny [rad] = '))
                n = int(input('liczba promieni = '))
                add_len = float(input('wydłużenie [d] =')) * d
                ray_list = paraxial_ray_equation.many_rays_angles(z, av, r0, n, rp0max, rp0min)
                raya_list = []
                for i in range(np.shape(ray_list)[0]):
                    raya_list.append(paraxial_ray_equation.add_ray_after(z, ray_list[i], add_len))
                    plt.plot(raya_list[i][0]/(2*radii[0]), raya_list[i][1]/(2*radii[0]))
                plt.grid()
                plt.xlabel('z/D', fontsize=24)
                plt.ylabel('r/D', fontsize=24)
                plt.show()

            else:
                print('b')
                r0 = float(input('r0 [d] = ')) * d
                rp0 = float(input('kąt [rad] = '))
                add_len = float(input('wydłużenie [d] =')) * d
                ray = paraxial_ray_equation.paraxial_ray_equation_RK4(z,av,r0,rp0)
                aray = paraxial_ray_equation.add_ray_after(z,ray,add_len)
                plt.plot(aray[0]/(2*radii[0]),aray[1]/(2*radii[0]))
                plt.xlabel('z/D', fontsize=24)
                plt.ylabel('r/D', fontsize=24)
                plt.grid()
                plt.show()
        elif e == 4:
            print('trajektoria z prawej wydłużona')
            print('a: Poglądowa ilość promieni (do wybrania r0, kąt_max, kąt_min, ilość)')
            print('b: Konktretny promień (r0, kąt)')
            choice = input('wybór: ')
            if choice == 'a':
                print('a')
                r0 = float(input('r0 [d] = ')) * d
                rp0max = float(input('kąt maksymalny [rad] = '))
                rp0min = float(input('kąt minimalny [rad] = '))
                n = int(input('liczba promieni = '))
                add_len = float(input('wydłużenie [d] =')) * d
                rray_list = paraxial_ray_equation.many_rays_angles_reverse(z, av, r0, n, rp0max, rp0min)
                rraya_list = []
                for i in range(np.shape(rray_list)[0]):
                    rraya_list.append(paraxial_ray_equation.add_ray_after_reverse(z, rray_list[i], add_len))
                    plt.plot(rraya_list[i][0]/(2*radii[0]),rraya_list[i][1]/(2*radii[0]))
                plt.grid()
                plt.xlabel('z/D', fontsize=24)
                plt.ylabel('r/D', fontsize=24)
                plt.show()
            else:
                print('b')
                r0 = float(input('r0 [d] = ')) * d
                rp0 = float(input('kąt [rad] = '))
                add_len = float(input('wydłużenie [d] =')) * d
                rray = paraxial_ray_equation.reverse_paraxial_ray_equation_RK4(z, av, r0, rp0)
                rraya = paraxial_ray_equation.add_ray_after_reverse(z, rray, add_len)
                plt.plot(rraya[0]/(2*radii[0]), rraya[1]/(2*radii[0]))
                plt.xlabel('z/D', fontsize=24)
                plt.ylabel('r/D', fontsize=24)
                plt.grid()
                plt.show()
        elif e == 5:
            print('druga pochodna potencjału')
            d2avdz2 = d2_dz2(av)
            plt.plot(z/(2*radii[0]),d2avdz2)
            plt.xlabel('z/D', fontsize=24)
            plt.show()


    return 0


def kontrola_2(z,avn,avi,radii):
    end = 'k'
    e = None
    d = 0.01

    d2_dz2 = findiff.FinDiff(0, z[1]-z[0], 2)
    while e != end:
        print('\n\nWYBIERZ CO ROBIĆ DALEJ:\n\n')
        print('1: Oblicz trajektrodię promieni elektronowych lecacą z lewej strony')
        print('2: Oblicz trajektrodię promieni elektronowych lecacą z prawej strony')
        print('3: Oblicz wydłużoną trajektorię promieni z lewej strony')
        print('4: Oblicz wydłużoną trajektorię promieni z prawej strony')
        print('5: Wyświetl druga pochodną potencjałów')
        print('6: Przytnij potencjał/y i zastąp stałymi wartościami - pomaga usunąć błędy na brzegach')
        print('k: Zakończ')
        e = input('wybór: ')
        if e=='k':
            exit()
        else:
            e = int(e)

        if e == 1:
            print('trajektoria z lewej')
            print('a: Po dwa promienie z danego potencjału (do wybrania r0, kąt_max, kąt_min, ilość)')
            print('b: Konktretny promień (r0, kąt), dla obu')
            choice = input('wybór:')
            if choice == 'a':
                print('a')
                r0 = float(input('r0 [d] = ')) * d
                rp01 = float(input('kąt_1  [rad] = '))
                rp02 = float(input('kąt_2  [rad] = '))

                rayn1 = paraxial_ray_equation.paraxial_ray_equation_RK4(z,avn,r0,rp01)
                rayn2 = paraxial_ray_equation.paraxial_ray_equation_RK4(z,avn,r0,rp02)
                rayi1 = paraxial_ray_equation.paraxial_ray_equation_RK4(z,avi,r0,rp01)
                rayi2 = paraxial_ray_equation.paraxial_ray_equation_RK4(z,avi,r0,rp02)

                plt.plot(z/(2*radii[0]),rayn1, c = 'm', label = 'NORMALNE CDM')
                plt.plot(z/(2*radii[0]),rayn2, c = 'm')
                plt.plot(z /(2*radii[0]), rayi1, c='g', label='ULEPSZONE CDM')
                plt.plot(z /(2*radii[0]), rayi2, c='g')
                plt.xlabel('z/D')
                plt.ylabel('r/D')
                plt.legend()
                plt.grid()
                plt.show()


            else:
                print('b')
                r0 = float(input('r0 [d] = ')) * d
                rp0 = float(input('kąt [rad] = '))
                rayn = paraxial_ray_equation.paraxial_ray_equation_RK4(z,avn,r0,rp0)
                rayi = paraxial_ray_equation.paraxial_ray_equation_RK4(z,avi,r0,rp0)
                plt.plot(z/(2*radii[0]),rayn/(2*radii[0]), label = 'NORMALNE CDM')
                plt.plot(z/(2*radii[0]),rayi/(2*radii[0]), label = 'ULEPSZONE CDM')
                plt.xlabel('z/D', fontsize=24)
                plt.ylabel('r/D', fontsize=24)
                plt.title(f'TRAJETORIA ELEKTRONU r0 = {r0} d, kąt = {rp0} rad')
                plt.grid()
                plt.legend()
                plt.show()
        elif e == 2:
            print('trajektroria z prawej')
            print('a: Dwa promienie')
            print('b: Konktretny promień (r0, kąt)')
            choice = input('wybór:')
            if choice == 'a':
                print('a')
                r0 = float(input('r0 [d] = ')) * d
                rp01 = float(input('kąt_1 [rad] = '))
                rp02 = float(input('kąt_2 [rad] = '))

                rayn1 = paraxial_ray_equation.reverse_paraxial_ray_equation_RK4(z, avn, r0, rp01)
                rayn2 = paraxial_ray_equation.reverse_paraxial_ray_equation_RK4(z, avn, r0, rp02)
                rayi1 = paraxial_ray_equation.reverse_paraxial_ray_equation_RK4(z, avi, r0, rp01)
                rayi2 = paraxial_ray_equation.reverse_paraxial_ray_equation_RK4(z, avi, r0, rp02)

                plt.plot(z[::-1] / 2*radii[0], rayn1, c='m', label='NORMALNE CDM')
                plt.plot(z[::-1] / 2*radii[0], rayn2, c='m')
                plt.plot(z[::-1] / 2*radii[0], rayi1, c='g', label='ULEPSZONE CDM')
                plt.plot(z[::-1]/ 2*radii[0], rayi2, c='g')
                plt.xlabel('z/D', fontsize=24)
                plt.ylabel('r/D', fontsize=24)
                plt.legend()
                plt.grid()
                plt.show()
            else:
                print('b')
                r0 = float(input('r0 [d] = ')) * d
                rp0 = float(input('kąt [rad] = '))
                rayn = paraxial_ray_equation.reverse_paraxial_ray_equation_RK4(z, avn, r0, rp0)
                rayi = paraxial_ray_equation.reverse_paraxial_ray_equation_RK4(z, avi, r0, rp0)
                plt.plot(z[::-1] / 2*radii[0], rayn / 2*radii[0], label='NORMALNE CDM')
                plt.plot(z[::-1] / 2*radii[0], rayi / 2*radii[0], label='ULEPSZONE CDM')
                plt.xlabel('z/D', fontsize=24)
                plt.ylabel('r/D', fontsize=24)
                plt.title(f'TRAJETORIA ELEKTRONU r0 = {r0} d, rp0 = {rp0} rad')
                plt.grid()
                plt.legend()
                plt.show()
        elif e == 3:
            print('trajektrodia z lewej wydłużona')
            print('a: Poglądowa ilość promieni (do wybrania r0, kąt_max, kąt_min, ilość)')
            print('b: Konktretny promień (r0, kąt)')
            choice = input('wybór: ')
            if choice == 'a':
                print('a')
                r0 = float(input('r0 [d] = ')) * d
                rp01 = float(input('kąt_1 [rad] = '))
                rp02 = float(input('kąt_2 [rad]= '))
                add_len = float(input('wydłużenie [d] =')) * d

                rayn1 = paraxial_ray_equation.paraxial_ray_equation_RK4(z, avn, r0, rp01)
                rayn2 = paraxial_ray_equation.paraxial_ray_equation_RK4(z, avn, r0, rp02)
                rayi1 = paraxial_ray_equation.paraxial_ray_equation_RK4(z, avi, r0, rp01)
                rayi2 = paraxial_ray_equation.paraxial_ray_equation_RK4(z, avi, r0, rp02)

                arayn1 = paraxial_ray_equation.add_ray_after(z,rayn1,add_len)
                arayn2 = paraxial_ray_equation.add_ray_after(z,rayn2,add_len)
                arayi1 = paraxial_ray_equation.add_ray_after(z,rayi1,add_len)
                arayi2 = paraxial_ray_equation.add_ray_after(z,rayi2,add_len)

                plt.plot(arayn1[0]/(2*radii[0]),arayn1[1]/2*radii[0], c = 'm', label = 'NORMALNE CDM')
                plt.plot(arayn2[0]/(2*radii[0]),arayn2[1]/2*radii[0], c = 'm')
                plt.plot(arayi1[0]/(2*radii[0]),arayi1[1]/2*radii[0], c = 'g', label = 'ULEPSZONE CDM')
                plt.plot(arayi1[0]/(2*radii[0]),arayi2[1]/2*radii[0], c = 'g')
                plt.xlabel('z/D', fontsize=24)
                plt.ylabel('r/D', fontsize=24)
                plt.grid()
                plt.legend()
                plt.show()

            else:
                print('b')
                r0 = float(input('r0 [d] = ')) * d
                rp0 = float(input('kąt [rad] = '))
                add_len = float(input('wydłużenie [d] =')) * d
                rayn = paraxial_ray_equation.paraxial_ray_equation_RK4(z,avn,r0,rp0)
                rayi = paraxial_ray_equation.paraxial_ray_equation_RK4(z,avi,r0,rp0)
                arayn = paraxial_ray_equation.add_ray_after(z,rayn,add_len)
                arayi = paraxial_ray_equation.add_ray_after(z,rayi,add_len)
                plt.plot(arayn[0]/(2*radii[0]),arayn[1]/(2*radii[0]), c='m', label = 'NORMALNE CDM')
                plt.plot(arayi[0]/(2*radii[0]),arayi[1]/(2*radii[0]), c='g', label = 'ULEPSZONE CDM')
                plt.xlabel('z/D', fontsize=24)
                plt.ylabel('r/D', fontsize=24)
                plt.grid()
                plt.legend()
                plt.show()
        elif e == 4:
            print('trajektoria z prawej wydłużona')
            print('a: Poglądowa ilość promieni (do wybrania r0, kąt_max, kąt_min, ilość)')
            print('b: Konktretny promień (r0, kąt)')
            choice = input('wybór: ')
            if choice == 'a':
                print('a')
                r0 = float(input('r0 [d] = ')) * d
                rp01 = float(input('kąt_1 [rad] = '))
                rp02 = float(input('kąt_2 [rad] = '))
                add_len = float(input('wydłużenie [d] =')) * d

                rayn1 = paraxial_ray_equation.reverse_paraxial_ray_equation_RK4(z, avn, r0, rp01)
                rayn2 = paraxial_ray_equation.reverse_paraxial_ray_equation_RK4(z, avn, r0, rp02)
                rayi1 = paraxial_ray_equation.reverse_paraxial_ray_equation_RK4(z, avi, r0, rp01)
                rayi2 = paraxial_ray_equation.reverse_paraxial_ray_equation_RK4(z, avi, r0, rp02)

                arayn1 = paraxial_ray_equation.add_ray_after_reverse(z, rayn1, add_len)
                arayn2 = paraxial_ray_equation.add_ray_after_reverse(z, rayn2, add_len)
                arayi1 = paraxial_ray_equation.add_ray_after_reverse(z, rayi1, add_len)
                arayi2 = paraxial_ray_equation.add_ray_after_reverse(z, rayi2, add_len)

                plt.plot(arayn1[0] / 2*radii[0], arayn1[1] / 2 * (radii[0]), c='m', label='NORMALNE CDM')
                plt.plot(arayn2[0] / 2*radii[0], arayn2[1] / 2 * (radii[0]), c='m')
                plt.plot(arayi1[0] / 2*radii[0], arayi1[1] / 2 * (radii[0]), c='g', label='ULEPSZONE CDM')
                plt.plot(arayi1[0] / 2*radii[0], arayi2[1] / 2 * (radii[0]), c='g')
                plt.xlabel('z/D', fontsize=24)
                plt.ylabel('r/D', fontsize=24)
                plt.grid()
                plt.legend()
                plt.show()
            else:
                print('b')
                r0 = float(input('r0 [d] = ')) * d
                rp0 = float(input('kąt [rad] = '))
                add_len = float(input('wydłużenie [d] =')) * d
                rayn = paraxial_ray_equation.reverse_paraxial_ray_equation_RK4(z, avn, r0, rp0)
                rayi = paraxial_ray_equation.reverse_paraxial_ray_equation_RK4(z, avi, r0, rp0)
                arayn = paraxial_ray_equation.add_ray_after_reverse(z, rayn, add_len)
                arayi = paraxial_ray_equation.add_ray_after_reverse(z, rayi, add_len)
                plt.plot(arayn[0] / 2*radii[0], arayn[1] / 2*radii[0], c='m', label='NORMALNE CDM')
                plt.plot(arayi[0] / 2*radii[0], arayi[1] / 2*radii[0], c='g', label='ULEPSZONE CDM')
                plt.xlabel('z/D', fontsize=24)
                plt.ylabel('r/D', fontsize=24)
                plt.grid()
                plt.legend()
                plt.show()
        elif e == 5:
            print('druga pochodna potencjału')
            d2avndz2 = d2_dz2(avn)
            d2avidz2 = d2_dz2(avi)
            plt.plot(z/(2*radii[0]),d2avndz2, c ='m', label='NORMALNE')
            plt.plot(z/(2*radii[0]),d2avidz2, c ='g', label='ULEPSZONE')
            plt.grid()
            plt.xlabel('z/D', fontsize=24)
            plt.legend()
            plt.show()
        elif e == 6:
            print('Przycięcie potencjałów')
            finish = 0
            while finish == 0:
                print('KTÓRY POTENCJAŁ CHCESZ PODCIĄĆ')
                print('1:NORMALNY')
                print('2:ULEPSZONY')
                c1 = int(input('wybór = '))
                if c1 == 1:
                    print('NORMALNY')
                    print('wybierz o ile punktów przyciąć potencjał')
                    cut_num = int(input('ilość punktów do ucięcia'))
                    new_av = potential_fix.fix_replace(avn, cut_num)
                    d2avndz2 = d2_dz2(avn)
                    d2_new_av_dz2 = d2_dz2(new_av)
                    plt.plot(d2avndz2, c='k', label='druga pochodna potencjału pierwotnego')
                    plt.plot(d2_new_av_dz2, c='r', label=f'druga pochodna przycięta {cut_num}')
                    plt.xlabel('numer punktu')
                    plt.legend()
                    plt.grid()
                    plt.show()
                    print('CZY ZASTOSOWAĆ ZMIANĘ?')
                    print('1: TAK')
                    print('2: NIE, zmieniam dalej')
                    print('3: KONIEC PROGRAMU')
                    e = input('wybór: ')
                    e = int(e)

                    if e == 1:
                        avn = new_av
                        finish = 1
                    elif e == 2:
                        finish = 0
                    elif e == 3:
                        exit()
                elif c1 == 2:
                    print('ULEPSZONY')
                    print('wybierz o ile punktów przyciąć potencjał')
                    cut_num = int(input('ilość punktów do ucięcia'))
                    new_avi = potential_fix.fix_replace(avi, cut_num)
                    d2avidz2 = d2_dz2(avi)
                    d2_new_av_dz2 = d2_dz2(new_avi)
                    plt.plot(d2avidz2, c='k', label='druga pochodna potencjału pierwotnego')
                    plt.plot(d2_new_av_dz2, c='r', label=f'druga pochodna przycięta {cut_num}')
                    plt.xlabel('numer punktu')
                    plt.legend()
                    plt.grid()
                    plt.show()
                    print('CZY ZASTOSOWAĆ ZMIANĘ?')
                    print('1: TAK')
                    print('2: NIE, zmieniam dalej')
                    print('3: KONIEC PROGRAMU')
                    e = input('wybór: ')
                    e = int(e)

                    if e == 1:
                        avi = new_avi
                        finish = 1
                    elif e == 2:
                        finish = 0
                    elif e == 3:
                        exit()
    return 0