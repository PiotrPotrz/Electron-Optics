import numpy as np
import matplotlib.pyplot as plt
import findiff

def fix_replace(potential: np.array, number: int):
    new_potential = np.array(np.ones(len(potential)))
    for i in range(len(potential)):
        if i<number:
            new_potential[i] = potential[number]
        elif i>=len(potential)-number:
            new_potential[i] = potential[len(potential) - number]
        else:
            new_potential[i] = potential[i]
    return new_potential


"""d = 0.01

av1 = np.genfromtxt('C:\\Users\Piotr\\PycharmProjects\\CIRCULAR APERTURE\\wyniki\\[0.05, 0.05],[0.001],[0.005, 0.005],[5, 6],200,400,-0.05,0.05,normal.txt',delimiter=',')
plt.plot(av1)
plt.show()
z = np.linspace( - 5 * d, 5 * d, 200)

d2_dz2 = findiff.FinDiff(0, z[1] - z[0], 2)
plt.plot(d2_dz2(av1))
plt.show()

av1_fixed = fix_replace(av1,25)

plt.plot(av1_fixed)
plt.show()



plt.plot(d2_dz2(av1_fixed))
plt.grid()
plt.show()

ava = two_cylinder_lens.axial_potential_general(5,6,0.5*d,z,0.1*d)

ray_a1 = paraxial_ray_equation.paraxial_ray_equation_RK4(z,ava,0.001 * d, 0)
ray_a2 = paraxial_ray_equation.paraxial_ray_equation_RK4(z,ava,0.001 * d, 0.01)

ray_1 = paraxial_ray_equation.paraxial_ray_equation_RK4(z,av1_fixed,0.001 * d, 0)
ray_2 = paraxial_ray_equation.paraxial_ray_equation_RK4(z,av1_fixed,0.001 * d, 0.01)

plt.plot(z/d, ray_a1/d,label = 'analityczny1')
plt.plot(z/d, ray_a2/d,label = 'analityczny2')
plt.plot(z/d, ray_1/d,label = 'numfix1')
plt.plot(z/d, ray_2/d,label = 'numfix2')
plt.legend()
plt.grid()
plt.show()


ap_raya1 = paraxial_ray_equation.add_ray_after(z,ray_a1,200*d)
ap_raya2 = paraxial_ray_equation.add_ray_after(z,ray_a2,200*d)
ap_ray1 = paraxial_ray_equation.add_ray_after(z,ray_1,200*d)
ap_ray2 = paraxial_ray_equation.add_ray_after(z,ray_2,200*d)


plt.plot(ap_ray1[0]/d, ap_raya1[1]/d,label = 'analityczny1')
plt.plot(ap_ray1[0]/d, ap_raya2[1]/d,label = 'analityczny2')
plt.plot(ap_ray1[0]/d, ap_ray1[1]/d,label = 'numfix1')
plt.plot(ap_ray1[0]/d, ap_ray2[1]/d,label = 'numfix2')
plt.legend()
plt.grid()
plt.show()"""