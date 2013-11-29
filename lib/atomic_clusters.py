import numpy
import math
from schrodinger import by_diagonalization
from mean_field import effective_potentials


'''
define and implement decent output

!!!search for 'disgusting'!!!

better naming here and there


implement energy2 function
implement automatic l_max selection
implement keeping of needed states and eigenvalues only in l_max != -1 case
implement filling of states -> (updated to) implement check on magick numbers -> (updated to) implemented some rough way
'''

class Rho:
    def __init__(self, mesh, step, beta = 0):
        self.mesh = mesh
        self.step = step
        self.beta = beta
    def compute(self, occupied_states, degeneracies, old_rho):
        rho = numpy.zeros(len(self.mesh))
        for i in range(len(self.mesh)):
            for j in range(len(degeneracies)):#len(degeneracies) is not that elegant
                rho[i] = rho[i] + degeneracies[j] * math.pow(occupied_states[:, j][i], 2)
        if self.beta == 0:
            return rho / (4 * math.pi)
        else:
            return self.beta * rho / (4 * math.pi) + (1 - self.beta) * old_rho

class Energy1:
    def __init__(self, mesh, step):
        self.mesh = mesh
        self.step = step
    def compute(self, occupied_eigenvalues, degeneracies, rho, hartree_potential, exchange_contribute, correlation_contribute):
        return numpy.sum(eigenvalues * degeneracies) - 2 * math.pi * numpy.sum(self.mesh * self.mesh * rho * hartree_potential) * self.step - 4 * math.pi * numpy.sum(self.mesh * self.mesh * rho * rho * (exchange_contribute + correlation_contribute)) * self.step# is the exchange and correlation potential correct?

'''
class Energy2:
    def __init__(self, mesh, step, external_potential):
        self.mesh = mesh
        self.step = step
        self.external_potential = external_potential
    def r_laplacian(r_function):# better naming needed
        d2r_r_function = numpy.empty(len(r_function))
        for i in range(1, len(state) - 1):
            d2r_state = (r_function[i + 1] - 2 * r_function[i] + r_function[i - 1])
        # 0th order approx
        d2r_state[0] = d2r_state[1]
        d2r_state[len(d2r_state) - 1] = d2r_state[len(d2r_state) - 2]
        return d2r_state / self.step**2

    def energy(self, occupied_states, angular_momentum, degeneracies):
        kinetic_energy_radial_density = numpy.zeros(len(self.mesh))
        for i in range(len(angular_momentum)):# len(angular_momentum) is not that elegant
            kinetic_energy_radial_density = kinetic_energy_radial_density - degeneracies[i] * 0.5 * (occupied_states[:, i] * self.r_laplacian(self.mesh * occupied_states[:, i]) / self.mesh - 1 / (self.mesh * self.mesh) * angular_momentum[i] * (angular_momentum[i] + 1))
        kinetic_energy = numpy.sum(kinetic_energy_radial_density * self.mesh * self.mesh) * self.step

        interaction_energy_radial_density
        for i in range(len(angular_momentum)):
            for j in range(i):
                sfd
# scipy.integrate and scipy.spherical_harmonics are your friends
'''


if __name__ == '__main__':
#SETTINGS
#mesh settings
    R = 20
    step = 0.02

#number of electrons and related quantities settings
    N = 40

#    electronic_configuration = numpy.array([1, 1])# 8 atoms
#    electronic_configuration = numpy.array([2, 1, 1])# 20 atoms
    electronic_configuration = numpy.array([2, 2, 1, 1])# 40 atoms

    l_max = - 1# set to - 1 if you give set the electronic configuration manually

    wigner_seitz_radius = 3.93
    R_c = math.pow(N, 1 / 3) * wigner_seitz_radius
    positive_charge_density = N / (4 / 3 * math.pi * math.pow(R_c, 3))

#initial density settings
    mu = 1

#stopping criteria
    single_energy_epsilon = 0.001
#    double_energies_epsilon = 0.1
#    ITERATIONS = 1

#self consistency
    beta = 0.25


#external potential
    def v_ext(mesh):
        external_potential = numpy.empty(len(mesh))
        for i in range(len(mesh)):
            if mesh[i] <= R_c:
                external_potential[i] = - 2 * math.pi * positive_charge_density * (math.pow(R_c, 2) - 1 / 3 * math.pow(mesh[i], 2))
            else:
                external_potential[i] = - 4 / 3 * math.pi * positive_charge_density * math.pow(R_c, 3) / mesh[i]
        return external_potential

#initial density
    def initial_density(mesh, normalization_constant = 1):
        rho = numpy.empty(len(mesh))
        for i in range(len(mesh)):
            rho[i] = normalization_constant / (1 + math.exp(- mu * (R_c - mesh[i])))
        return rho

    def density_normalization(mesh, step, rho):# generalization to arbitrarly spaced meshes
        return 1 / (4 * math.pi * numpy.sum(rho * mesh * mesh) * step)


    mesh = numpy.arange(step, R, step)#careful, the last point is missing

    external_potential = v_ext(mesh)

    rho = initial_density(mesh, 1)
    initial_density_normalization_constant = density_normalization(mesh, step, rho)
    rho = initial_density(mesh, N * initial_density_normalization_constant)

    lda_kohn_sham = effective_potentials.LdaKohnSham(mesh, step, external_potential)# think about naming
    effective_potential = lda_kohn_sham.kohn_sham_potential# effective_potential(rho, l)

    rho_compute = Rho(mesh, step, beta).compute

#    j = ITERATIONS
#    while j > 0:

    energy_compute = Energy1(mesh,step).compute
    old_energy = 0# i doubled the condition on the while cycle to allow this
    energy_change = 0
    j = 0
    while energy_change > single_energy_epsilon or j < 2:

        eigenvalues = numpy.empty(0)
        eigenvectors = numpy.empty( (len(mesh), 0) )
        angular_momentum = numpy.empty(0)

        if l_max != - 1:
            for l in range(l_max + 1):
                tmp_eigenvalues, tmp_eigenvectors = by_diagonalization.solve(mesh, step, effective_potential(rho, l))

                eigenv_needed = int(N / (2 * (2 * l + 1)))
                tmp_eigenvalues = tmp_eigenvalues[0:eigenv_needed]#is int() good enough, you might want to keep int() + 1
                tmp_eigenvectors = tmp_eigenvectors[:,0:eigenv_needed]

                for i in range(len(tmp_eigenvalues)):#not elegant
                    tmp_eigenvectors[:, i] = tmp_eigenvectors[:, i] / mesh

                eigenvalues = numpy.append(eigenvalues, tmp_eigenvalues)
                eigenvectors = numpy.append(eigenvectors, tmp_eigenvectors, 1)
                angular_momentum = numpy.append(angular_momentum, l * numpy.ones(len(tmp_eigenvalues)))

        else:
            for l in range(len(electronic_configuration)):
                tmp_eigenvalues, tmp_eigenvectors = by_diagonalization.solve(mesh, step, effective_potential(rho, l))

                tmp_eigenvalues = tmp_eigenvalues[0:electronic_configuration[l]]#is int() good enough, you might want to keep int() + 1
                tmp_eigenvectors = tmp_eigenvectors[:,0:electronic_configuration[l]]

                for i in range(len(tmp_eigenvalues)):#not elegant
                    tmp_eigenvectors[:, i] = tmp_eigenvectors[:, i] / mesh

                eigenvalues = numpy.append(eigenvalues, tmp_eigenvalues)
                eigenvectors = numpy.append(eigenvectors, tmp_eigenvectors, 1)
                angular_momentum = numpy.append(angular_momentum, l * numpy.ones(len(tmp_eigenvalues)))

        #this block is needed only in the l_max != -1 case, left it here for clarity
        idx = eigenvalues.argsort()
        eigenvalues = eigenvalues[idx]
        eigenvectors = eigenvectors[:, idx]
        angular_momentum = angular_momentum[idx]

        #this could be optimized allocating the space needed instead of appending things later
        unassigned_particles = N
        occupied_eigenvalues = numpy.empty(0)
        occupied_eigenvectors = numpy.empty( (len(mesh), 0) )
        degeneracies = numpy.empty(0)

        i = 0
        while unassigned_particles > 0:
            occupied_eigenvalues = numpy.append(occupied_eigenvalues, eigenvalues[i])
            occupied_eigenvectors = numpy.append(occupied_eigenvectors, numpy.reshape(eigenvectors[:, i], (len(mesh), 1)), 1)
            degeneracies = numpy.append(degeneracies, 2 * (2 * angular_momentum[i] + 1))
            unassigned_particles = unassigned_particles - degeneracies[i]
            i = i + 1

        if unassigned_particles != 0:
            print('not a closed shell configuration!')
            exit()

        rho = rho_compute(occupied_eigenvectors, degeneracies, rho)

#        j = j - 1

        current_energy = energy_compute(occupied_eigenvalues, degeneracies, rho, lda_kohn_sham.hartree_potential(rho), lda_kohn_sham.lda_exchange_density_rho_derivative(rho), lda_kohn_sham.lda_correlation_density_rho_derivative(rho))
        try:
            energy_change = abs(old_energy - current_energy)# is it possible to exploit the fact that GS energy is lowest?
        except:
            pass
        old_energy = current_energy
        print(current_energy)

        j = j + 1

    output = open(str(N) + 'at' + '_epsilon' + str(single_energy_epsilon) + '_beta' + str(beta) + '.dat', 'w')
    output.write('# energy [hartree]:' + str(current_energy) + '\n')
    output.write('# iterations before stopping: ' + str(j) + '\n')
    output.write('#Runtime settings\n')
    output.write('#number of atoms: ' + str(N) + '\n')
    output.write('#electronic configuration: ' + str(electronic_configuration) + '\n')
    output.write('#wiegner seitz radious: ' + str(wigner_seitz_radius) + '\n')
    output.write('#R_c: ' + str(R_c) + '\n')
    output.write('#beta: ' + str(beta) + '\n\n')
    output.write('#r\trho\tpotential\n')
    potential = effective_potential(rho, 0)
    for i in range(len(mesh)):
        output.write(str(mesh[i]) + '\t' + str(rho[i]) + '\t' + str(potential[i]) + '\n')

#    print(density_normalization(mesh, step, rho))
