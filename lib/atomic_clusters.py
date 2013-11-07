import numpy
import math
from schrodinger import by_diagonalization
from mean_field import effective_potentials


'''
change variable name local_rho
change name correlation and exchange potentials

!!!implement averaging of rhos
!!!implement energy functions

implement automatic l_max selection
implement keeping of needed states and eigenvalues only
implement filling of states -> (updated to) implement check on magick numbers
'''


#SETTINGS

#mesh settings
R = 16
step = 0.016

#number of electrons and related quantities settings
N = 8

#shell_structure = numpy.array( [1,1] )

wigner_seitz_radius = 4
R_c = math.pow(N, 1 / 3) * wigner_seitz_radius
positive_charge_density = N / (4 / 3 * math.pi * math.pow(R_c, 3))

l_max = 1#WILL YOU TRY TO FIND A WAY TO CALCULATE THIS FROM N? (Not sure you can though, eigenstates' order in terms of eigenvalues depends on the potential).

#initial density settings
mu = 1

#stopping criteria
#epsilon = 0.01
ITERATIONS = 2

#self consistency
beta = 0.5

'''
def zero(x):
    return 0
    '''

class Rho:
    def __init__(self, mesh, step, occupied_states, degeneracies, beta = 0, old_rho = 'you dont draw a shit lebowsky'):#how do i give a function as a default?
        self.mesh = mesh
        self.step = step
        self.occupied_states = occupied_states
        self.degeneracies = degeneracies
        self.beta = beta
        self.old_rho = old_rho
    def compute(self, x):
        local_rho = 0
        for i in range(len(self.degeneracies)):#len(degeneracies) is not that elegant
            local_rho = local_rho + self.degeneracies[i] * math.pow(self.occupied_states[:, i][x / self.step - 1], 2)
        if self.beta == 0:
            return local_rho / (4 * math.pi)
        else:
            return self.beta * local_rho / (4 * math.pi) + (1 - self.beta) * self.old_rho(x)

'''
class Energy1:
    def __init__(self, eigenvalues, degeneracies, rho, hartree_potential, exchange_potential, correlation_potential):
        self.eigenvalues = eigenvalues
    def compute(self):
        return numpy.sum(self.eigenvalues * self.degeneracies) #to finish

class Energy2:
    def energy():
        pass
'''

if __name__ == '__main__':
    class LdaExchangePotential:
        def __init__(self, rho):
            self.rho = rho
        def compute(self, x):
            return - math.pow(3 * self.rho(x) / math.pi, 1 / 3)

    class LdaCorrelationPotential:
        def __init__(self, rho):
            self.rho = rho
        def compute(self, x):#to optimize
            return - 0.44 / (4 * math.pi) * 1 / math.pow(7.8 + math.pow(3 / (4 * math.pi * self.rho(x)), 1 / 3), 2) * 1 / math.pow(3 / (4 * math.pi * self.rho(x)), 2 / 3) * 1 / self.rho(x) - 0.44 / (7.8 + math.pow(3 / (4 * math.pi * self.rho(x)), 1 / 3))

    class ExternalPotential:
        def compute(x):
            if x <= R_c:
                return - 2 * math.pi * positive_charge_density * (math.pow(R_c, 2) - 1 / 3 * math.pow(x, 2))
            else:
                return - 4 / 3 * math.pi * positive_charge_density * math.pow(R_c, 3) / x

    class InitialDensity:
        def __init__(self, normalization_constant):
            self.normalization_constant = normalization_constant
        def compute(self, x):
            return self.normalization_constant / (1 + math.exp(- mu * (R_c - 0.5 - x)))

    def density_normalization(mesh, step, rho):#better naming needed
        density = numpy.empty(len(mesh))#find some better naming
        for i in range(len(mesh)):
            density[i] = rho(mesh[i])
        return 1 / (4 * math.pi * numpy.sum(density * mesh * mesh) * step)#4 pi thing


    mesh = numpy.arange(step, R, step)#careful, the last point is missing

    rho = InitialDensity(1)
    density_normalization_constant = density_normalization(mesh, step, rho.compute)
    rho = InitialDensity(N * density_normalization_constant)

    external_potential = ExternalPotential

#    while energy_discrepancy > epsilon:
    j = ITERATIONS
    while j > 0:


        eigenvalues = numpy.empty(0)
        eigenvectors = numpy.empty( (len(mesh), 0) )
        angular_momentum = numpy.empty(0)

        exchange_potential = LdaExchangePotential(rho.compute)
        correlation_potential = LdaCorrelationPotential(rho.compute)

        for l in range(l_max + 1):
            effective_potential = effective_potentials.KohnSham(mesh, step, rho.compute, exchange_potential.compute, correlation_potential.compute, external_potential.compute, l)
            tmp_eigenvalues, tmp_eigenvectors = by_diagonalization.solve(mesh, step, effective_potential.compute)

            eigenv_needed = int(N / (2 * (2 * l + 1)))

            tmp_eigenvalues = tmp_eigenvalues[0:eigenv_needed]#is int() good enough, you might want to keep int() + 1
            tmp_eigenvectors = tmp_eigenvectors[:,0:eigenv_needed]

            for i in range(len(tmp_eigenvalues)):#not elegant
                tmp_eigenvectors[:, i] = tmp_eigenvectors[:, i] / mesh

            eigenvalues = numpy.append(eigenvalues, tmp_eigenvalues)
            eigenvectors = numpy.append(eigenvectors, tmp_eigenvectors, 1)
            angular_momentum = numpy.append(angular_momentum, l * numpy.ones(len(tmp_eigenvalues)))

        idx = eigenvalues.argsort()
        eigenvalues = eigenvalues[idx]
        eigenvectors = eigenvectors[:, idx]
        angular_momentum = angular_momentum[idx]

        #you should check wether N is a magick number at some stage
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


#        rho = Rho(mesh, step, occupied_eigenvectors, degeneracies, beta, rho.compute)

        j = j - 1

#        energy_discrepancy = 
#    rho = Rho(mesh, step, occupied_eigenvectors, degeneracies)
    for x in mesh:
        print(x, rho.compute(x))


#    print(density_normalization(mesh, step, rho.compute))


    '''
    exchange_potential = LdaExchangePotential(rho.compute)
    correlation_potential = LdaCorrelationPotential(rho.compute)
    effective_potential = effective_potentials.KohnSham(mesh, step, rho.compute, exchange_potential.compute, correlation_potential.compute, external_potential.compute, 0)
    for x in mesh:
        print(x, effective_potential.compute(x))
        '''
