import numpy
import math

'''
 - INTERACTION POTENTIAL THING!
 - the way the centrifugal potential is computed is confusing
 - some potential() and some potential.compute()
'''
'''
class Centrifugal:
    def __init__(self, l):
        self.l = l
    def compute(self, x):
        return self.l * (self.l + 1) / math.pow(x, 2)

class Hartree:
    def __init__(self, mesh, step, rho, l):
        self.mesh = mesh
        self.step = step
        self.rho = rho
#        self.interaction_potential = interaction_potential
        self.l = l
    def compute(self, x):# it is not implemented in the most general case yet, it just holds for coulomb interaction and spherical symmetry
        """
        Returns the Hartree potential integrated some rough way
        """
        centrifugal_potential = Centrifugal(self.l)
        hartree_potential = 0
        for y in self.mesh:
            hartree_potential = hartree_potential + self.rho(y) * (abs(x - y) - (x + y)) * y / x * self.step
        return - 2 * math.pi * hartree_potential + centrifugal_potential.compute(x)

class KohnSham:
    def __init__(self, mesh, step, rho, exchange_potential, correlation_potential, external_potential, l):
        self.mesh = mesh
        self.step = step
        self.rho = rho
#        self.interaction_potential = interaction_potential
        self.exchange_potential = exchange_potential
        self.correlation_potential = correlation_potential
        self.external_potential = external_potential
        self.l = l
    def compute(self, x):
        hartree_potential = Hartree(self.mesh, self.step, self.rho, self.l)
        return hartree_potential.compute(x) + self.exchange_potential(x) + self.correlation_potential(x) + self.external_potential(x)
'''

class LdaKohnSham:
    def __init__(self, mesh, step, external_potential):
        '''
        mesh: numpy array of pts on which things are calculated
        rho: numpy array of the density calculated on the pts of the mesh
        external_potential: array with the external potential calculated on the pts of the mesh
        '''
        self.mesh = mesh
        self.step = step# only needed in the hartree_potential() which could easily be generalized for arbitrarly spaced meshes
        self.external_potential = external_potential

    def lda_correlation_potential(self, rho):#to optimize
        correlation_potential = numpy.empty(len(self.mesh))
        for i in range(len(rho)):
#            correlation_potential[i] = - 0.44 / (4 * math.pi) * 1 / math.pow(7.8 + math.pow(3 / (4 * math.pi * rho[i]), 1 / 3), 2) * 1 / math.pow(3 / (4 * math.pi * rho[i]), 2 / 3) * 1 / rho[i] - 0.44 / (7.8 + math.pow(3 / (4 * math.pi * rho[i]), 1 / 3))
            correlation_potential[i] = - 0.48 * 4 / 3 * math.pow(rho[i], 1 / 3)
        return correlation_potential
    def lda_correlation_density_rho_derivative(self, rho):#not needen in the kohn sham equation
        correlation_density_derivative = numpy.empty(len(self.mesh))
        for i in range(len(rho)):
#            correlation_density_derivative[i] = - 0.44 / (4 * math.pi) * 1 / math.pow(7.8 + math.pow(3 / (4 * math.pi * rho[i]), 1 / 3), 2) * 1 / math.pow(3 / (4 * math.pi * rho[i]), 2 / 3) * 1 / rho[i]
            correlation_density_derivative[i] = - 0.48 / 3 * math.pow(rho[i], 1 / 3)
        return correlation_density_derivative

    def lda_exchange_potential(self, rho):
        exchange_potential = numpy.empty(len(self.mesh))
        for i in range(len(self.mesh)):
            exchange_potential[i] = - math.pow(3 * rho[i] / math.pi, 1 / 3)
        return exchange_potential
    def lda_exchange_density_rho_derivative(self, rho):#not needen in the kohn sham equation
        exchange_density_derivative = numpy.empty(len(self.mesh))
        for i in range(len(self.mesh)):
            exchange_density_derivative[i] = - 3 / 4 * math.pow(3 / math.pi, 1 / 3) * math.pow(rho[i], - 2 / 3)
        return exchange_density_derivative

    def centrifugal_potential(self, l):
        return 0.5 * l * (l + 1) / (self.mesh * self.mesh)

    def hartree_potential(self, rho):
        hartree_potential = numpy.zeros(len(self.mesh))
        for i in range(len(self.mesh)):
            for j in range(len(self.mesh)):
                hartree_potential[i] = hartree_potential[i] + rho[j] * (abs(self.mesh[i] - self.mesh[j]) - (self.mesh[i] + self.mesh[j])) * self.mesh[j] / self.mesh[i] * self.step
        return - 2 * math.pi * hartree_potential

    def kohn_sham_potential(self, rho, l):
        return self.hartree_potential(rho) + self.lda_exchange_potential(rho) + self.lda_correlation_potential(rho) + self.external_potential + self.centrifugal_potential(l)
