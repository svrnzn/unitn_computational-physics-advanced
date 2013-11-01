import numpy
import math

'''
 - INTERACTION POTENTIAL THING!
 - the way the centrifugal potential is computed is confusing
 - some potential() and some potential.compute()
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
        '''
        Returns the Hartree potential integrated some rough way
        '''
        centrifugal_potential = Centrifugal(self.l)
        hartree_potential = 0
        for y in self.mesh:
            if y <= x:
                hartree_potential = hartree_potential + self.rho(y) * math.pow(y, 2) * self.step / x
            else:
                hartree_potential = hartree_potential + self.rho(y) * y * self.step
        return 4 * math.pi * hartree_potential + centrifugal_potential.compute(x)

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
        return hartree_potential.compute(x) + self.exchange_potential(self.rho(x)) + self.correlation_potential(self.rho(x)) + self.external_potential(x)
