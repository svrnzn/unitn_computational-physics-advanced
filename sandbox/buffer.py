R = 7
Dx = 0.01
N = int(R / Dx) + 1
Dx_inverse_square = 1 / Dx**2

scattering_length = 1

class EffectivePotential:
    def effective_potential(x, previous_eigenvector = numpy.zeros(N), l = 0):
        '''
        THE PART THAT FOLLOWS THE else IS WRONG
        '''
        if x != 0:
            return 0.5 * x**2 + l * (l + 1) / x**2 + scattering_length * (previous_eigenvector[round(x / Dx)])**2
        else:#missing the potential due to l
            return scattering_length * (previous_eigenvector[1])**2# DIRTY WORKAROUND
