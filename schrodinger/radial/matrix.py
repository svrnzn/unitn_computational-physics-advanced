import numpy

R = 7
Dx = 0.01
N = int(R / Dx) + 1
Dx_inverse_square = 1 / Dx**2

scattering_length = 1

def effective_potential(x, previous_eigenvector = numpy.zeros(N), l = 0):
    '''
    THE PART THAT FOLLOWS THE else IS WRONG
    '''
    if x != 0:
        return 0.5 * x**2 + l * (l + 1) / x**2 + scattering_length * (previous_eigenvector[round(x / Dx)])**2
    else:#missing the potential due to l
        return scattering_length * (previous_eigenvector[1])**2# DIRTY WORKAROUND

def radial_r(previous_eigenvector = numpy.zeros(N), l = 0, pt_derivative = 3):
    M = numpy.zeros( (N, N) )
    '''
    if pt_derivative = 3:
    '''
    for i in range(N):
        if i > 0:
            M[i][i - 1] = 0 - 0.5 * Dx_inverse_square
        M[i][i] = Dx_inverse_square + effective_potential(i * Dx, previous_eigenvector, l)
        if i < N - 1:
            M[i][i + 1] = 0 - 0.5 * Dx_inverse_square
    '''
    elif pt_derivative = 5:
        for i in range(N):
            if i > 1:
                M[i][i - 2] = 
            if i > 0:
                M[i][i - 1] = 
            M[i][i] = 
            if i < N - 1:
                M[i][i + 1] = 
            if i < N - 2:
                M[i][i + 2] = 
    '''
    eigenvalues, eigenvectors  = numpy.linalg.eigh(M)
    idx = eigenvalues.argsort()
    eigenvalues = eigenvalues[idx]
    eigenvectors = eigenvectors[:, idx]
#    for i in range(N):
#        eigenvectors[:, i] = eigenvectors[:, i] / (Dx * numpy.sum(eigenvectors[:, i]))
    return eigenvalues, eigenvectors


if __name__ == "__main__":
    print('please use one of the scripts to run this code!')
