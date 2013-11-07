import numpy
import math

def solve(mesh, step, potential, boundary_conditions = 'zero'):
    """Returns a tuple with eigenvalues and eigenvectors (eigenvectros are the columns of the 2D numpy array)"""
    step_inverse_square = 1 / step**2
    no_of_points = len(mesh)
    M = numpy.zeros( (no_of_points, no_of_points) )
    if boundary_conditions == 'zero':
        for i in range(no_of_points):
            if i > 0:
                M[i][i - 1] = - 0.5 * step_inverse_square
            M[i][i] = step_inverse_square + potential(mesh[i])
            if i < no_of_points - 1:
                M[i][i + 1] = - 0.5 * step_inverse_square
    elif boudary_conditions == 'periodic':
        print('periodic boundary conditions are not implemented yet')
        exit
    eigenvalues, eigenvectors  = numpy.linalg.eigh(M)
    idx = eigenvalues.argsort()
    eigenvalues = eigenvalues[idx]
    eigenvectors = eigenvectors[:, idx]
    eigenvectors = eigenvectors / math.sqrt(step)
    return eigenvalues, eigenvectors
