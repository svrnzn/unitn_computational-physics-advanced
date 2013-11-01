import numpy

def solve(mesh, step, potential, boundary_conditions = 'zero'):
    step_inverse_square = 1 / step**2
    no_of_points = len(mesh)
    M = numpy.zeros( (no_of_points, no_of_points) )
    if boundary_conditions == 'zero':
        for i in range(no_of_points):
            if i > 0:
                M[i][i - 1] = 0 - 0.5 * step_inverse_square
            M[i][i] = step_inverse_square + potential(mesh[i])
            if i < no_of_points - 1:
                M[i][i + 1] = 0 - 0.5 * step_inverse_square
    elif boudary_conditions == 'periodic':
        print('periodic boundary conditions are not implemented yet')
        exit
    eigenvalues, eigenvectors  = numpy.linalg.eigh(M)
    idx = eigenvalues.argsort()
    eigenvalues = eigenvalues[idx]
    eigenvectors = eigenvectors[:, idx]
    eigenvestors = eigenvectors / step
    return eigenvalues, eigenvectors

if __name__ == "__main__":
    step = 0.005
    R = 5
    mesh = numpy.arange(0, R, step)
    def potential(x):
        return 0.5 * (x)**2
    eigenvalues, eigenvectors = solve(mesh, step, potential)
#    for eigenvalue in eigenvalues:
#        print(eigenvalue)
    for i in eigenvectors[:, 0]:
        print(i)
