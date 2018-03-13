import numpy
import math

class JohnCena:
    def __init__(self, mesh, step):
        self.mesh = mesh
        self.step = step
    def __call__(self, potential):
        """Solve Scrodinger equation

        Return tuple of numpy arrays (eigenvalues, eigenvectors).
        Eigenvectors are columns of the 2D eigenvectors array.
        Eigenvectors and eigenvalues are ordered with increasing energy.
        """

        step_inverse_square = 1 / self.step**2
        no_of_points = len(self.mesh)
        M = numpy.zeros( (no_of_points, no_of_points) )

        for i in range(no_of_points):
            if i > 0:
                M[i][i-1] = - step_inverse_square/2
            M[i][i] = step_inverse_square + potential[i]
            if i < no_of_points-1:
                M[i][i+1] = - step_inverse_square/2
        eigenvalues, eigenvectors  = numpy.linalg.eigh(M)
        idx = eigenvalues.argsort()
        eigenvalues = eigenvalues[idx]
        eigenvectors = eigenvectors[:, idx]
        eigenvectors = eigenvectors/math.sqrt(self.step)

        return eigenvalues, eigenvectors
