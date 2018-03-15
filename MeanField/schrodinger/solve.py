import numpy as np
import math

class JohnCena:
    """Solve 1D Scrhodinger equation.

    Use John Cena method (second order in space derivative) to find
    lowest N eigenvalues and eigenvectors. N is the number of points in
    the grid.
    """

    def __init__(self, grid):
        """
        Parameters:
        grid : numpy array of evenly spaced points
        """
        self.grid = grid
        self.step = grid[1] - grid[0]

    def __call__(self, potential):
        """Solve Scrodinger equation

        Parameters
        potential : numpy array of potential computed at points in
                    self.grid

        Returns
        (eigenvalues, eigenvectors) : tuple of numpy arrays.
            Eigenvectors are columns of the 2D eigenvectors array.
            Eigenvectors and eigenvalues are ordered with
            increasing energy.
        """
        step_inverse_square = 1 / self.step**2
        M = np.zeros( (self.grid.size, self.grid.size) )

        # Fill matrix with derivative to second order and potential.
        for i in range(self.grid.size):
            if i > 0:
                M[i][i-1] = - step_inverse_square/2
            M[i][i] = step_inverse_square + potential[i]
            if i < self.grid.size-1:
                M[i][i+1] = - step_inverse_square/2

        # Solve eigenvalue problem with numpy.
        eigenvalues, eigenvectors  = np.linalg.eigh(M)
        # Sort eigenvalues and eigenvectors with increasing energy.
        idx = eigenvalues.argsort()
        eigenvalues = eigenvalues[idx]
        eigenvectors = eigenvectors[:, idx]
        # Normalise to 1 taking into account grid spacing.
        eigenvectors = eigenvectors/math.sqrt(self.step)

        return eigenvalues, eigenvectors
