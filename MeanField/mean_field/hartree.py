import numpy as np


class EffectivePotential:
    """Compute Hartree effective potential from density."""

    def __init__(self, grid):
        """
        Arguments
        grid : numpy array of evenly spaced points.
        """
        self.grid = grid
        self.step = np.abs(grid[1]-grid[0])

    def __call__(self, density):
        """Return Hartree effective potential.

        The interaction used here is an approximation to the Coulomb
        potential to avoid divercences of the mean field theory.

        Argument
        density : density computed at point in self.grid.

        Returns
        hartree_potential : Hartree effective potential computed at
            points in self.grid .
        """
        hartree_potential = np.empty(self.grid.size)
        for i in range(self.grid.size):
            numerator = (density
                         *(np.abs(self.grid[i]-self.grid)
                            -(self.grid[i]+self.grid))*self.grid)
            denominator = self.grid[i]
            hartree_potential[i] = self.step*np.sum(numerator)/denominator
        # Note that a factor 2 is already in hartree_potential.
        return - 2*np.pi*hartree_potential
