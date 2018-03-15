import numpy as np
import mean_field.hartree as hartree


def degeneracies(shell_structure):
    """Compute degeneracies from shell structure.

    Argument
    shell_structure : numpy array
        shell_structure[i] : number of orbitals with L = i involved.

    Return
    degeneracies : numpy array
        degeneracy of states forming the shell structure ordered with
        increasing L and energy.
    """
    degeneracies = np.empty(0)
    for L in range(shell_structure.size):
        degeneracies = np.append(degeneracies,
                                    ((2 * (2*L + 1))
                                      * np.ones(shell_structure[L])))
    return degeneracies

def lda_x_energy_density(density):
    """Return LDA exchange energy density.

    See Modern Many-Particle Physics, Enrico Lipparini, chapter 2.
    """
    lda_x_energy_density = - 3*np.power(3*density/np.pi, 1/3)/4
    return lda_x_energy_density

def lda_x_energy_density_d(density):
    """Derivative with respect to density of lda_x_energy_density."""
    lda_x_energy_density_d = - np.power(3/np.pi, 1/3)/np.power(density, 2/3)/4
    return np.zeros(density.size)

def lda_c_energy_density(density):
    """Return LDA exchange correlation density.

    See Modern Many-Particle Physics, Enrico Lipparini, chapter 2self.
    """
    lda_c_energy_density = - 0.44 / (7.8+np.power(3/4/np.pi/density, 1/3))
    return lda_c_energy_density

def lda_c_energy_density_d(density):
    """Derivative with respect to density of lda_c_energy_density."""
    lda_c_energy_density_d = - (0.44*np.power(3/4/np.pi, 1/4)/3
                                /np.power(7.8+np.power(3/4/np.pi/density, 1/3),
                                          2)
                                /np.power(density, 2/3))
    return lda_c_energy_density_d

def centrifugal_potential(grid, L):
    """Return centrifugal potential.

    Arguments
    grid : numpy array of points where to compute the tpotential
    L : angular momentum quantum number.
    """
    return L*(L+1)/np.power(grid, 2)/2


class LDAEffectivePotential:

    def __init__(self, grid, external_potential):
        """
        Arguments
        grid : numpy array of evenly spaced points
        external_potential : numpy array of external potential
                             computed at points in grid.
        """
        self.grid = grid
        self.step = np.abs(grid[1] - grid[0])
        self.external_potential = external_potential
        self.hartree_potential = hartree.EffectivePotential(grid)
        self.lda_x_energy_density = lda_x_energy_density
        self.lda_c_energy_density = lda_c_energy_density
        self.lda_x_energy_density_d = lda_x_energy_density_d
        self.lda_c_energy_density_d = lda_c_energy_density_d
        self.centrifugal_potential = centrifugal_potential

    def __call__(self, density, L):
        """Return effective potential in LDA approximation.

        Arguments
        density : numpy array of density at points in self.grid
        L : orbital angular momentum quantum number.

        Returns
        : numpy array of effective potential computed at points in
          self.grid.
        """
        centrifugal_potential = self.centrifugal_potential(self.grid, L)
        hartree_potential = self.hartree_potential(density)
        lda_x_energy_density = self.lda_x_energy_density(density)
        lda_c_energy_density = self.lda_c_energy_density(density)
        lda_x_energy_density_d = self.lda_x_energy_density_d(density)
        lda_c_energy_density_d = self.lda_c_energy_density_d(density)
        return (centrifugal_potential + self.external_potential
                + hartree_potential
                + lda_x_energy_density + lda_c_energy_density
                + lda_x_energy_density_d*density
                + lda_c_energy_density_d*density)


class LDAEnergy:
    """Compute Kohn-Sham energy in LDA."""

    def __init__(self, grid, shell_structure):
        """
        Arguments
        grid : numpy array of evenly spaces points
        shell_structure : numpy array
            shell_structure[i] : number of orbitals with L = i involved.
        """
        self.grid = grid
        self.step = np.abs(grid[1]-grid[0])
        self.shell_structure = shell_structure
        self.degeneracies = degeneracies(shell_structure)
        self.hartree_potential = hartree.EffectivePotential(grid)
        self.lda_x_energy_density_d = lda_x_energy_density_d
        self.lda_c_energy_density_d = lda_c_energy_density_d

    def __call__(self, occupied_eigenvalues, density):
        """Return Kohn-Sham energy in LDA.

        Arguments
        occupied_eigenvalues : numpy array of eigenvalues of occupied
            states
        density : numpy array of density computed at points in
            self.grid.

        Returns
        : numpy array of Kohn-Sham effective potential in LDA
            computed at points in self.grid.
        """
        eigenvalues_contrib = np.sum(occupied_eigenvalues*self.degeneracies)
        hartree_potential = self.hartree_potential(density)
        hartree_contrib = -(2*np.pi*self.step
                             *np.sum(self.grid*self.grid*density
                                     *hartree_potential))
        x_contrib_density = - self.lda_x_energy_density_d(density)*density
        c_contrib_density = - self.lda_c_energy_density_d(density)*density
        xc_contrib = (4*np.pi*self.step
                      * np.sum(self.grid*self.grid*density
                               * (x_contrib_density+c_contrib_density)))
        return eigenvalues_contrib + hartree_contrib + xc_contrib


class Density:
    """Compute density from shell structure."""

    def __init__(self, grid, shell_structure):
        """
        Arguments
        grid : numpy array of evenly spaces points
        shell_structure : numpy array
            shell_structure[i] : number of orbitals with L = i involved.
        """
        self.grid = grid
        self.step = np.abs(grid[1] - grid[0])
        self.shell_structure = shell_structure
        self.degeneracies = degeneracies(shell_structure)

    def __call__(self, occupied_states):
        """Return density computed from shell structure.

        Argument
        occupied_states : numpy 2D array, columns are eigenvectors
            they must be ordered with increasing L.

        Returns
        : radial density.
        """
        density = np.zeros(self.grid.size)
        for i in range(occupied_states.shape[1]):
            density = (density
                       + self.degeneracies[i]*np.power(occupied_states[:, i],
                                                       2))
        # Divide by 4*pi to normalise radially.
        return density/4/np.pi
