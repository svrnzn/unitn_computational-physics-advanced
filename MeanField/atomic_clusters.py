import numpy as np
import matplotlib.pyplot as plt
import schrodinger.solve
from mean_field import kohn_sham
from collections import OrderedDict# Only needed to remove multiple
                                   # labels in plot's legend.


# ========
# SETTINGS
# ========

# grid settings.
R = 20
step = 0.04

# Number of electrons and shell structure.
N = 40
# Occupied orbitals. shell_structure[i] is the number of lowest
# orbitals with angular momentum quantum number L = i to be filled.
#shell_structure = np.array([1])# 2 atoms
#shell_structure = np.array([1, 1])# 8 atoms
#shell_structure = np.array([2, 1, 1])# 20 atoms
shell_structure = np.array([2, 2, 1, 1])# 40 atoms
#shell_structure = np.array([2, 2, 1, 1, 1])# 58 atoms

# Ion's core properties.
# Radious of a single core ion in units Bohr radious.
wigner_seitz_radius = 3.93
# Initial guess density settings. Disatnce over which the initial
# density drops in units Bohr radious.
healing_lenght = 1

# Self consistency and stopping criteria.
# Self consistency update factor. The density at a successive step is
# the weighted average of the current and previous densities with
# weights beta and (1-beta) respectively.
beta = 0.25
# Relative energy change between successive iterations of the self
# consistent procedure.
energy_rel_delta = 0.001


# ===========
# DEFINITIONS
# ===========

# Compute core radious and density from Wiegner-Seitz radious.
core_radious = np.power(N, 1/3) * wigner_seitz_radius
core_density = 3*N/4/np.pi/np.power(core_radious, 3)

# External potential.
def uniform_sphere_potential(grid, core_radious, core_density):
    """Return uniform sphere electrostatic potential.

    Electrostatic potential sourced by a uniformly (positively) charged
    sphere of radious core_radious and density core_density. Length in
    unit Bohr radious, charge in unit e=1.

    Arguments
    grid : numpy array of points where potential is computed
    core_radious : radious of charged sphere
    core_density : charge density of sphere.

    Returns
    external_potential : numpy array of electrostatic potential
        computed at points in grid.
    """
    external_potential = np.empty(grid.size)
    for i in range(len(grid)):
        if grid[i] <= core_radious:
            external_potential[i] = -(2*np.pi*core_density
                                      *(np.power(core_radious, 2)
                                      - np.power(grid[i], 2)/3))
        else:
            external_potential[i] = -(4*np.pi*core_density
                                      *np.power(core_radious, 3)/3/grid[i])
    return external_potential

# Initial guess density.
def initial_guess_density(grid, core_radious, healing_lenght, normalisation):
    """Return initial guess density.

    Fermi-Dirac like initial guess density suitable for atomic cluster
    systems as in Modern Many-Particle Physics, Enrico Lipparini,
    chapter 2. Arbitrary units.

    Arguments
    grid : numpy array of points where density is computed
    core_radious : radious of positively charged ions core
    healing_lenght : lenght over which the initial density drops
    normalisation : normalisation of initial guess densityself.

    Returns
    density : numpy array of density computed at points in grid, note
        that it is normalised radially i.e. to normalisation/4/np.pi
        when integrating along the grid.
    """
    step = np.abs(grid[1]-grid[0])
    density = np.empty(grid.size)
    for i in range(grid.size):
        density[i] = 1 / (1 + np.exp(-(core_radious-grid[i])/healing_lenght))
    density = (normalisation*density
               /(4*np.pi*np.sum(density*grid*grid)*step))
    return density

# Set up grid.
grid = np.arange(step, R, step)# starts at step to avoid zero division
# Compute external potential.
external_potential = uniform_sphere_potential(grid, core_radious, core_density)

# Instantiate classes to compute Kohn-Sham effective potential, solve
# Schrodinger equation, Kohn-Sham energy and density.
effective_potential_compute = kohn_sham.LDAEffectivePotential(grid,
                                                              external_potential)
schrodinger_solve = schrodinger.solve.JohnCena(grid)
energy_compute = kohn_sham.LDAEnergy(grid, shell_structure)
density_compute = kohn_sham.Density(grid, shell_structure)

# Compute initial guess density.
density = initial_guess_density(grid, core_radious, healing_lenght, N)


# ================
# SELF CONSISTENCY
# ================

# Self consistent cycle.
j = 0 # Needed to enter cycle at least twice.
while j < 2 or np.abs(old_energy-energy)/np.abs(energy) > energy_rel_delta:
    # Update density for all iterations after the first.
    try:
        density = (1-beta)*old_density + beta*density_compute(eigenvectors)
    except:
        pass

    # Define variables to store eigenstars.
    eigenvalues = np.empty(0)
    eigenvectors = np.empty( (len(grid), 0) )

    # Iterate over angular momenum quantum numbers involved.
    for L in range(shell_structure.size):
        # Compute effective potential for given density and angular
        # momentum quantum number.
        effective_potential = effective_potential_compute(density, L)
        # Solve Schrodinger euqation.
        tmp_eigenstars = schrodinger_solve(effective_potential)

        # Separate eigenvalues from eigenvectors.
        tmp_eigenvalues, tmp_eigenvectors = tmp_eigenstars
        # Keep only needed eigenstars as defined in shell_structure.
        tmp_eigenvalues = tmp_eigenvalues[0:shell_structure[L]]
        tmp_eigenvectors = tmp_eigenvectors[:, 0:shell_structure[L]]
        # Convert from reduced to radial wavefunction.
        for i in range(shell_structure[L]):
            tmp_eigenvectors[:, i] = tmp_eigenvectors[:, i]/grid

        # Append eigenstars for a given L to previously computed ones.
        eigenvalues = np.append(eigenvalues, tmp_eigenvalues)
        eigenvectors = np.append(eigenvectors, tmp_eigenvectors, 1)

    # Store previous density to use in next iteration and compute new
    # one.
    old_density = density
    density = density_compute(eigenvectors)

    # Store previous energy, after the first iteration, for stopping
    # criterion and compute new one.
    try:
        old_energy = energy
    except:
        pass
    energy = energy_compute(eigenvalues, density)

    # Feedback in stdout.
    print('Energy per particle (Hartree):', energy/N)

    j = j + 1


# ======
# OUTPUT
# ======

# Compute effective potentials for the various Ls.
potentials = np.empty( (shell_structure.size, grid.size) )
for L in range(shell_structure.size):
    potentials[L] = effective_potential_compute(density, L)

# Write file in /tmp/ with data.
output = open('/tmp/' + str(N) + 'at' + '_epsilon' + str(energy_rel_delta)
              + '_beta' + str(beta) + '.dat', 'w')
output.write('# energy [Hartree]:' + str(energy) + '\n')
output.write('# iterations before stopping: ' + str(j) + '\n')
output.write('# Runtime settings\n')
output.write('# number of atoms: ' + str(N) + '\n')
output.write('# shell structure: ' + str(shell_structure) + '\n')
output.write('# wiegner seitz radious: ' + str(wigner_seitz_radius) + '\n')
output.write('# R_c: ' + str(core_radious) + '\n')
output.write('# beta: ' + str(beta) + '\n\n')

title_line = '# r\tdensity\t'
for L in range(shell_structure.size):
    title_line = title_line + '\t' + 'potential L=' + str(L)
title_line = title_line + '\n'
output.write(title_line)

for i in range(grid.size):
    data_line = str(grid[i]) + '\t' + str(density[i])
    for L in range(shell_structure.size):
        data_line = data_line + '\t' + str(potentials[L, i])
    data_line = data_line + '\n'
    output.write(data_line)


# Plot straight to display.
# Set figure size.
plt.figure(figsize=(6,12))

# Define list to assign colors.
colors = ['g', 'r', 'c', 'm', 'y', 'k', 'w']
# Define angular momenta associated to the shell structure to assign
# labels.
angular_momenta = np.empty(0)
for L in range(shell_structure.size):
    angular_momenta = np.append(angular_momenta,
                                    L*np.ones(shell_structure[L]))

# First subplot: Density.
plt.subplot(3, 1, 1)
plt.plot(grid, density, lw=1, c='b', label='Kohn-Sham density')
plt.plot(grid, initial_guess_density(grid, core_radious, healing_lenght, N),
         lw=0.75, ls='--', c='b', label='Initial guess density')
plt.legend()
plt.ylabel('Density')

# Second subplot: Orbitals.
plt.subplot(3, 1, 2)#, sharex=plt.subplot(3, 1, 1))
for i in range(eigenvectors.shape[1]):
    plt.plot(grid, grid*np.power(eigenvectors[:, i], 2), lw=1,
             c=colors[int(angular_momenta[i])],
             label='L = '+str(int(angular_momenta[i])))
# Remove multiple labels.
handles, labels = plt.gca().get_legend_handles_labels()
by_label = OrderedDict(zip(labels, handles))
plt.legend(by_label.values(), by_label.keys())
plt.ylabel('Orbitals Density')

# Third subplot: Effective Potentials.
plt.subplot(3, 1, 3)#, sharex=plt.subplot(3, 1, 1))
plt.ylim( (-0.25, 0.3) )
for L in range(shell_structure.size):
    plt.plot(grid, potentials[L], lw=1, c=colors[L], label='L = '+str(L))
plt.legend()
plt.ylabel('Effective Potential')
plt.xlabel(r'Radious [$a_0$]')

plt.show()
