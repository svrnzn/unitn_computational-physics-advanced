import numpy

'''
ERROR MANAGEMENT TO BE IMPLEMENTED SOME PROPER WAY
'''

def ftcs(current_solution, speed, x_step, t_step, boundary_conditions = 'periodic'):
    if boundary_conditions = 'periodic':
        return current_solution - 0.5 * speed * t_step / x_step * (numpy.roll(current_solution, - 1) - numpy.roll(current_solution, 1))
    else:
        print('boudary conditions not implemented')
        exit

def lax_friedrichs(current_solution, speed, x_step, t_step, boundary_conditions = 'periodic'):
    if boundary_conditions = 'periodic':
        return 0.5 * (numpy.roll(current_solution, 1) + numpy.roll(current_solution, - 1)) - 0.5 * speed * t_step / x_step * (numpy.roll(current_solution, - 1) - numpy.roll(current_solution, 1))
    else:
        print('boudary conditions not implemented')
        exit

def leapfrog(current_solution, previous_solution, speed, x_step, t_step, boundary_conditions = 'periodic'):
    if boundary_conditions = 'periodic':
        return previous_solution - speed * t_step / x_step * (numpy.roll(current_solution, - 1) - numpy.roll(current_solution, 1))
    else:
        print('boudary conditions not implemented')
        exit

def lax_wendroff(current_solution, speed, x_step, t_step, boundary_conditions = 'periodic'):
    if boundary_conditions = 'periodic':
        return current_solution - 0.5 * speed * t_step / x_step * (numpy.roll(current_solution, - 1) - numpy.roll(current_solution, 1)) + 0.5 * (speed * t_step / x_step)**2 * (numpy.roll(current_solution, - 1) - 2 * current_solution  + numpy.roll(current_solution, 1))
    else:
        print('boudary conditions not implemented')
        exit
