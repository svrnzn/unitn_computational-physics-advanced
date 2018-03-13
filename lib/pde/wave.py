import numpy

'''
YOU SHOULD IMPLEMENT THE ADVECTION ALGORITHM BETTER AND IMPORT THAT ONE HERE!
'''

def lax_friedrichs(current_solution, speed, x_step, t_step):
    current_r = 
    current_s = 
    '''
    r = 0.5 * (numpy.roll(current_r, - 1) + numpy.roll(current_r, 1)) + 0.5 * t_step / x_step * (numpy.roll(current_s, - 1) - numpy.roll(current_s, 1))
    s = 0.5 * (numpy.roll(current_s, - 1) + numpy.roll(current_s, 1)) + 0.5 * t_step / x_step * (numpy.roll(current_r, - 1) - numpy.roll(current_r, 1))
    '''
    return 

def lax_wendroff():
    pass

def leapfrog(current_solution, previous_solution, speed, x_step, t_step):
    return (speed * t_step / x_step)**2 * numpy.roll(current_solution, - 1) + 2 * (1 - (speed * t_step / x_step)**2) * current_solution + (speed * t_step / x_step)**2 * numpy.roll(current_solution, 1) - previous_solution
