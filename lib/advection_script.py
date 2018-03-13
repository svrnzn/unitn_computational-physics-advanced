import numpy
import math
from pde import advection

results = 'solution'
method = 'lax_wendroff'
initial_function = 'gaussian'
currant_factor = 0.5
t_max = 20

speed = 1

left_boundary = 0
right_boundary = 10
J = 101
x_step = (right_boundary - left_boundary) / (J - 1)
t_step = currant_factor * x_step / speed


if initial_function == 'box':
    def initial_f(x):# box
        if x <= 2.5 or x >= 7.5:
            return 0
        else:
            return 1
elif initial_function == 'gaussian':
    def initial_f(x):# gaussian
        x_0 = 5
        return math.exp(- (x - x_0)**2)

mesh = numpy.linspace(left_boundary, right_boundary, J)

current_solution = numpy.empty(J)
for i in range(J):
    current_solution[i] = initial_f(mesh[i])

if method == 'leapfrog':
    previous_solution = numpy.empty(J)
    for i in range(J):
        previous_solution = initial_f(mesh[i] + speed * t_step)

t = 0
while t < t_max + 2 * t_step:#why this?
    if method == 'lax_wendroff':
        current_solution = advection.lax_wendroff(current_solution, speed, x_step, t_step)
    elif method == 'lax_friedrichs':
        current_solution = advection.lax_friedrichs(current_solution, speed, x_step, t_step)
    elif method == 'leapfrog':
        previous_solution, current_solution = current_solution, advection.leapfrog(current_solution, previous_solution, speed, x_step, t_step)
    elif method == 'ftcs':
        current_solution = advection.ftcs(current_solution, speed, x_step, t_step)
    else:
        print('MACHEOOOHHHH!!!')
        quit
    if results == 'norm':
        print(t, numpy.linalg.norm(current_solution) / math.pow(J, 0.5))
    t = t + t_step

if results == 'solution':
    for i in range(J):
        print(mesh[i], current_solution[i])
