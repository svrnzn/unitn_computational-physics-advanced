import sys
import numpy
import matrix

if sys.argv[1] == '--auto':
    epsilon = float(sys.argv[2])
    difference = epsilon + 1
    ITERATIONS = 0
    eigenvectors = numpy.zeros( (matrix.N, matrix.N) )
    eigenvalues = numpy.zeros(matrix.N)
    while difference > epsilon:
        old_eigenvalue = eigenvalues[0]
        eigenvalues, eigenvectors = matrix.radial_r(eigenvectors[:, 0])
        for i in range(1, matrix.N):
            eigenvectors[i] = eigenvectors[i] / (i * matrix.Dx)
        eigenvectors[0] = eigenvectors[1]#DIRT
        for i in range(matrix.N):
            eigenvectors[:, i] = eigenvectors[:,i] / (numpy.sum(eigenvectors[:,i] * matrix.Dx))
        difference = abs(old_eigenvalue - eigenvalues[0])
        ITERATIONS = ITERATIONS + 1
        print(eigenvalues[0])
    print('Stopped in {} iterations.'.format(ITERATIONS))

else:
    ITERATIONS = int(sys.argv[1])
    eigenvectors = numpy.zeros( (matrix.N, matrix.N) )
    for i in range(1, matrix.N):
        eigenvectors[i] = eigenvectors[i] / (i * matrix.Dx)
    eigenvectors[0] = eigenvectors[1]#DIRT
    for i in range(matrix.N):
        eigenvectors[:, i] = eigenvectors[:,i] / (numpy.sum(eigenvectors[:,i] * matrix.Dx))
    for i in range(ITERATIONS):
        eigenvalues, eigenvectors = matrix.radial_r(eigenvectors[:, 0])
        print(eigenvalues[0])
